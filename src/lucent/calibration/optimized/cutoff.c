/* C11, binary64. Reference calls come from the same SciPy build as the CAD sampler. */
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

typedef double (*special_fn)(double, int);
static special_fn normal_cdf, normal_quantile;

void bind_normal(special_fn cdf, special_fn quantile) {
    normal_cdf = cdf;
    normal_quantile = quantile;
}

/* Use the original rounded CDF inside the guard interval, including its ends.
   This covers the plateau of ndtr near 0.5 as well as ordinary rounded ties. */
void prepare_guards(int n, const double *cum, double *lo, double *hi) {
    const double margin = 64.0 * DBL_EPSILON;
    for (int j=0; j<n; ++j) {
        double p = cum[j];
        lo[j] = -INFINITY;
        hi[j] = INFINITY;
        if (p > margin && p < 1.0-margin) {
            double a = normal_quantile(p-margin,1);
            double b = normal_quantile(p+margin,1);
            if (normal_cdf(a,1) < p && normal_cdf(b,1) > p) {
                lo[j] = a;
                hi[j] = b;
            }
        }
    }
}

static inline int below(double z, int j, const double *cum, const double *lo, const double *hi) {
    if (z < lo[j]) return 1;
    if (z > hi[j]) return 0;
    return normal_cdf(z,1) <= cum[j];
}

static inline int bin(double z, int ncut, const double *cum, const double *lo, const double *hi) {
    if (ncut == 0 || below(z,0,cum,lo,hi)) return 0;
    int l=1, r=ncut;
    while (l<r) {
        int mid=l+(r-l)/2;
        if (below(z,mid,cum,lo,hi)) r=mid;
        else l=mid+1;
    }
    return l;
}

void labels(int n, const double *z, int ncut, const double *cum, const double *lo,
            const double *hi, int64_t *out) {
    for (int i=0; i<n; ++i) out[i]=bin(z[i],ncut,cum,lo,hi);
}

/* Read one projected field, count its section-wide categories, then produce either
   a scalar contrast or the receiver-by-field matrix consumed by all contrasts. */
void transform(int n, int batches, int k, int ncut, int64_t si, int64_t sb,
               const double *z, const double *cum, const double *lo, const double *hi,
               const int64_t *mapping, const int64_t *outside, const double *contrast,
               int64_t *bins, int64_t *counts, double *ranks, double *out) {
    for (int b=0; b<batches; ++b) {
        for (int j=0; j<k; ++j) counts[j]=outside[j];
        for (int i=0; i<n; ++i) {
            bins[i]=mapping[bin(z[i*si+b*sb],ncut,cum,lo,hi)];
            ++counts[bins[i]];
        }
        int64_t total=0, prefix=0;
        for (int j=1; j<k; ++j) total+=counts[j];
        ranks[0]=0.0;
        for (int j=1; j<k; ++j) {
            prefix+=counts[j];
            ranks[j]=(prefix-0.5*counts[j]+0.5)/(total>0 ? total : 1);
        }
        if (contrast) {
            double s=0.0;
            for (int i=0; i<n; ++i) s+=contrast[i]*ranks[bins[i]];
            out[b]=s;
        } else {
            for (int i=0; i<n; ++i) out[(size_t)i*batches+b]=ranks[bins[i]];
        }
    }
}
