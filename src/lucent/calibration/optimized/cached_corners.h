/* Cache finite corners within one probability-grid invocation. Preserve endpoint
   order because swapping the two Owen subtractions changes rounding. */
#include "corners.h"
#include <string.h>

typedef struct {
    uint64_t h, k, r;
    double value;
    int used;
} Entry;

static uint64_t bits(double x) {
    uint64_t b;
    memcpy(&b,&x,sizeof b);
    return b;
}

static double reuse(Entry *cache, double h, double k, double ph, double pk, double r, double sd) {
    if (!isfinite(h) || !isfinite(k)) return corner(h,k,ph,pk,r,sd);
    uint64_t x=bits(h), y=bits(k), z=bits(r);
    uint64_t hash=x*0x9e3779b97f4a7c15ULL+y*0xbf58476d1ce4e5b9ULL+z*0x94d049bb133111ebULL;
    Entry *e=&cache[hash >> 54];
    if (e->used && e->h==x && e->k==y && e->r==z) return e->value;
    double value;
    if (h==k && h!=0.) {
        double term=owen(h,(1.-r)/sd,1);
        value=.5*(ph+pk)-term-term;
    } else {
        value=corner(h,k,ph,pk,r,sd);
    }
    *e=(Entry){x,y,z,value,1};
    return value;
}

void cached_rectangles(int rows,int n,int threads,const double *z,const double *cdf,
                       const double *rho,double *result) {
    #pragma omp parallel if(threads>1) num_threads(threads)
    {
        Entry cache[1024]={0};
        #pragma omp for schedule(static)
        for(int i=0;i<rows*n;++i) {
            int j=i%n;
            const double *a=z+j,*p=cdf+j;
            double r=rho[i],sd=sqrt(1.-r*r);
            result[i]=reuse(cache,a[n],a[3*n],p[n],p[3*n],r,sd)
                     -reuse(cache,a[0],a[3*n],p[0],p[3*n],r,sd)
                     -reuse(cache,a[n],a[2*n],p[n],p[2*n],r,sd)
                     +reuse(cache,a[0],a[2*n],p[0],p[2*n],r,sd);
        }
    }
}
