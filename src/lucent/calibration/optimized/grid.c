/* Bounded Hermite expansions on the unchanged Eq. 12 probability grid. */
#include "cached_corners.h"
#include <stdlib.h>

enum { MAX_ORDER=96, CACHE_SIZE=512 };
typedef struct {
    uint64_t key[4];
    int index,used;
} CategoryKey;
typedef struct {
    double base,norm,term[MAX_ORDER];
} Polynomial;

static void basis(double x,const double *roots,double *v) {
    if(!isfinite(x)) { memset(v,0,MAX_ORDER*sizeof(*v));return; }
    v[0]=exp(-.5*x*x)/sqrt(2.*acos(-1.));
    v[1]=x*v[0];
    for(int k=2;k<MAX_ORDER;++k) v[k]=(x*v[k-1]-roots[k-1]*v[k-2])/roots[k];
}

/* Knock out eight terms per dependent step. */
static double evaluate(const double *term,int k,double r,double r8) {
    double lane[8];
    memcpy(lane,term+k-8,sizeof lane);
    for(int m=k-16;m>=0;m-=8) {
        for(int j=0;j<8;++j) lane[j]=term[m+j]+r8*lane[j];
    }
    double out=lane[7];
    for(int j=6;j>=0;--j) out=lane[j]+r*out;
    return out;
}

void adaptive_rectangles(int rows,int n,int threads,const double *z,const double *cdf,
                         const double *rho,double *result) {
    CategoryKey seen[CACHE_SIZE]={0};
    Polynomial *poly=malloc((size_t)n*sizeof(*poly));
    int *index=malloc((size_t)n*sizeof(*index));
    if(!poly || !index) {
        free(poly);free(index);
        cached_rectangles(rows,n,threads,z,cdf,rho,result);return;
    }
    double roots[MAX_ORDER];
    for(int k=0;k<MAX_ORDER;++k) roots[k]=sqrt((double)k);
    int count=0;
    for(int j=0;j<n;++j) {
        uint64_t key[4],hash=0;
        for(int k=0;k<4;++k) {
            key[k]=bits(z[k*n+j]);
            hash=(hash^(key[k]>>24)^key[k])*0x9e3779b97f4a7c15ULL;
        }
        CategoryKey *entry=&seen[hash>>55];
        if(entry->used && memcmp(key,entry->key,sizeof key)==0) {index[j]=entry->index;continue;}
        memcpy(entry->key,key,sizeof key);
        index[j]=entry->index=count++; entry->used=1;
        Polynomial *p=&poly[index[j]];
        double v[4][MAX_ORDER];
        for(int k=0;k<4;++k) basis(z[k*n+j],roots,v[k]);
        p->base=(cdf[n+j]-cdf[j])*(cdf[3*n+j]-cdf[2*n+j]);
        p->norm=sqrt(p->base);
        for(int k=0;k<MAX_ORDER;++k)
            p->term[k]=(v[0][k]-v[1][k])*(v[2][k]-v[3][k])/(k+1.);
    }
    #pragma omp parallel if(threads>1) num_threads(threads)
    {
        Entry cache[1024]={0};
        #pragma omp for schedule(static)
        for(int i=0;i<rows*n;++i) {
            int j=i%n;
            double r=rho[i],a_r=fabs(r);
            if(a_r<.76) {
                const Polynomial *p=&poly[index[j]];
                double r2=r*r,r4=r2*r2,r8=r4*r4;
                double power=r8,limit=5.e-15*p->norm;
                int k=8;
                while(k<64 && power*a_r>limit) {power*=power;k*=2;}
                if(k==64 && power*a_r>limit) {
                    double r16=r8*r8;
                    power*=r16*r16; k=96;
                }
                if(power*a_r<=limit) {
                    double est=p->base+r*evaluate(p->term,k,r,r8);
                    /* Parseval bound, checked against the computed mass. */
                    if(est>1.01e-10 && power*a_r*p->norm<=1.e-14*est) {
#ifdef REPORT_DEGREE
                        result[i]=k;
#else
                        result[i]=est;
#endif
                        continue;
                    }
                }
            }
#ifdef REPORT_DEGREE
            result[i]=0.;
#else
            const double *a=z+j,*p=cdf+j;
            double sd=sqrt(1.-r*r);
            result[i]=reuse(cache,a[n],a[3*n],p[n],p[3*n],r,sd)
                     -reuse(cache,a[0],a[3*n],p[0],p[3*n],r,sd)
                     -reuse(cache,a[n],a[2*n],p[n],p[2*n],r,sd)
                     +reuse(cache,a[0],a[2*n],p[0],p[2*n],r,sd);
#endif
        }
    }
    free(index);free(poly);
}
