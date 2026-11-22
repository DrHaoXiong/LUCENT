/* Keep the row-bound ABI for paired comparison. Use the four original Owen
   corners when the coarse interval alone cannot rule out the extreme tail. */
#include "cached_corners.h"
#include <stdlib.h>

enum { DEGREE=16 };
typedef struct { double c[DEGREE],mass,variance,lo,hi; } Marginal;
typedef struct { double c[DEGREE],base,scale,upper; } Cell;

int row_bounds(int rows,int n,int categories,const double *z,const double *cdf,
               const int64_t *left,const int64_t *right,const int64_t *weights,
               const double *rho,double *low,double *high,const double *inner_l,
               const double *inner_u,int tighten) {
    (void)inner_l;(void)inner_u;(void)tighten;
    Marginal *bins=calloc(categories,sizeof(*bins));
    Cell *cells=malloc((size_t)n*sizeof(*cells));
    if(!bins||!cells) {free(bins);free(cells);return 0;}
    double last[DEGREE]={0},roots[DEGREE+1];
    for(int k=0;k<=DEGREE;++k) roots[k]=sqrt((double)k);
    for(int i=0;i<=categories;++i) {
        double v[DEGREE]={0},x=z[i];
        if(isfinite(x)) {
            v[0]=exp(-.5*x*x)/sqrt(2.*acos(-1.));v[1]=x*v[0];
            for(int k=2;k<DEGREE;++k) v[k]=(x*v[k-1]-roots[k-1]*v[k-2])/roots[k];
        }
        if(i) {
            Marginal *b=&bins[i-1];
            b->lo=phi(z[i-1],1);b->hi=phi(z[i],1);
            b->mass=cdf[i]-cdf[i-1]; b->variance=b->mass*(1.-b->mass);
            for(int k=0;k<DEGREE;++k) {
                b->c[k]=(last[k]-v[k])/roots[k+1];
                b->variance-=b->c[k]*b->c[k];
            }
            b->variance=fmax(0.,b->variance)+1e-14;
        }
        memcpy(last,v,sizeof last);
    }
    for(int j=0;j<n;++j) {
        const Marginal *a=&bins[left[j]],*b=&bins[right[j]];
        Cell *s=&cells[j];
        s->base=a->mass*b->mass;s->scale=sqrt(a->variance*b->variance);
        s->upper=fmin(a->mass,b->mass)+5e-14;
        for(int k=0;k<DEGREE;++k) s->c[k]=a->c[k]*b->c[k];
    }
    #pragma omp parallel if(rows*n>=16384) num_threads(4)
    {
        Entry cache[1024]={0};
        #pragma omp for schedule(static)
        for(int i=0;i<rows;++i) {
            double lsum=0.,hsum=0.;
            for(int j=0;j<n;++j) {
                const Cell *s=&cells[j];double r=rho[i*n+j];
                double r2=r*r,r4=r2*r2,r8=r4*r4;
                double lane[4];
                for(int t=0;t<4;++t)
                    lane[t]=s->c[t]+r4*(s->c[t+4]+r4*(s->c[t+8]+r4*s->c[t+12]));
                double p=s->base+r*(lane[0]+r*(lane[1]+r*(lane[2]+r*lane[3])));
                double err=fabs(r)*r8*r8*s->scale+5e-14;
                if(p-err<=1.01e-10) {
                    int a=(int)left[j],b=(int)right[j];
                    double sd=sqrt(1.-r2);
                    p=reuse(cache,z[a+1],z[b+1],bins[a].hi,bins[b].hi,r,sd)
                      -reuse(cache,z[a],z[b+1],bins[a].lo,bins[b].hi,r,sd)
                      -reuse(cache,z[a+1],z[b],bins[a].hi,bins[b].lo,r,sd)
                      +reuse(cache,z[a],z[b],bins[a].lo,bins[b].lo,r,sd);
                    err=5e-14;
                    if(!(p>1.02e-10)) {lsum=-INFINITY;hsum=INFINITY;break;}
                }
                if(err<=p*1e-6) {
                    double val=weights[j]*log(p),pad=weights[j]*1.00001e-6;
                    lsum+=val-pad;hsum+=val+pad;
                } else {
                    lsum+=weights[j]*log(p-err);
                    hsum+=weights[j]*log(fmin(p+err,s->upper));
                }
            }
            low[i]=lsum-1e-8;high[i]=hsum+1e-8;
        }
    }
    free(cells);free(bins);return 1;
}
