/* Preserve the 256-point check, including its original summation order.
   Elide erfc(+/-infinity) calls when one bound spans an open tail. */
#include <immintrin.h>
#include <math.h>
#include <stdint.h>

extern __m256d _ZGVdN4v_erfc(__m256d);

void integrate(int rows,int n,int order,int threads,const int64_t *index,
               const double *nodes,const double *weights,const double *width,
               const double *other,const double *lower,const double *upper,
               const double *rho,double *out) {
    const __m256d scale=_mm256_set1_pd(-M_SQRT1_2),half=_mm256_set1_pd(.5);
    #pragma omp parallel for if(threads>1) num_threads(threads) schedule(static)
    for(int i=0;i<rows*n;++i) {
        int j=i%n;
        double r=rho[i],sum=0.;
        if(r==0.) {out[i]=width[j]*other[j];continue;}
        __m256d sd=_mm256_set1_pd(sqrt(1.-r*r)),corr=_mm256_set1_pd(r);
        __m256d lb=_mm256_set1_pd(lower[j]),ub=_mm256_set1_pd(upper[j]);
        const double *v=nodes+index[j]*order;
        for(int k=0;k<order;k+=4) {
            __m256d shift=corr*_mm256_loadu_pd(v+k),p;
            if(lower[j]==-INFINITY) {
                __m256d hi=(ub-shift)/sd;
                p=_ZGVdN4v_erfc(hi*scale)*half;
            } else if(upper[j]==INFINITY) {
                __m256d lo=(lb-shift)/sd;
                __m256d pos=_mm256_cmp_pd(lo,_mm256_setzero_pd(),_CMP_GT_OQ);
                __m256d val=_ZGVdN4v_erfc(_mm256_blendv_pd(lo,-lo,pos)*scale);
                p=_mm256_blendv_pd(_mm256_set1_pd(2.)-val,val,pos)*half;
            } else {
                __m256d lo=(lb-shift)/sd,hi=(ub-shift)/sd;
                __m256d pos=_mm256_cmp_pd(lo,_mm256_setzero_pd(),_CMP_GT_OQ);
                __m256d a=_mm256_blendv_pd(hi,-lo,pos);
                __m256d b=_mm256_blendv_pd(lo,-hi,pos);
                p=(_ZGVdN4v_erfc(a*scale)-_ZGVdN4v_erfc(b*scale))*half;
            }
            double terms[4];
            _mm256_storeu_pd(terms,p*_mm256_loadu_pd(weights+k));
            for(int t=0;t<4;++t) sum+=terms[t];
        }
        out[i]=width[j]*sum;
    }
}

#include <immintrin.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
enum { D=96, Q=256, NS=512, NE=512 };
typedef struct {uint64_t key; int used; double v[D];} Endpoint;
typedef struct {int used;double nodes[Q],mu[D],bound;} Span;
typedef struct {double roots[D+1],w[Q];int weights_set;Span spans[NS];Endpoint ends[NE];
                uint64_t calls,expanded,direct,hits,misses;} Context;
typedef struct {double mu[D],bound;} Moment;
typedef struct {double v[D],tail_scale;} Terms;
typedef struct {uint64_t a,b;int span,index,used;} BinKey;

void *create_ctx(void) {
    Context *c=calloc(1,sizeof(Context));
    if(c) for(int k=0;k<=D;++k)c->roots[k]=sqrt((double)k);
    return c;
}
void free_ctx(void *c) {free(c);}
void context_stats(Context *c,uint64_t *s) {
    s[0]=c->calls;s[1]=c->expanded;s[2]=c->direct;s[3]=c->hits;s[4]=c->misses;
}
static uint64_t as_bits(double x){uint64_t b;memcpy(&b,&x,8);return b;}
static uint64_t hash(uint64_t x){return (x^(x>>31))*0x9e3779b97f4a7c15ULL;}

static void moments(Context *c,const double *nodes,const double *w,Moment *dst) {
    Span *s=&c->spans[hash(as_bits(nodes[0])+3*as_bits(nodes[Q-1]))>>55];
    if(!s->used || memcmp(s->nodes,nodes,Q*sizeof(double))) {
        ++c->misses;
        double prev[Q],cur[Q],sum=0.,wsum=0.;
        for(int j=0;j<Q;++j) {prev[j]=1.;cur[j]=nodes[j];sum+=w[j]*nodes[j];wsum+=fabs(w[j]);}
        s->mu[0]=sum;
        for(int k=2;k<=D;++k) {
            sum=0.;
            for(int j=0;j<Q;++j) {
                double v=(nodes[j]*cur[j]-c->roots[k-1]*prev[j])/c->roots[k];
                prev[j]=cur[j];cur[j]=v;sum+=w[j]*v;
            }
            s->mu[k-1]=sum;
        }
        double z=fmax(fabs(nodes[0]),fabs(nodes[Q-1]));
        const double t=.875;
        s->bound=wsum*pow(1.-t*t,-.25)*exp(t*z*z/(2.*(1.+t)));
        memcpy(s->nodes,nodes,Q*sizeof(double));s->used=1;
    } else ++c->hits;
    memcpy(dst->mu,s->mu,sizeof dst->mu);dst->bound=s->bound;
}
static void endpoint(Context *c,double x,double *out) {
    uint64_t key=as_bits(x);Endpoint *e=&c->ends[hash(key)>>55];
    if(!e->used || e->key!=key) {
        if(!isfinite(x)) memset(e->v,0,sizeof e->v);
        else {
            double prev=0.,cur=exp(-.5*x*x)/sqrt(2.*acos(-1.));
            for(int k=1;k<=D;++k) {
                e->v[k-1]=cur/c->roots[k];
                double next=(x*cur-c->roots[k-1]*prev)/c->roots[k];
                prev=cur;cur=next;
            }
        }
        e->used=1;e->key=key;
    }
    memcpy(out,e->v,sizeof e->v);
}
static double polynomial(const double *v,int degree,double r) {
    double r2=r*r,r4=r2*r2,r8=r4*r4,lane[8];
    memcpy(lane,v+degree-8,sizeof lane);
    for(int k=degree-16;k>=0;k-=8)
        for(int j=0;j<8;++j) lane[j]=v[k+j]+r8*lane[j];
    double val=lane[7];
    for(int j=6;j>=0;--j) val=lane[j]+r*val;
    return val;
}
/* integrate() is the original 256-node fallback, built in this translation unit. */
void integrate(int,int,int,int,const int64_t*,const double*,const double*,const double*,
               const double*,const double*,const double*,const double*,double*);

void moment_integrate(Context *c,int rows,int n,int order,int threads,const int64_t *index,
                      const double *nodes,const double *w,const double *width,const double *other,
                      const double *lower,const double *upper,const double *rho,double *out) {
    int num=0;
    for(int j=0;j<n;++j) if(index[j]>=num) num=(int)index[j]+1;
    Moment *m=malloc((size_t)num*sizeof(*m));
    Terms *terms=malloc((size_t)n*sizeof(*terms));
    int *ids=malloc((size_t)n*sizeof(*ids));
    if(!m || !terms || !ids || order!=Q) {
        integrate(rows,n,order,threads,index,nodes,w,width,other,lower,upper,rho,out);
        free(m);free(terms);free(ids);return;
    }
    if(!c->weights_set || memcmp(c->w,w,sizeof c->w)) {
        memset(c->spans,0,sizeof c->spans);memcpy(c->w,w,sizeof c->w);c->weights_set=1;
    }
    for(int s=0;s<num;++s) moments(c,nodes+s*Q,w,m+s);
    BinKey keys[512]={0};
    int count=0;
    for(int j=0;j<n;++j) {
        uint64_t a=as_bits(lower[j]),b=as_bits(upper[j]);
        int s=(int)index[j];
        BinKey *key=&keys[hash(a+7*hash(b)+31*(uint64_t)s)>>55];
        if(key->used && key->a==a && key->b==b && key->span==s) {ids[j]=key->index;continue;}
        ids[j]=count++;*key=(BinKey){a,b,s,ids[j],1};
        double l[D],h[D]; endpoint(c,lower[j],l);endpoint(c,upper[j],h);
        Terms *v=&terms[ids[j]];
        for(int k=0;k<D;++k)v->v[k]=(l[k]-h[k])*m[s].mu[k];
        v->tail_scale=m[s].bound*sqrt(other[j]*(1.-other[j]));
    }
    uint64_t used=0,original=0;
    const __m256d scale=_mm256_set1_pd(-M_SQRT1_2),half=_mm256_set1_pd(.5);
    #pragma omp parallel for if(threads>1) num_threads(threads) schedule(static) reduction(+:used,original)
    for(int i=0;i<rows*n;++i) {
        int j=i%n;double r=rho[i],sum=0.;
        if(r==0.) {out[i]=width[j]*other[j];continue;}
        if(fabs(r)<.72 && width[j]>=1.e-6 && other[j]>=1.e-6) {
            const Terms *v=&terms[ids[j]];
            double q=fabs(r)/sqrt(.875),q2=q*q,q4=q2*q2,q8=q4*q4;
            double power=q8,limit=2.e-15*other[j];
            int k=8;
            while(k<64 && power*q*v->tail_scale>limit) {power*=power;k*=2;}
            if(k==64 && power*q*v->tail_scale>limit) {double q16=q8*q8;power*=q16*q16;k=96;}
            double err=power*q*v->tail_scale;
            if(err<=limit) {
                double p=other[j]+r*polynomial(v->v,k,r);
                if(width[j]*p>1.01e-10 && err<=5e-15*p) {out[i]=width[j]*p;++used;continue;}
            }
        }
        ++original;
                __m256d sd=_mm256_set1_pd(sqrt(1.-r*r)),corr=_mm256_set1_pd(r);
        __m256d lb=_mm256_set1_pd(lower[j]),ub=_mm256_set1_pd(upper[j]);
        const double *v=nodes+index[j]*order;
        for(int k=0;k<order;k+=4) {
            __m256d shift=corr*_mm256_loadu_pd(v+k),p;
            if(lower[j]==-INFINITY) {
                __m256d hi=(ub-shift)/sd;
                p=_ZGVdN4v_erfc(hi*scale)*half;
            } else if(upper[j]==INFINITY) {
                __m256d lo=(lb-shift)/sd;
                __m256d pos=_mm256_cmp_pd(lo,_mm256_setzero_pd(),_CMP_GT_OQ);
                __m256d val=_ZGVdN4v_erfc(_mm256_blendv_pd(lo,-lo,pos)*scale);
                p=_mm256_blendv_pd(_mm256_set1_pd(2.)-val,val,pos)*half;
            } else {
                __m256d lo=(lb-shift)/sd,hi=(ub-shift)/sd;
                __m256d pos=_mm256_cmp_pd(lo,_mm256_setzero_pd(),_CMP_GT_OQ);
                __m256d a=_mm256_blendv_pd(hi,-lo,pos);
                __m256d b=_mm256_blendv_pd(lo,-hi,pos);
                p=(_ZGVdN4v_erfc(a*scale)-_ZGVdN4v_erfc(b*scale))*half;
            }
            double terms[4];
            _mm256_storeu_pd(terms,p*_mm256_loadu_pd(w+k));
            for(int t=0;t<4;++t) sum+=terms[t];
        }
        out[i]=width[j]*sum;
    }
    c->calls++;c->expanded+=used;c->direct+=original;
    free(m);free(terms);free(ids);
}
