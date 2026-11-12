"""Thread independent probability rows; preserve quadrature, fit and Gaussian draws."""
from concurrent.futures import ThreadPoolExecutor
import time
import numpy as np
from scipy.special import ndtr
from threadpoolctl import ThreadpoolController
from .rectangles import fast_rectangle
from .copula import rectangle

class ThreadedBackend:
    def __init__(self,geom,threads=4):
        self.geom=geom;self.threads=threads;self.executor=ThreadPoolExecutor(max_workers=threads);self.controller=ThreadpoolController()
        self.distance_ids=np.unique(geom.distance,return_inverse=True)[1];self.fit_seconds=self.bootstrap_seconds=0.;self.fits=self.fallbacks=0
    def probability(self,func,al,ah,bl,bh,r,order=None):
        args=np.broadcast_arrays(al,ah,bl,bh,r);shape=args[-1].shape;size=args[-1].size
        threshold=1024 if order else 16384
        if self.threads<=1 or size<threshold:return func(*args) if order is None else func(*args,order)
        flat=[a.ravel() for a in args];cuts=np.linspace(0,size,self.threads+1,dtype=int)
        def part(lo,hi):
            a=[x[lo:hi] for x in flat];return func(*a) if order is None else func(*a,order)
        futures=[self.executor.submit(part,int(lo),int(hi)) for lo,hi in zip(cuts[:-1],cuts[1:])]
        return np.concatenate([f.result() for f in futures]).reshape(shape)
    def fit(self,raw):
        t=time.perf_counter();self.fits+=1;values,labels,counts=np.unique(raw,return_inverse=True,return_counts=True)
        if len(values)<2:raise ValueError('RESPONSE_CONSTANT')
        cumulative=np.r_[0.,np.cumsum(counts)/len(raw)]
        left,right=labels[self.geom.pairs[:,0]],labels[self.geom.pairs[:,1]]
        keys=np.column_stack([self.distance_ids,np.minimum(left,right),np.maximum(left,right)])
        _,first,weight=np.unique(keys,axis=0,return_index=True,return_counts=True);l,r=left[first],right[first]
        al,ah,bl,bh=cumulative[l],cumulative[l+1],cumulative[r],cumulative[r+1];rho=self.geom.rho[:,first]
        ll=np.log(self.probability(fast_rectangle,al[None,:],ah[None,:],bl[None,:],bh[None,:],rho))@weight
        best=np.argsort(-ll,kind='stable')[:2]
        check=np.log(self.probability(rectangle,al[None,:],ah[None,:],bl[None,:],bh[None,:],rho[best],256))@weight
        drift=np.abs(check-ll[best])
        if drift.max()>.02 or best[np.argmax(check)]!=best[0] or ll[best[0]]-ll[best[1]]<1e-7:
            self.fallbacks+=1;result=self.geom.fit(raw)
        else:
            chosen=int(best[0]);shape=self.geom.model.shapes[chosen]
            result={'shape':chosen,'fraction':shape.fraction,'range_multiplier':shape.range_multiplier,'values':values,'counts':counts,'cumulative':cumulative,'log_likelihood':float(ll[chosen]),'pair_count':len(left),'quadrature_order':0,'quadrature_loglik_drift':float(drift.max())}
        self.fit_seconds+=time.perf_counter()-t;return result
    def bootstrap(self,fit,eps):
        t=time.perf_counter();shape=self.geom.model.shapes[fit['shape']]
        if shape.group==0:z=eps
        else:
            with self.controller.limit(limits=8):z=self.geom.model.bases[shape.group]@(np.sqrt(shape.eigenvalues)[:,None]*eps)
        labels=np.searchsorted(fit['cumulative'][1:],ndtr(z),side='left');result=fit['values'][np.minimum(labels,len(fit['values'])-1)]
        self.bootstrap_seconds+=time.perf_counter()-t;return result
    def close(self):self.executor.shutdown(wait=True)
