"""Grid screening and 256-node verification for the fitted count copula."""
from collections import OrderedDict
import ctypes as ct
import time
import weakref

import numpy as np
from scipy.special import ndtr, ndtri, cython_special
from ..copula import quadrature, rectangle
from ..rectangles import fast_rectangle


def ptr(arr):
    return ct.c_void_p(arr.ctypes.data)

def capsule(name, sig):
    cp=cython_special.__pyx_capi__[name]
    get=ct.pythonapi.PyCapsule_GetPointer
    get.argtypes=[ct.py_object,ct.c_char_p];get.restype=ct.c_void_p
    return get(cp,sig.encode())

def make_probability(lib):
    def probability(func,al,ah,bl,bh,r,order=None):
        bounds=np.ascontiguousarray(np.stack([x.ravel() for x in (al,ah,bl,bh)]))
        r=np.ascontiguousarray(r)
        rows,n=r.shape
        out=np.empty_like(r)
        if order is None:
            z=ndtri(bounds);cdf=ndtr(z)
            lib.rectangles(rows,n,4 if r.size>=16384 else 1,ptr(z),ptr(cdf),ptr(r),ptr(out))
            bad=~np.isfinite(out)|(out<1e-10)
            if bad.any():
                ii,jj=np.nonzero(bad)
                out[ii,jj]=rectangle(*(a[jj] for a in bounds),r[ii,jj],256)
        else:
            al,ah,bl,bh=bounds
            swap=(ah-al)>(bh-bl)
            al,ah,bl,bh=(np.where(swap,bl,al),np.where(swap,bh,ah),np.where(swap,al,bl),np.where(swap,ah,bh))
            spans,inv=np.unique(np.stack((al,ah),axis=1),axis=0,return_inverse=True)
            u,w=quadrature(order)
            latent=ndtri(spans[:,0,None]+(spans[:,1]-spans[:,0])[:,None]*u)
            inputs=[np.ascontiguousarray(x) for x in (inv.astype(np.int64),latent,w,ah-al,bh-bl,ndtri(bl),ndtri(bh),r,out)]
            lib.integrate(rows,n,order,4 if r.size>=1024 else 1,*map(ptr,inputs))
            if not np.isfinite(out).all() or (out<=0).any(): raise ArithmeticError("invalid probability")
        return out
    return probability


class Prepared:
    def __init__(self,n,scalar,vector):
        self.n=n
        self.base=make_probability(scalar)
        self.vector=vector
        self.bound=np.arange(n+1,dtype=float)/n
        self.z=ndtri(self.bound)
        self.u,self.w=quadrature(256)
        self.nodes=OrderedDict()

    def __call__(self,func,al,ah,bl,bh,r,order=None):
        if order!=256:
            return self.base(func,al,ah,bl,bh,r,order)
        # The RC1 fitted empirical marginal has exactly n observations.
        bounds=np.stack([x.ravel() for x in (al,ah,bl,bh)])
        ids=np.rint(bounds*self.n).astype(np.int64)
        if not np.array_equal(self.bound[ids],bounds):
            return self.base(func,al,ah,bl,bh,r,order)
        ia,ib,ic,id_=ids
        # Preserve the original floating span comparison, including exact ties.
        swap=(bounds[1]-bounds[0])>(bounds[3]-bounds[2])
        a,b=np.where(swap,ic,ia),np.where(swap,id_,ib)
        c,d=np.where(swap,ia,ic),np.where(swap,ib,id_)
        stride=self.n+1
        unique,inv=np.unique(a*stride+b,return_inverse=True)
        latent=[]
        for item in unique:
            key=int(item)
            row=self.nodes.get(key)
            if row is None:
                lo,hi=divmod(key,stride)
                start,stop=self.bound[lo],self.bound[hi]
                row=ndtri(start+(stop-start)*self.u)
                self.nodes[key]=row
                if len(self.nodes)>2048:
                    self.nodes.popitem(last=False)
            else:
                self.nodes.move_to_end(key)
            latent.append(row)
        r=np.ascontiguousarray(r)
        out=np.empty_like(r)
        args=[np.ascontiguousarray(x) for x in
              (inv.astype(np.int64),np.stack(latent),self.w,self.bound[b]-self.bound[a],
               self.bound[d]-self.bound[c],self.z[c],self.z[d],r,out)]
        self.vector.integrate(*r.shape,256,4 if r.size>=1024 else 1,*map(ptr,args))
        if not np.isfinite(out).all() or (out<=0).any():
            raise ArithmeticError("invalid rectangle")
        return out


def make_fit(bounds,counters,tighten=False):
    def fit(self,raw):
        start=time.perf_counter()
        self.fits+=1
        values,labels,counts=np.unique(raw,return_inverse=True,return_counts=True)
        n=len(values)
        if n<2: raise ValueError("RESPONSE_CONSTANT")
        cumulative=np.r_[0.,np.cumsum(counts)/len(raw)]
        l,r=labels[self.geom.pairs[:,0]],labels[self.geom.pairs[:,1]]
        keys=(self.distance_ids*n+np.minimum(l,r))*n+np.maximum(l,r)
        _,first,w=np.unique(keys,return_index=True,return_counts=True)
        left,right=l[first],r[first]
        args=[cumulative[left][None,:],cumulative[left+1][None,:],
              cumulative[right][None,:],cumulative[right+1][None,:]]
        rho=np.ascontiguousarray(self.geom.rho[:,first])
        lo=np.empty(len(rho));hi=np.empty_like(lo)
        if tighten:
            half=.25*np.diff(cumulative)
            inner_l=ndtri(cumulative[:-1]+half)
            inner_u=ndtri(cumulative[1:]-half)
        else:
            inner_l=inner_u=cumulative
        ok=bounds(*rho.shape,n,*map(ptr,(ndtri(cumulative),cumulative,left,right,w,rho,lo,hi,
                                             inner_l,inner_u)),int(tighten))
        all_ids=np.arange(len(rho))
        keep=all_ids if not ok else np.flatnonzero(hi>=np.partition(lo,-2)[-2])
        if len(keep)<2: keep=all_ids
        def evaluate(ids):
            return np.log(self.probability(fast_rectangle,*args,rho[ids]))@w
        ll=np.full(len(rho),-np.inf)
        ll[keep]=evaluate(keep)
        best=np.argsort(-ll,kind="stable")[:2]
        if len(keep)<len(rho) and (np.any(ll[keep]<lo[keep]) or np.any(ll[keep]>hi[keep])
                or ll[best[0]]-ll[best[1]]<1.e-7):
            counters["rechecks"]+=1
            ll=evaluate(all_ids)
            best=np.argsort(-ll,kind="stable")[:2]
        counters["evaluated"]+=len(keep)
        counters["fits"]+=1
        check=np.log(self.probability(rectangle,*args,rho[best],256))@w
        drift=np.abs(check-ll[best])
        gap=ll[best[0]]-ll[best[1]]
        if drift.max()>.02 or best[np.argmax(check)]!=best[0] or gap<1.e-7:
            self.fallbacks+=1
            result=self.geom.fit(raw)
        else:
            chosen=int(best[0])
            shape=self.geom.model.shapes[chosen]
            result=dict(shape=chosen,fraction=shape.fraction,range_multiplier=shape.range_multiplier,
                        values=values,counts=counts,cumulative=cumulative,log_likelihood=float(ll[chosen]),
                        pair_count=len(l),quadrature_order=0,quadrature_loglik_drift=float(drift.max()))
        self.fit_seconds+=time.perf_counter()-start
        return result
    return fit

class Moments:
    def __init__(self,lib):
        self.lib=lib;self.ctx=lib.create_ctx()
        if not self.ctx: raise MemoryError("quadrature moment context")
        self.close=weakref.finalize(self,lib.free_ctx,self.ctx)
    def integrate(self,*args):
        self.lib.moment_integrate(self.ctx,*args)
    def stats(self):
        data=np.zeros(5,np.uint64)
        self.lib.context_stats(self.ctx,ptr(data))
        return data.tolist()
