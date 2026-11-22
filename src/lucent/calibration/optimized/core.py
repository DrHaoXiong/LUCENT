"""Fixed-maximum screening and deep-tail sampling. Counts use inclusive tie handling."""
from fractions import Fraction
import hashlib
import json
import numpy as np

def seed(role,*key):
    return int.from_bytes(hashlib.sha256(json.dumps(["RC1_DEEP_V1",role,*key]).encode()).digest()[:16],"little")

def screen_limit(budget,alpha="0.02"):
    return int(Fraction(alpha)*(budget+1))

def screen(stream,budget=1999,alpha="0.02"):
    limit=screen_limit(budget,alpha)
    hits=0
    for done,event in enumerate(stream,1):
        hits+=int(event)
        if hits>=limit:
            return False,done
        if done==budget:
            return True,done
    raise ValueError("incomplete screening stream")

def observed(raw,v,o,c,be,fit,hseed,r,rc,score):
    t=float(rc._score(c,rc.positive_ecdf(raw[:,None],v,o))[0])
    if t<=1e-12: raise ValueError("expected eligible positive contrast")
    rng=rc._rng(hseed,"observed_inner")
    k=ties=0
    for offset in range(0,r,64):
        m=min(64,r-offset)
        x=score(be,fit,rng.normal(size=(m,len(raw))).T,v,o,c)
        k+=int(np.count_nonzero(x>=t))
        ties+=int(np.count_nonzero(x==t))
    return t,k,ties

def rounds(raw,v,o,coeff,be,fit,gseed,cutoffs,r,begin,end,rc,score):
    k=np.asarray(cutoffs)
    hits=np.zeros(len(k),np.int64)
    equals=hits.copy()
    # Equality of the final inner count is only known when all R draws were used.
    complete=hits.copy()
    physical=0
    for u in range(begin,end):
        outer=be.bootstrap(fit,rc._rng(gseed,"outer",u).normal(size=(len(raw),1)))
        y=rc.positive_ecdf(outer,v,o)
        t=np.array([rc._score(c,y)[0] for c in coeff])
        active=np.flatnonzero((t>1e-12)&(np.ptp(y)>1e-15))
        hits+=((t<=1e-12)|(np.ptp(y)<=1e-15))&(k==r)
        if not len(active): continue
        fb=be.fit(outer[:,0])
        rng=rc._rng(gseed,"inner",u)
        counts=np.zeros(len(k),np.int64)
        done=0
        while len(active):
            size=min(64,r-done)
            eps=rng.normal(size=(size,len(raw))).T
            if len(active)==1:
                j=active[0]
                x=score(be,fb,eps,v,o,coeff[j])[None,:]
            else:
                field=np.ascontiguousarray(rc.positive_ecdf(be.bootstrap(fb,eps),v,o))
                x=np.einsum("ij,jk->ik",coeff[active],field,optimize=False)
            counts[active]+=np.count_nonzero(x>=t[active,None],axis=1)
            physical+=size
            done+=size
            stop=(counts[active]>k[active])|(counts[active]+r-done<=k[active])
            idx=active[stop]
            hits[idx]+=counts[idx]<=k[idx]
            if done==r:
                complete[idx]+=1
                equals[idx]+=counts[idx]==k[idx]
            active=active[~stop]
    return hits.tolist(),equals.tolist(),complete.tolist(),physical
