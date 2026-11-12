"""Exact spatial covariance geometry with cached eigendecompositions."""
import numpy as np
from scipy.linalg import eigh
from scipy.spatial.distance import cdist
from threadpoolctl import threadpool_limits

class LazyShape:
    def __init__(self,model,multiplier,fraction,group):
        self.model=model;self.range_multiplier=multiplier;self.fraction=fraction;self.group=group
    @property
    def eigenvalues(self):
        if self.group==0:return np.ones(self.model.n)
        values,_=self.model.decompose(self.group)
        return 1-self.fraction+self.fraction*values


class LazyBases:
    def __init__(self,model):self.model=model
    def __getitem__(self,i):
        if i==0:raise RuntimeError('IID bootstrap must use exact identity shortcut')
        return self.model.decompose(i)[1]


class LazyGeometryModel:
    def __init__(self,xy,fractions,ranges):
        self.xy=np.asarray(xy,float);self.n=len(xy);self.ranges=tuple(ranges)
        distance=cdist(self.xy,self.xy)
        self.scale=float(np.median(np.partition(distance,8,axis=1)[:,8]))
        assert self.n>8 and self.scale>0
        self._decompositions={};self.bases=LazyBases(self)
        self.shapes=[LazyShape(self,0.,0.,0)]
        for group,multiplier in enumerate(ranges,1):
            self.shapes.extend(LazyShape(self,multiplier,f,group) for f in fractions if f!=0)
    def decompose(self,group):
        if group not in self._decompositions:
            distance=cdist(self.xy,self.xy)
            kernel=np.exp(-distance/(self.scale*self.ranges[group-1]));del distance
            with threadpool_limits(limits=1):values,basis=eigh(kernel,overwrite_a=True)
            assert values.min()>=-1e-9
            self._decompositions[group]=(values,basis)
        return self._decompositions[group]
