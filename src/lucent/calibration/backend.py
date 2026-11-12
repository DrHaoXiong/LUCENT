"""CPU-budgeted native sampler for the bounded CAD Visium receiver pools."""
import time
import numpy as np
from scipy.special import ndtr
from .threaded import ThreadedBackend

class CPUBackend(ThreadedBackend):
    def __init__(self,geom,threads=4):
        super().__init__(geom,threads);self.spectral={}
    def bootstrap(self,fit,eps):
        t=time.perf_counter();index=fit['shape'];shape=self.geom.model.shapes[index]
        if shape.group==0:z=eps
        else:
            if index not in self.spectral:self.spectral[index]=(self.geom.model.bases[shape.group],np.sqrt(shape.eigenvalues))
            basis,root=self.spectral[index]
            # Small native matrix products use one BLAS thread; fit quadrature rows use four.
            with self.controller.limit(limits=1 if self.geom.n<=768 else 8):z=basis@(root[:,None]*eps)
        labels=np.searchsorted(fit['cumulative'][1:],ndtr(z),side='left');result=fit['values'][np.minimum(labels,len(fit['values'])-1)]
        self.bootstrap_seconds+=time.perf_counter()-t;return result
