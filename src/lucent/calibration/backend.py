"""CPU-budgeted native sampler for the bounded CAD Visium receiver pools."""
import time
import numpy as np
from scipy.special import ndtr
from .threaded import ThreadedBackend

class CPUBackend(ThreadedBackend):
    def __init__(self,geom,threads=4,*,large_pool_blas_threads=8):
        if type(large_pool_blas_threads) is not int or large_pool_blas_threads < 1:
            raise ValueError("large_pool_blas_threads must be a positive integer")
        super().__init__(geom,threads);self.spectral={}
        self.projection_threads = 1 if geom.n <= 768 else large_pool_blas_threads
    def project(self,fit,eps):
        index=fit['shape'];shape=self.geom.model.shapes[index]
        if shape.group==0:z=eps
        else:
            if index not in self.spectral:self.spectral[index]=(self.geom.model.bases[shape.group],np.sqrt(shape.eigenvalues))
            basis,root=self.spectral[index]
            # Small native matrix products use one BLAS thread; fit quadrature rows use four.
            with self.controller.limit(limits=self.projection_threads, user_api="blas"):z=basis@(root[:,None]*eps)
        return z
    def bootstrap(self,fit,eps):
        t=time.perf_counter();z=self.project(fit,eps)
        labels=np.searchsorted(fit['cumulative'][1:],ndtr(z),side='left');result=fit['values'][np.minimum(labels,len(fit['values'])-1)]
        self.bootstrap_seconds+=time.perf_counter()-t;return result
