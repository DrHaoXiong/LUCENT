"""Optional float64 inner field kernel for a frozen receiver pool."""
import ctypes
from pathlib import Path
import subprocess

import numpy as np
from scipy.special import cython_special


def _pointer(array):
    return array.ctypes.data


def build(directory):
    source = Path(__file__).with_name("cutoff.c")
    path = Path(directory) / "layer4.so"
    subprocess.run(["cc", "-std=c11", "-O3", "-ffp-contract=off", "-fno-fast-math",
                    "-fPIC", "-shared", str(source), "-o", str(path), "-lm"], check=True)
    lib = ctypes.CDLL(str(path))
    p, i, long = ctypes.c_void_p, ctypes.c_int, ctypes.c_int64
    for name, args in [
        ("bind_normal", [p,p]), ("prepare_guards", [i,p,p,p]),
        ("labels", [i,p,i,p,p,p,p]), ("transform", [i]*4+[long]*2+[p]*11)
    ]:
        function = getattr(lib, name)
        function.argtypes = args
        function.restype = None
    get = ctypes.pythonapi.PyCapsule_GetPointer
    get.argtypes = [ctypes.py_object, ctypes.c_char_p]
    get.restype = p
    lib.bind_normal(*(get(cython_special.__pyx_capi__[name], b"double (double, int __pyx_skip_dispatch)")
                      for name in ("__pyx_fuse_1ndtr","ndtri")))
    return lib


class Cutoffs:
    def __init__(self, lib, probabilities):
        self.lib = lib
        self.cum = np.ascontiguousarray(probabilities, dtype=np.float64)
        self.lo = np.empty_like(self.cum)
        self.hi = np.empty_like(self.cum)
        self.args = (len(self.cum), *map(_pointer, (self.cum, self.lo, self.hi)))
        lib.prepare_guards(*self.args)

    def labels(self, z):
        flat = np.ascontiguousarray(z, dtype=np.float64).ravel()
        result = np.empty(len(flat), dtype=np.int64)
        self.lib.labels(len(flat), _pointer(flat), *self.args, _pointer(result))
        return result.reshape(np.shape(z))


class FieldKernel:
    """One instance per receiver group/worker. project supplies the original Gaussian transform."""
    def __init__(self, lib, project, values, outside, *, small_first=False):
        self.lib = lib
        self.project = project
        self.values = np.asarray(values)
        self.outside = np.ascontiguousarray(outside, dtype=np.int64)
        if (self.values.ndim != 1 or len(self.values) == 0 or self.values[0] != 0
                or self.outside.shape != self.values.shape or np.any(np.diff(self.values) <= 0)
                or np.any(self.outside < 0) or not np.array_equal(self.outside, outside)):
            raise ValueError("Expected ordered full-section counts starting with zero")
        self._source_outside = outside
        self.small_first = small_first
        self._fit = None
        self._bins = np.empty(0, dtype=np.int64)
        self._counts = np.empty_like(self.outside)
        self._ranks = np.empty(len(self.values))

    def _prepare(self, fit):
        if fit is self._fit:
            return
        support = fit["values"]
        cum = np.asarray(fit["cumulative"])
        if (len(support) == 0 or cum.shape != (len(support)+1,) or
                not np.isfinite(cum).all() or np.any(np.diff(cum)<0) or
                cum[0] != 0 or np.any(cum[1:-1]>1)):
            raise ValueError("Invalid fitted count distribution")
        mapping = np.searchsorted(self.values, support)
        if np.any(mapping >= len(self.values)) or not np.array_equal(self.values[mapping], support):
            raise ValueError("Receiver count outside full-section support")
        self._map = np.ascontiguousarray(mapping,dtype=np.int64)
        # The original native search clamps indices to the last support value.
        self._cuts = Cutoffs(self.lib, cum[1:-1])
        self._fit = fit

    def transform(self, fit, z, contrast=None):
        self._prepare(fit)
        z = np.asarray(z, dtype=np.float64)
        n, batches = z.shape
        if len(self._bins) != n:
            self._bins = np.empty(n, dtype=np.int64)
        if contrast is None:
            out = np.empty((n, batches))
            c = None
        else:
            contrast = np.ascontiguousarray(contrast,dtype=np.float64)
            if contrast.shape != (n,):
                raise ValueError("Contrast does not match receiver field")
            c = _pointer(contrast)
            out = np.empty(batches)
        self.lib.transform(n,batches,len(self.values),len(self._cuts.cum),
                           z.strides[0]//8,z.strides[1]//8,_pointer(z),*self._cuts.args[1:],
                           _pointer(self._map),_pointer(self.outside),c,_pointer(self._bins),
                           _pointer(self._counts),_pointer(self._ranks),_pointer(out))
        return out

    def _pool(self, values, outside):
        if ((values is not self.values and not np.array_equal(values,self.values)) or
                (outside is not self._source_outside and not np.array_equal(outside,self.outside))):
            raise ValueError("Kernel is bound to a different receiver pool")

    def score(self, backend, fit, eps, values, outside, contrast):
        self._pool(values,outside)
        return self.transform(fit,self.project(backend,fit,eps),contrast)

    def rounds(self, raw, values, outside, coeff, backend, fit, group_seed, cutoffs,
               r, begin, end, rc, score):
        self._pool(values,outside)
        k = np.asarray(cutoffs)
        hits = np.zeros(len(k),dtype=np.int64)
        equals = np.zeros_like(hits)
        complete = np.zeros_like(hits)
        physical = 0
        for u in range(begin,end):
            eps = rc._rng(group_seed,"outer",u).normal(size=(len(raw),1))
            outer = backend.bootstrap(fit,eps)
            field = rc.positive_ecdf(outer,values,outside)
            t = np.array([rc._score(c,field)[0] for c in coeff])
            inactive = (t<=1e-12) | (np.ptp(field)<=1e-15)
            hits += inactive & (k==r)
            active = np.flatnonzero(~inactive)
            if not len(active):
                continue
            refit = backend.fit(outer[:,0])
            rng = rc._rng(group_seed,"inner",u)
            counts = np.zeros_like(k)
            done = 0
            while len(active):
                size = 64
                if self.small_first and len(active)==1 and k[active[0]]==0 and done<56:
                    size = 8 if done==0 else (16 if done==8 else 32)
                size = min(size,r-done)
                eps = rng.normal(size=(size,len(raw))).T
                projected = self.project(backend,refit,eps)
                if len(active)==1:
                    x = self.transform(refit,projected,coeff[active[0]])[None,:]
                else:
                    ranks = self.transform(refit,projected)
                    x = np.einsum("ij,jk->ik",coeff[active],ranks,optimize=False)
                counts[active] += np.count_nonzero(x>=t[active,None],axis=1)
                done += size
                physical += size
                stop = (counts[active]>k[active]) | (counts[active]+r-done<=k[active])
                idx = active[stop]
                hits[idx] += counts[idx]<=k[idx]
                if done==r:
                    complete[idx] += 1
                    equals[idx] += counts[idx]==k[idx]
                active = active[~stop]
        return hits.tolist(),equals.tolist(),complete.tolist(),physical
