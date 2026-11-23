"""Grouped two-stage RC1 calibration on Linux x86-64 with AVX2 and glibc."""
from contextlib import contextmanager
import ctypes as ct
import platform
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace

import numpy as np
from .. import refit
from ..backend import CPUBackend
from ..copula import CountGeometry, GRID_FRACTIONS, GRID_RANGES
from ..geometry import LazyGeometryModel
from . import field_kernel
from .fit import Moments, Prepared, capsule, make_fit, make_probability
from .sequential import two_stage
from .sequential_family import adjust_family
from .batch import calibrate_batch

__all__ = ["optimized_session", "adjust_family", "calibrate_batch"]


def _library(directory, name, *, vector=False):
    src = Path(__file__).with_name(name + ".c")
    output = directory / (name + ".so")
    options = ["-mavx2"] if vector else []
    libraries = ["-lmvec", "-lm"] if vector else ["-lm"]
    subprocess.run(["cc", "-O3", "-fPIC", "-shared", "-fopenmp", "-fno-fast-math",
                    "-ffp-contract=off", *options, str(src), *libraries, "-o", str(output)],
                   check=True)
    return ct.CDLL(str(output))


def _function(lib, name, args, result=None):
    fn = getattr(lib, name)
    fn.argtypes = args
    fn.restype = result
    return fn


class _Session:
    def __init__(self, directory):
        self.grid = _library(directory, "grid")
        self.bounds = _library(directory, "bounds")
        self.quad = _library(directory, "quadrature", vector=True)
        self.fields = field_kernel.build(directory)
        ptr, integer = ct.c_void_p, ct.c_int
        for lib in (self.grid, self.bounds):
            bind = _function(lib, "bind_special", [ptr, ptr])
            bind(capsule("__pyx_fuse_1ndtr", "double (double, int __pyx_skip_dispatch)"),
                 capsule("owens_t", "double (double, double, int __pyx_skip_dispatch)"))
        self.row_bounds = _function(self.bounds, "row_bounds", [integer]*3+[ptr]*10+[integer], integer)
        grid = _function(self.grid, "adaptive_rectangles", [integer]*3+[ptr]*4)
        abi = [integer]*4+[ptr]*9
        _function(self.quad, "moment_integrate", [ptr]+abi)
        _function(self.quad, "create_ctx", [], ptr)
        _function(self.quad, "free_ctx", [ptr])
        _function(self.quad, "context_stats", [ptr, ptr])
        direct = _function(self.quad, "integrate", abi)
        self.direct = SimpleNamespace(rectangles=grid, integrate=direct)
        self.grid_prob = make_probability(self.direct)

    def calibrate_group(self, data, section, receptor, receiver_label, hypotheses, *,
                        large_pool_blas_threads=8, **settings):
        """hypotheses maps distinct positive hypothesis IDs to (ligand, sender label)."""
        with self._prepare(data, receptor, receiver_label, hypotheses, large_pool_blas_threads) as group:
            raw, values, outside, contrasts, backend, fitted, fields = group
            return two_stage(raw, values, outside, contrasts, list(hypotheses),
                             (section, receptor, receiver_label), backend, fitted, refit,
                             fields.score, round_sampler=fields.rounds, **settings)

    @contextmanager
    def _prepare(self, data, receptor, receiver_label, hypotheses, large_pool_blas_threads):
        from ...spatial import prepare_hypothesis

        if not hypotheses:
            raise ValueError("Expected at least one positive hypothesis")
        rows = [prepare_hypothesis(data, ligand, receptor, sender, receiver_label)
                for ligand, sender in hypotheses.values()]
        first = rows[0]
        if len(first["raw"]) <= 8 or any((r["summary"]["centered_score"] or 0) <= 1e-12
                                         for r in rows):
            raise ValueError("Grouped calibration requires a full pool of at least 9 and positive contrasts")
        raw, values, outside = (first[k] for k in ("raw", "values", "outside"))
        xy = first["coordinates"]
        model = LazyGeometryModel(xy, GRID_FRACTIONS, GRID_RANGES)
        backend = CPUBackend(CountGeometry(xy, model), threads=1,
                             large_pool_blas_threads=large_pool_blas_threads)
        moments = Moments(self.quad)
        try:
            high = Prepared(len(raw), self.direct, moments)
            def probability(func, al, ah, bl, bh, rho, order=None):
                if order is None:
                    return self.grid_prob(func, al, ah, bl, bh, rho)
                return high(func, al, ah, bl, bh, rho, order)
            backend.probability = probability
            fit = make_fit(self.row_bounds, dict(fits=0, evaluated=0, rechecks=0))
            backend.fit = fit.__get__(backend)
            fields = field_kernel.FieldKernel(self.fields, CPUBackend.project, values, outside,
                                               small_first=len(rows)==1)
            contrasts = np.stack([row["contrast"] for row in rows])
            yield raw, values, outside, contrasts, backend, backend.fit(raw), fields
        finally:
            moments.close()
            backend.close()


@contextmanager
def optimized_session():
    """Build the native kernels in a temporary directory once per worker."""
    if platform.system() != "Linux" or platform.machine() not in {"x86_64", "AMD64"}:
        raise RuntimeError("Native RC1 requires x86-64 Linux (including WSL2)")
    if "avx2" not in Path("/proc/cpuinfo").read_text().split():
        raise RuntimeError("Native RC1 requires AVX2")
    with tempfile.TemporaryDirectory(prefix="lucent-rc1-") as directory:
        yield _Session(Path(directory))
