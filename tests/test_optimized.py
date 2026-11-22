"""Replay the public optimized entry point against the reference count backend."""
import platform
from pathlib import Path
import unittest

import numpy as np

from lucent.calibration import refit
from lucent.calibration.backend import CPUBackend
from lucent.calibration.copula import CountGeometry, GRID_FRACTIONS, GRID_RANGES
from lucent.calibration.geometry import LazyGeometryModel
from lucent.calibration.optimized import adjust_family, optimized_session
from lucent.calibration.optimized.sequential import two_stage
from lucent.spatial import SpatialData, prepare_hypothesis


@unittest.skipUnless(platform.system() == "Linux" and platform.machine() == "x86_64"
                     and "avx2" in Path("/proc/cpuinfo").read_text().split(),
                     "Linux x86-64 AVX2 kernels")
class OptimizedReplay(unittest.TestCase):
    def test_two_stage_group(self):
        n = 36
        x = np.arange(n) % 6
        y = np.arange(n) // 6
        data = SpatialData(np.column_stack((x, x//2, x)), ["L1", "L2", "R"],
                           np.column_stack((x, y)), ["region"]*n,
                           np.arange(n), (np.arange(n)+6) % n)
        settings = dict(shallow=(15, 19), deep=(31, 39), target=4, cutoff=1.)
        with optimized_session() as session:
            for ids in (["L1"], ["L1", "L2"]):
                rows = [prepare_hypothesis(data, i, "R", "region", "region") for i in ids]
                first = rows[0]
                raw, v, o, xy = (first[k] for k in ("raw", "values", "outside", "coordinates"))
                model = LazyGeometryModel(xy, GRID_FRACTIONS, GRID_RANGES)
                be = CPUBackend(CountGeometry(xy, model), threads=1)
                def score(backend, fit, eps, values, outside, c):
                    return refit._score(c, refit.positive_ecdf(backend.bootstrap(fit, eps), values, outside))
                try:
                    ref = two_stage(raw, v, o, np.stack([row["contrast"] for row in rows]),
                                    ids, ("test", "R", "region"), be, be.fit(raw), refit, score, **settings)
                finally:
                    be.close()
                out = session.calibrate_group(data, "test", "R", "region",
                                              {i: (i, "region") for i in ids}, **settings)
                self.assertEqual(out, ref)
                family = adjust_family([*ids, "unavailable"], ids, out)
                self.assertEqual(family["p"][-1], 1.)
