"""Replay the public optimized entry point against the reference count backend."""
import platform
import json
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from lucent.calibration import refit
from lucent.calibration.backend import CPUBackend
from lucent.calibration.copula import CountGeometry, GRID_FRACTIONS, GRID_RANGES
from lucent.calibration.geometry import LazyGeometryModel
from lucent.calibration.optimized import adjust_family, optimized_session
from lucent.calibration.optimized.sequential import two_stage
from lucent.calibration.optimized.batch import _deep_results
from lucent.cli import main
from lucent.spatial import SpatialData, prepare_hypothesis


@unittest.skipUnless(platform.system() == "Linux" and platform.machine() == "x86_64"
                     and "avx2" in Path("/proc/cpuinfo").read_text().split(),
                     "Linux x86-64 AVX2 kernels")
class OptimizedReplay(unittest.TestCase):
    def test_retirement_invalidates_lookahead_in_the_same_block(self):
        # Deliberately make the second contrast's score depend on the active set.
        # Its valid hits occur only after the first contrast has stopped.
        def block(key, indices, k, r, begin, end):
            flags = [int(i == 0 or len(indices) == 1) for i in indices]
            return np.tile(flags, (end-begin, 1))

        with ThreadPoolExecutor(2) as pool, patch("lucent.calibration.optimized.batch._block", block):
            results = _deep_results(pool, ("R", "region"), np.array([0, 1]), np.array([0, 0]),
                                    [], (15, 9), 2, 2, 4)
        self.assertEqual([(r.hits, r.draws) for r in results], [(2, 2), (2, 4)])

    def test_batch_cli_matches_serial_with_short_final_block(self):
        n = 36
        x = np.arange(n) % 6
        labels = ["region"]*n
        data = SpatialData(np.column_stack((x, x//2, x, x//2)), ["L1", "L2", "R", "R2"],
                           np.column_stack((x, np.arange(n)//6)), labels,
                           np.arange(n), (np.arange(n)+6) % n)
        groups = {"R": {3: ("L1", "region"), 4: ("L2", "region")}, "R2": {7: ("L1", "region")}}
        settings = dict(shallow=(15, 19), deep=(31, 39), target=3, cutoff=1.)
        with optimized_session() as session:
            reference = [row for receptor, hs in groups.items() for row in session.calibrate_group(
                data, "test", receptor, "region", hs, **settings)]
        config = dict(section="test", **settings, hypotheses=[
            dict(id=i, receptor=r, receiver_label="region", ligand=l, sender_label=s)
            for r, hs in groups.items() for i, (l, s) in hs.items()])
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            np.savez(root/"data.npz", **vars(data))
            (root/"config.json").write_text(json.dumps(config))
            self.assertEqual(main(["batch","--input",str(root/"data.npz"),"--config",str(root/"config.json"),
                                   "--output",str(root/"rows.json"),"--workers","2","--block-size","4"]), 0)
            rows = json.loads((root/"rows.json").read_bytes())["rows"]
        self.assertEqual(rows, [dict(r, shallow=asdict(r["shallow"]), deep=asdict(r["deep"])) for r in reference])

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
