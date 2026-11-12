"""Main raw-count-refit calibration with explicit SC baseline, version RC1.

The geometry backend fits counts and simulates counts. It receives no truth.
Decision-only acceleration returns exact bounds on the prescribed finite-MC p.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Any, Literal, Protocol

import numpy as np
from numpy.typing import NDArray

Array = NDArray[Any]


class CountBackend(Protocol):
    def fit(self, response: Array) -> dict[str, Any]: ...
    def bootstrap(self, fit: dict[str, Any], innovations: Array) -> Array: ...


@dataclass(frozen=True)
class CountCalibrationConfig:
    method: Literal["sc", "raw-count-refit"] = "raw-count-refit"
    inner_draws: int = 999
    outer_draws: int = 1999
    decision_only: bool = False
    alphas: tuple[float, ...] = (0.05, 0.01)
    inner_batch: int = 64

    def __post_init__(self) -> None:
        if self.method not in ("sc", "raw-count-refit"):
            raise ValueError("Unknown calibration method")
        if min(self.inner_draws, self.outer_draws, self.inner_batch) < 1:
            raise ValueError("MC budgets must be positive")
        if not self.alphas or any(not 0 < a < 1 for a in self.alphas):
            raise ValueError("alphas must lie in (0,1)")


def positive_ecdf(raw: Array, values: Array, outside: Array) -> Array:
    """Full receiver plus fixed outside; positive average ranks; zeros stay zero."""
    ids = np.searchsorted(values, raw)
    columns = raw.shape[1]
    counts = (
        np.bincount(
            (ids + np.arange(columns)[None, :] * len(values)).ravel(),
            minlength=len(values) * columns,
        )
        .reshape(columns, len(values))
        .T
        + outside[:, None]
    )
    positive = counts[1:]
    levels = np.vstack(
        [
            np.zeros(columns),
            (np.cumsum(positive, axis=0) - positive / 2 + 0.5)
            / np.maximum(positive.sum(axis=0), 1),
        ]
    )
    return levels[ids, np.arange(columns)[None, :]]


def _rng(seed: int, *identity: Any) -> np.random.Generator:
    digest = hashlib.sha256(json.dumps([seed, *identity]).encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:16], "little"))


def _score(c: Array, y: Array) -> Array:
    return np.asarray(np.einsum("i,ij->j", c, y, optimize=False))


def calibrate_count_reference(
    raw: Array,
    contrast: Array,
    values: Array,
    outside: Array,
    backend: CountBackend,
    *,
    seed: int,
    config: CountCalibrationConfig | None = None,
    fitted: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compute RC1 p or finite-MC bounds; SC is an explicit baseline option."""
    cfg = config or CountCalibrationConfig()
    raw, c, values, outside = map(np.asarray, (raw, contrast, values, outside))
    if raw.ndim != 1 or c.shape != raw.shape or len(raw) < 2:
        raise ValueError("Complete receiver and contrast must be matching vectors")
    if values.ndim != 1 or len(values) < 1 or values[0] != 0:
        raise ValueError("Count categories must start at zero")
    if not np.all(np.diff(values) > 0) or outside.shape != values.shape:
        raise ValueError("Invalid category ordering or outside histogram")
    if not all(np.isfinite(x).all() for x in (raw, c, values, outside)):
        raise ValueError("Nonfinite input")
    if any(np.any(x < 0) or np.any(x != np.rint(x)) for x in (raw, values, outside)):
        raise ValueError("Counts and histogram must be nonnegative integers")
    if not np.isin(raw, values).all() or abs(float(c.sum())) > 1e-9:
        raise ValueError("Unmatched count category or noncentered contrast")
    y = positive_ecdf(raw[:, None], values, outside)
    t = float(_score(c, y)[0])
    result: dict[str, Any] = {
        "method": cfg.method,
        "version": "RC1",
        "status": "FROZEN_MAIN" if cfg.method == "raw-count-refit" else "REFERENCE_BASELINE",
        "implementation_revision": "RC1_MAIN_20261013",
        "validation_status": "SCOPE_LIMITED",
        "score": t,
        "p": None,
        "p_interval": None,
        "reject": {},
        "inner_draws": cfg.inner_draws,
        "outer_budget": cfg.outer_draws,
        "outer_done": 0,
        "outer_hits": 0,
        "fits": 0,
        "inner_generated": 0,
    }
    if np.linalg.norm(c) <= 1e-12 or np.ptp(y) <= 1e-15:
        return {**result, "evaluable": False, "reason": "DEGENERATE_CONTRAST_OR_RESPONSE"}
    result["evaluable"] = True
    if t <= 1e-12:
        return {
            **result,
            "p": 1.0,
            "p_interval": [1.0, 1.0],
            "p_sc": 1.0,
            "reject": {str(a): False for a in cfg.alphas},
        }
    fit = fitted if fitted is not None else backend.fit(raw)
    e = _rng(seed, "observed_inner").normal(size=(cfg.inner_draws, len(raw))).T
    field = positive_ecdf(backend.bootstrap(fit, e), values, outside)
    observed_exceed = int(np.count_nonzero(_score(c, field) >= t))
    p_sc = (observed_exceed + 1) / (cfg.inner_draws + 1)
    result["p_sc"] = p_sc
    if cfg.method == "sc":
        return {
            **result,
            "p": p_sc,
            "p_interval": [p_sc, p_sc],
            "reject": {str(a): p_sc < a for a in cfg.alphas},
        }
    hits = 0
    for b in range(cfg.outer_draws):
        outer = backend.bootstrap(fit, _rng(seed, "outer", b).normal(size=(len(raw), 1)))
        outer_y = positive_ecdf(outer, values, outside)
        tb = float(_score(c, outer_y)[0])
        if np.ptp(outer_y) <= 1e-15 or tb <= 1e-12:
            qualifies = p_sc >= 1.0
        else:
            fb = backend.fit(outer[:, 0])
            result["fits"] += 1
            inner_rng = _rng(seed, "inner", b)
            count = done = 0
            # Inner p <= p_sc iff exceedance count <= observed_exceed.
            while done < cfg.inner_draws:
                size = min(cfg.inner_batch, cfg.inner_draws - done)
                eps = inner_rng.normal(size=(size, len(raw))).T
                yy = positive_ecdf(backend.bootstrap(fb, eps), values, outside)
                count += int(np.count_nonzero(_score(c, yy) >= tb))
                done += size
                if (
                    count > observed_exceed or count + cfg.inner_draws - done <= observed_exceed
                ):
                    break
            result["inner_generated"] += done
            qualifies = count <= observed_exceed
        hits += int(qualifies)
        completed = b + 1
        low = (hits + 1) / (cfg.outer_draws + 1)
        high = (hits + cfg.outer_draws - completed + 1) / (cfg.outer_draws + 1)
        resolved = all(high < a or low >= a for a in cfg.alphas)
        if cfg.decision_only and resolved:
            break
    result.update(outer_done=completed, outer_hits=hits, p_interval=[low, high])
    result["p"] = low if completed == cfg.outer_draws else None
    result["reject"] = {str(a): high < a for a in cfg.alphas}
    result["decision_exact_for_frozen_MC_budget"] = True
    return result
