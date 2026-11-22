"""Truncated Besag–Clifford sampling for groups of positive RC1 contrasts.

One call covers one receiver group. For full-family adjustment, retain all
original family rows; untested or nonpromoted rows receive adjustment input 1.
The Monte Carlo stopping rule does not establish exchangeability of the
underlying fitted parametric bootstrap.
"""
from dataclasses import dataclass

import numpy as np

from .core import observed, rounds, seed


@dataclass(frozen=True)
class Stopped:
    hits: int
    draws: int
    p: float


def stopped_p(hits, draws, maximum, target):
    if not (1 <= target <= maximum and 1 <= draws <= maximum and 0 <= hits <= min(target, draws)):
        raise ValueError("Invalid stopping counts")
    if hits == target:
        return hits / draws
    if draws == maximum:
        return (hits + 1) / (maximum + 1)
    raise ValueError("Neither stopping boundary reached")


def sample(n, event, maximum, target):
    """event(u, active) returns one inclusive outer hit per active contrast."""
    if n < 1 or not 1 <= target <= maximum:
        raise ValueError("Invalid group size or stopping budget")
    counts = np.zeros(n, dtype=int)
    results = [None] * n
    active = np.arange(n)
    for u in range(maximum):
        values = np.asarray(event(u, active))
        if values.shape != active.shape or not np.isin(values, [0, 1]).all():
            raise ValueError("Expected one binary outer outcome per active contrast")
        counts[active] += values.astype(int)
        done = u + 1
        stopped = (counts[active] == target) | (done == maximum)
        for j in active[stopped]:
            h = int(counts[j])
            results[j] = Stopped(h, done, stopped_p(h, done, maximum, target))
        active = active[~stopped]
        if len(active) == 0:
            break
    return results


def two_stage(raw, values, outside, coeff, ids, group_key, backend, fit, rc, score,
              shallow=(999, 1999), deep=(39999, 19999), target=20, cutoff=.02,
              *, round_sampler=rounds):
    """Return shallow and optional deep results in the supplied hypothesis order.

    Use a stable group_key identifying section, receptor and receiver pool.
    """
    ids = tuple(ids)
    coeff = np.asarray(coeff)
    if not ids or len(set(ids)) != len(ids) or len(ids) != len(coeff):
        raise ValueError("Expected distinct hypothesis IDs matching the contrast rows")
    if not 0 <= cutoff <= 1 or any(r < 1 or b < target or target < 1 for r, b in (shallow, deep)):
        raise ValueError("Invalid cutoff or sampling budgets")

    def stage(indices, budget, name):
        r, b = budget
        contrasts = coeff[indices]
        gs = seed("bc-group", name, *group_key)
        ks = np.array([observed(raw, values, outside, c, backend, fit,
                                seed("bc-observed", name, *group_key, ids[i]),
                                r, rc, score)[1] for i, c in zip(indices, contrasts)])

        def event(u, active):
            return round_sampler(raw, values, outside, contrasts[active], backend, fit, gs,
                                 ks[active], r, u, u + 1, rc, score)[0]

        return sample(len(indices), event, b, target)

    first = stage(np.arange(len(ids)), shallow, "screen")
    promoted = [i for i, result in enumerate(first) if result.p <= cutoff]
    second = dict(zip(promoted, stage(promoted, deep, "deep"))) if promoted else {}
    return [dict(id=h, shallow=first[i], deep=second.get(i),
                 family_input=second[i].p if i in second else 1.)
            for i, h in enumerate(ids)]
