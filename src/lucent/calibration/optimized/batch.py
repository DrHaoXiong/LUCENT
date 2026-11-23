"""Process-parallel calibration of positive hypotheses on one complete section."""
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
from multiprocessing.util import Finalize

import numpy as np

from ...spatial import SpatialData
from .. import refit
from .core import observed, seed
from .sequential import sample


def _initialize(path, section, groups, threads):
    from . import optimized_session
    global _data, _section, _groups, _threads, _session, _contexts, _key, _prepared
    _data = SpatialData.load(path)
    _section, _groups, _threads = section, groups, threads
    resources = ExitStack()
    Finalize(None, resources.close, exitpriority=10)
    _session = resources.enter_context(optimized_session())
    _contexts = resources.enter_context(ExitStack())
    _key = None


def _group(key):
    global _key, _prepared
    if key != _key:
        _contexts.close()
        _key = None
        _prepared = _contexts.enter_context(_session._prepare(_data, *key, _groups[key], _threads))
        _key = key
    return _prepared


def _cutoffs(key, indices, r, name):
    raw, v, o, coeff, be, fit, fields = _group(key)
    ids = list(_groups[key])
    return np.array([observed(raw, v, o, coeff[i], be, fit,
                              seed("bc-observed", name, _section, *key, ids[i]),
                              r, refit, fields.score)[1] for i in indices])


def _events(key, indices, k, r, name):
    raw, v, o, coeff, be, fit, fields = _group(key)
    contrasts = coeff[indices]
    gs = seed("bc-group", name, _section, *key)

    def event(u, active):
        return fields.rounds(raw, v, o, contrasts[active], be, fit, gs, k[active],
                             r, u, u+1, refit, fields.score)[0]
    return event


def _screen(key, shallow, deep, target, cutoff, block_size):
    indices = np.arange(len(_groups[key]))
    k = _cutoffs(key, indices, shallow[0], "screen")
    first = sample(len(indices), _events(key, indices, k, shallow[0], "screen"), shallow[1], target)
    promoted = np.array([i for i, row in enumerate(first) if row.p <= cutoff], dtype=int)
    deep_k = _cutoffs(key, promoted, deep[0], "deep")
    event = _events(key, promoted, deep_k, deep[0], "deep")
    # Finish short deep trajectories in the process that already holds the group.
    hits = np.zeros(len(promoted), dtype=int)
    prefix = []
    for u in range(min(block_size, deep[1])):
        active = np.flatnonzero(hits < target)
        if not len(active):
            break
        values = event(u, active)
        prefix.append(values)
        hits[active] += values
    return first, promoted, deep_k, prefix


def _block(key, indices, k, r, begin, end):
    event = _events(key, indices, k, r, "deep")
    active = np.arange(len(k))
    return np.asarray([event(u, active) for u in range(begin, end)], dtype=np.uint8)


def _deep_results(pool, key, promoted, k, prefix, deep, target, workers, block_size):
    # Consume in outer-index order. A retirement can change scalar vs. grouped
    # floating-point scoring, so invalidate all lookahead made for the old set.
    batch = np.empty((0, 0), dtype=np.uint8)
    start = 0
    members = None

    def event(u, active):
        nonlocal batch, start, members
        if u < len(prefix):
            return prefix[u]
        if u >= start + len(batch) or not np.array_equal(active, members):
            start, members = u, active.copy()
            stop = min(deep[1], u + workers * block_size)
            futures = [pool.submit(_block, key, promoted[active], k[active], deep[0], a,
                                   min(a + block_size, stop)) for a in range(u, stop, block_size)]
            batch = np.concatenate([f.result() for f in futures])
        return batch[u-start]

    return sample(len(promoted), event, deep[1], target)


def calibrate_batch(section_path, section, groups, *, workers=1, large_pool_blas_threads=1,
                    block_size=128, shallow=(999, 1999), deep=(39999, 19999),
                    target=20, cutoff=.02):
    """Return ordered two-stage rows for {(receptor, region): {id: (ligand, sender)}}.

    Workers load the same full-section NPZ. Positive hypothesis IDs must be
    unique across this section. block_size counts outer draws per work unit.
    """
    from . import optimized_session
    if not isinstance(section, str) or not section.strip():
        raise ValueError("Supply a nonempty section name")
    if any(not isinstance(budget, (tuple, list)) or len(budget) != 2 for budget in (shallow, deep)):
        raise ValueError("Supply an inner and an outer draw budget for each stage")
    sizes = (workers, large_pool_blas_threads, block_size, target, *shallow, *deep)
    if any(type(n) is not int or n < 1 for n in sizes):
        raise ValueError("Budgets, block size and thread/process counts must be positive integers")
    if target > min(shallow[1], deep[1]) or not 0 <= cutoff <= 1:
        raise ValueError("Invalid sequential calibration settings")
    ids = [i for hypotheses in groups.values() for i in hypotheses]
    if any(not h for h in groups.values()) or len(set(ids)) != len(ids):
        raise ValueError("Supply nonempty receiver groups and distinct section-wide hypothesis IDs")
    if not ids:
        return []
    if workers == 1:
        data = SpatialData.load(section_path)
        with optimized_session() as session:
            return [row for key, hs in groups.items() for row in session.calibrate_group(
                data, section, *key, hs, large_pool_blas_threads=large_pool_blas_threads,
                shallow=shallow, deep=deep, target=target, cutoff=cutoff)]
    with ProcessPoolExecutor(workers, mp_context=mp.get_context("spawn"), initializer=_initialize,
                             initargs=(section_path, section, groups, large_pool_blas_threads)) as pool:
        tasks = {key: pool.submit(_screen, key, shallow, deep, target, cutoff, block_size)
                 for key in groups}
        rows = []
        for key, task in tasks.items():
            first, promoted, k, prefix = task.result()
            second = {}
            if len(promoted):
                results = _deep_results(pool, key, promoted, k, prefix, deep, target, workers, block_size)
                second = dict(zip(promoted, results, strict=True))
            for i, identity in enumerate(groups[key]):
                result = second.get(i)
                rows.append(dict(id=identity, shallow=first[i], deep=result,
                                 family_input=result.p if result else 1.))
    return rows
