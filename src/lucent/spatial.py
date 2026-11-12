"""Prepare ordered tissue hypotheses from full-section counts and geometry."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import rankdata
from .calibration import CountCalibrationConfig, calibrate_count_reference
from .calibration.copula import CountGeometry, GRID_FRACTIONS, GRID_RANGES
from .calibration.geometry import LazyGeometryModel
from .calibration.backend import CPUBackend


@dataclass
class SpatialData:
    """Whole-section raw counts and directed edges; rows are spots or cells."""
    counts: np.ndarray
    genes: np.ndarray
    coordinates: np.ndarray
    labels: np.ndarray
    sender: np.ndarray
    receiver: np.ndarray

    def __post_init__(self):
        self.counts = np.asarray(self.counts)
        self.genes = np.asarray(self.genes, dtype=str)
        self.coordinates = np.asarray(self.coordinates, dtype=float)
        self.labels = np.asarray(self.labels, dtype=str)
        self.sender = np.asarray(self.sender)
        self.receiver = np.asarray(self.receiver)
        n = len(self.coordinates)
        if self.coordinates.shape != (n, 2) or not np.isfinite(self.coordinates).all():
            raise ValueError('coordinates must be a finite n-by-2 array')
        if self.counts.shape != (n, len(self.genes)) or self.labels.shape != (n,):
            raise ValueError('counts, genes, coordinates and labels must align')
        if self.genes.ndim != 1 or len(set(self.genes)) != len(self.genes):
            raise ValueError('genes must be a unique one-dimensional array')
        if any(not x.strip() for x in self.genes) or any(not x.strip() for x in self.labels):
            raise ValueError('genes and labels must be nonempty')
        if (not np.issubdtype(self.counts.dtype, np.number) or
                not np.isfinite(self.counts).all() or np.any(self.counts < 0) or
                np.any(self.counts != np.rint(self.counts))):
            raise ValueError('counts must be finite nonnegative integer-valued measurements')
        if self.sender.ndim != 1 or self.receiver.shape != self.sender.shape:
            raise ValueError('sender and receiver must be aligned edge-index vectors')
        for index in (self.sender, self.receiver):
            if not np.issubdtype(index.dtype, np.integer) or np.any(index < 0) or np.any(index >= n):
                raise ValueError('edges must contain valid integer row indices')
        pairs = np.column_stack((self.sender, self.receiver))
        if np.any(self.sender == self.receiver) or len(np.unique(pairs, axis=0)) != len(pairs):
            raise ValueError('edges must be distinct directed non-self edges')

    @classmethod
    def load(cls, path: str | Path):
        with np.load(path, allow_pickle=False) as z:
            return cls(**{key: z[key] for key in
                         ('counts', 'genes', 'coordinates', 'labels', 'sender', 'receiver')})


def directed_knn(coordinates, k=8):
    """Build directed k-nearest-neighbour edges over a complete section."""
    xy = np.asarray(coordinates, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all():
        raise ValueError('coordinates must be finite n-by-2 values')
    if not isinstance(k, int) or not 1 <= k < len(xy):
        raise ValueError('k must be an integer in [1, n-1]')
    _, neighbours = cKDTree(xy).query(xy, k=k+1)
    receiver = np.concatenate([row[row != i][:k] for i, row in enumerate(neighbours)])
    return np.repeat(np.arange(len(xy)), k), receiver


def positive_ranks(counts):
    """Average positive ranks divided by the positive count; zeros stay zero."""
    counts = np.asarray(counts)
    result = np.zeros(len(counts), dtype=float)
    positive = counts > 0
    result[positive] = rankdata(counts[positive], method='average') / positive.sum()
    return result


def prepare_hypothesis(data, ligand, receptor, sender_label, receiver_label):
    """Construct edge contrasts while retaining the entire receiver label pool."""
    indices = {gene: i for i, gene in enumerate(data.genes)}
    if ligand not in indices or receptor not in indices:
        raise ValueError('both endpoint genes must occur in the expression matrix')
    if sender_label not in data.labels or receiver_label not in data.labels:
        raise ValueError('both labels must occur in the section')
    raw_l = data.counts[:, indices[ligand]]
    raw_r = data.counts[:, indices[receptor]]
    ligand_field, receptor_field = positive_ranks(raw_l), positive_ranks(raw_r)
    pool = np.flatnonzero(data.labels == receiver_label)
    mask = ((data.labels[data.sender] == sender_label) &
            (data.labels[data.receiver] == receiver_label))
    s, r = data.sender[mask], data.receiver[mask]
    contrast = np.zeros(len(pool))
    if len(s):
        u, v = ligand_field[s], receptor_field[r]
        contrast = np.bincount(np.searchsorted(pool, r), weights=(u-u.mean())/len(s), minlength=len(pool))
        product = float(np.mean(u*v))
        marginal_product = float(u.mean()*v.mean())
        score = float(np.einsum('i,i->', contrast, receptor_field[pool], optimize=False))
    else:
        product = marginal_product = score = None
    # Include a zero category even when every observed receptor count is positive.
    values = np.union1d([0], np.unique(raw_r))
    full_counts = np.bincount(np.searchsorted(values, raw_r), minlength=len(values))
    raw = raw_r[pool]
    outside = full_counts - np.bincount(np.searchsorted(values, raw), minlength=len(values))
    return dict(raw=raw, contrast=contrast, values=values, outside=outside,
                coordinates=data.coordinates[pool], pool_indices=pool,
                sender_indices=s, receiver_indices=r,
                summary=dict(ligand=ligand, receptor=receptor, sender=sender_label,
                             receiver=receiver_label, n_section=len(data.labels),
                             n_receiver=len(pool), n_edges=len(s),
                             n_joint_positive_edges=int(np.count_nonzero((raw_l[s]>0)&(raw_r[r]>0))),
                             product_score=product, marginal_product=marginal_product,
                             centered_score=score))


def test_hypothesis(data, ligand, receptor, sender_label, receiver_label, *, seed,
                    config=None, threads=4):
    """Evaluate one ordered hypothesis using bootstrap-refitted count calibration."""
    if not isinstance(seed, int) or seed < 0:
        raise ValueError('seed must be a nonnegative integer')
    if not isinstance(threads, int) or threads < 1:
        raise ValueError('threads must be a positive integer')
    prepared = prepare_hypothesis(data, ligand, receptor, sender_label, receiver_label)
    summary = prepared['summary']
    if not summary['n_edges'] or summary['n_receiver'] <= 8:
        return {'hypothesis': summary, 'seed': str(seed),
                'calibration': {'evaluable': False, 'score': summary['centered_score'], 'p': None,
                                'reason': 'NO_EDGES' if not summary['n_edges'] else 'RECEIVER_POOL_TOO_SMALL'}}
    xy = prepared['coordinates']
    model = LazyGeometryModel(xy, GRID_FRACTIONS, GRID_RANGES)
    backend = CPUBackend(CountGeometry(xy, model), threads)
    try:
        result = calibrate_count_reference(prepared['raw'], prepared['contrast'], prepared['values'],
                                           prepared['outside'], backend, seed=seed, config=config)
    finally:
        backend.close()
    return {'hypothesis': summary, 'seed': str(seed), 'calibration': result}
