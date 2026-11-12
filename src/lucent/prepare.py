"""Convert an annotated whole section to LUCENT's portable NPZ input."""
from pathlib import Path

import numpy as np

from .spatial import SpatialData


def prepare_h5ad_section(h5ad_path, genes, graph_path, output, *,
                         label_key, coordinates_key='spatial', layer=None):
    """Save selected measured genes while retaining every section observation.

    The graph must carry observation IDs in H5AD order and directed indices;
    it is reused, not reconstructed. H5AD input is opened read-only.
    """
    try:
        import anndata
    except ImportError as exc:
        raise ValueError('H5AD conversion requires anndata; from the repository root run: '
                         'python -m pip install ".[h5ad]"') from exc
    genes = [str(g).strip() for g in genes]
    if not genes or any(not g for g in genes) or len(set(genes)) != len(genes):
        raise ValueError('gene list must contain unique, nonempty identifiers')
    h5ad_path, graph_path, output = map(Path, (h5ad_path, graph_path, output))
    if output.suffix.lower() != '.npz':
        raise ValueError('spatial output must have the .npz extension')
    if output.resolve() in (h5ad_path.resolve(), graph_path.resolve()):
        raise ValueError('output must differ from both input files')
    if output.exists():
        raise ValueError('output already exists; select a new path')
    a = anndata.read_h5ad(h5ad_path, backed='r')
    try:
        if not a.obs_names.is_unique or not a.var_names.is_unique:
            raise ValueError('H5AD observation and gene identifiers must be unique')
        unit_ids = a.obs_names.to_numpy(dtype=str)
        if any(not x.strip() for x in unit_ids):
            raise ValueError('observation identifiers must be nonempty')
        missing = [g for g in genes if g not in a.var_names]
        if missing:
            raise ValueError('requested genes absent from H5AD: ' + ', '.join(missing[:10]))
        if label_key not in a.obs or a.obs[label_key].isna().any():
            raise ValueError('label column must exist and contain no missing values')
        if coordinates_key not in a.obsm:
            raise ValueError('coordinate key is absent from H5AD obsm')
        if layer is not None and layer not in a.layers:
            raise ValueError('requested raw-count layer is absent from H5AD')
        with np.load(graph_path, allow_pickle=False) as graph:
            if 'unit_ids' not in graph:
                raise ValueError('graph must include unit_ids in H5AD order')
            graph_ids = np.asarray(graph['unit_ids'], dtype=str)
            if not np.array_equal(unit_ids, graph_ids):
                raise ValueError('graph unit_ids do not match H5AD observation order')
            if 'sender' in graph and 'receiver' in graph:
                sender, receiver = graph['sender'], graph['receiver']
            elif 'sender_index' in graph and 'receiver_index' in graph:
                sender, receiver = graph['sender_index'], graph['receiver_index']
            else:
                raise ValueError('graph requires sender/receiver indices')
        subset = a[:, genes]
        counts = subset.X if layer is None else subset.layers[layer]
        counts = counts.toarray() if hasattr(counts, 'toarray') else np.asarray(counts)
        data = SpatialData(counts, np.asarray(genes, dtype=str),
                           np.asarray(a.obsm[coordinates_key]),
                           a.obs[label_key].to_numpy(dtype=str), sender, receiver)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('xb') as handle:
            np.savez_compressed(handle, counts=data.counts, genes=data.genes,
                                coordinates=data.coordinates, labels=data.labels,
                                sender=data.sender, receiver=data.receiver,
                                unit_ids=unit_ids)
    finally:
        a.file.close()
    return output
