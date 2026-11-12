"""Core positional mapping and endpoint-resolved genetic provenance."""
import math


def core_mapping(variants, genes, *, window_bp=50_000):
    """Map QC/LD-qualified GWS variants using body/window and locus-level fallback.

    Coordinates are 1-based, closed, and on a common genome assembly.
    Input variants retain allele-aware IDs and locus qualification from upstream.
    """
    if not isinstance(window_bp, int) or window_bp < 0:
        raise ValueError('window_bp must be a nonnegative integer')
    by_chrom = {}
    gene_ids = set()
    for gene in genes:
        gene = dict(gene)
        gene['start'], gene['end'] = int(gene['start']), int(gene['end'])
        if gene['start'] < 1 or gene['end'] < gene['start'] or gene['gene_id'] in gene_ids:
            raise ValueError('invalid or duplicate gene interval')
        gene_ids.add(gene['gene_id'])
        by_chrom.setdefault(str(gene['chrom']), []).append(gene)
    groups, seen = {}, set()
    for variant in variants:
        variant = dict(variant)
        if variant['variant_id'] in seen:
            raise ValueError('duplicate variant identity')
        seen.add(variant['variant_id'])
        variant['position'] = int(variant['position'])
        variant['p_value'] = float(variant['p_value'])
        if variant['position'] < 1 or not 0 < variant['p_value'] < 5e-8:
            raise ValueError('mapping requires GWS variants at valid positions')
        if variant.get('pip') in ('', None):
            variant['pip'] = None
        else:
            variant['pip'] = float(variant['pip'])
            if not math.isfinite(variant['pip']) or not 0 <= variant['pip'] <= 1:
                raise ValueError('PIP must be missing or in [0,1]')
        groups.setdefault(variant['locus_id'], []).append(variant)
    output = []
    for locus_id, members in groups.items():
        if len({str(v['chrom']) for v in members}) != 1:
            raise ValueError('a locus must lie on one chromosome')
        candidates = []
        for variant in members:
            position = variant['position']
            for gene in by_chrom.get(str(variant['chrom']), []):
                distance = max(gene['start']-position, position-gene['end'], 0)
                candidates.append((distance, variant, gene))
        selected = [x for x in candidates if x[0] <= window_bp]
        fallback = not selected
        if fallback:
            # Nearest gene(s) for each eligible variant in a locus with no window link.
            for variant in members:
                options = [x for x in candidates if x[1]['variant_id'] == variant['variant_id']]
                if options:
                    minimum = min(x[0] for x in options)
                    selected.extend(x for x in options if x[0] == minimum)
        for distance, variant, gene in selected:
            output.append({**variant, 'gene_id': gene['gene_id'],
                           'gene_symbol': gene.get('gene_symbol', gene['gene_id']),
                           'distance_bp': distance,
                           'mapping_basis': 'nearest_gene' if fallback else
                           ('gene_body' if distance == 0 else 'positional_window')})
    return output


def nominate_channels(mapping, catalogue):
    """Join direct single-gene channels to every supporting mapping record."""
    by_gene = {}
    for record in mapping:
        by_gene.setdefault(record['gene_id'], []).append(record)
    channels, ledger, seen = [], [], set()
    for channel in catalogue:
        channel = dict(channel)
        if channel['channel_id'] in seen:
            raise ValueError('catalogue channel IDs must be unique')
        seen.add(channel['channel_id'])
        if channel.get('evidence', 'Direct') != 'Direct':
            continue
        links = []
        for endpoint in ('ligand', 'receptor'):
            links.extend({**row, 'channel_id': channel['channel_id'], 'endpoint_role': endpoint}
                         for row in by_gene.get(channel[endpoint], []))
        if links:
            channels.append(channel)
            ledger.extend(links)
    return channels, ledger
