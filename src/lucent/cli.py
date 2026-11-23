"""Command-line entry points for prepared genetic and spatial inputs."""
import argparse
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
from . import __version__, SpatialData, test_hypothesis, CountCalibrationConfig
from .provenance import core_mapping, nominate_channels


def read_table(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')


def main(argv=None):
    parser = argparse.ArgumentParser(description='LUCENT: locus-traceable spatial communication evidence')
    parser.add_argument('--version', action='version', version=__version__)
    commands = parser.add_subparsers(dest='command', required=True)
    test = commands.add_parser('test', help='Calibrate one ordered hypothesis in a full section')
    test.add_argument('--input', required=True, type=Path, help='Full-section NPZ')
    test.add_argument('--config', required=True, type=Path, help='Hypothesis JSON')
    test.add_argument('--output', required=True, type=Path)
    batch = commands.add_parser('batch', help='Run grouped two-stage calibration on one section')
    batch.add_argument('--input', required=True, type=Path, help='Full-section NPZ')
    batch.add_argument('--config', required=True, type=Path, help='Positive hypotheses and budgets JSON')
    batch.add_argument('--output', required=True, type=Path)
    batch.add_argument('--workers', type=int, default=1, help='Outer worker processes')
    batch.add_argument('--large-pool-blas-threads', type=int, default=1,
                       help='BLAS projection threads per worker for pools above 768 points')
    batch.add_argument('--block-size', type=int, default=128, help='Deep outer draws per block')
    mapping = commands.add_parser('map', help='Construct a channel provenance ledger')
    for name in ('variants', 'genes', 'catalogue', 'output'):
        mapping.add_argument('--'+name, required=True, type=Path)
    mapping.add_argument('--window-bp', type=int, default=50000)
    prepare = commands.add_parser('prepare-spatial', help='Convert a whole H5AD section and saved graph to NPZ')
    prepare.add_argument('--h5ad', required=True, type=Path)
    prepare.add_argument('--genes-file', required=True, type=Path, help='One measured gene ID per line, no header')
    prepare.add_argument('--graph', required=True, type=Path, help='Graph NPZ including unit_ids')
    prepare.add_argument('--label-key', required=True, help='H5AD obs label column')
    prepare.add_argument('--coordinates-key', default='spatial', help='H5AD obsm key (default: spatial)')
    prepare.add_argument('--layer', help='H5AD raw-count layer (default: X)')
    prepare.add_argument('--output', required=True, type=Path)
    report = commands.add_parser('report', help='Export an offline interactive report')
    report.add_argument('--input', required=True, type=Path)
    report.add_argument('--output', required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == 'test':
            config = json.loads(args.config.read_text(encoding='utf-8-sig'))
            allowed = {'ligand','receptor','sender_label','receiver_label','seed','threads','inner_draws','outer_draws'}
            if set(config)-allowed:
                raise ValueError('unknown config keys: '+', '.join(sorted(set(config)-allowed)))
            calibration = CountCalibrationConfig(inner_draws=config.pop('inner_draws',999),
                                                  outer_draws=config.pop('outer_draws',1999))
            result = test_hypothesis(SpatialData.load(args.input), config=calibration, **config)
            result['software_version'] = __version__
            result['input_sha256'] = hashlib.sha256(args.input.read_bytes()).hexdigest()
            result['config_sha256'] = hashlib.sha256(args.config.read_bytes()).hexdigest()
            write_json(args.output, result)
        elif args.command == 'batch':
            from .calibration.optimized.batch import calibrate_batch
            config_bytes = args.config.read_bytes()
            cfg = json.loads(config_bytes)
            section = cfg.pop('section')
            groups = {}
            for h in cfg.pop('hypotheses'):
                key = (h['receptor'], h['receiver_label'])
                hypotheses = groups.setdefault(key, {})
                if h['id'] in hypotheses:
                    raise ValueError(f"duplicate hypothesis ID {h['id']}")
                hypotheses[h['id']] = (h['ligand'], h['sender_label'])
            if set(cfg) - {'shallow', 'deep', 'target', 'cutoff'}:
                raise ValueError('unknown batch calibration settings')
            with args.input.open('rb') as handle:
                input_sha = hashlib.file_digest(handle, 'sha256').hexdigest()
            rows = calibrate_batch(args.input, section, groups, workers=args.workers,
                                   large_pool_blas_threads=args.large_pool_blas_threads,
                                   block_size=args.block_size, **cfg)
            records = [dict(row, shallow=asdict(row['shallow']),
                            deep=asdict(row['deep']) if row['deep'] else None) for row in rows]
            write_json(args.output, dict(section=section, rows=records, workers=args.workers,
                       large_pool_blas_threads=args.large_pool_blas_threads, block_size=args.block_size,
                       software_version=__version__, input_sha256=input_sha,
                       config_sha256=hashlib.sha256(config_bytes).hexdigest()))
        elif args.command == 'map':
            mapped = core_mapping(read_table(args.variants), read_table(args.genes), window_bp=args.window_bp)
            channels, ledger = nominate_channels(mapped, read_table(args.catalogue))
            write_json(args.output, {'channels': channels, 'genetic_links': ledger, 'mapping_records': mapped})
        elif args.command == 'prepare-spatial':
            from .prepare import prepare_h5ad_section
            genes = [line.strip() for line in args.genes_file.read_text(encoding='utf-8-sig').splitlines()
                     if line.strip() and not line.lstrip().startswith('#')]
            prepare_h5ad_section(args.h5ad, genes, args.graph, args.output,
                                 label_key=args.label_key, coordinates_key=args.coordinates_key,
                                 layer=args.layer)
        else:
            from .report.exporter import main as export
            return export(['--input',str(args.input),'--output',str(args.output)])
    except (ValueError, KeyError, OSError) as exc:
        print(f'lucent: {exc}', file=sys.stderr)
        return 2
    print(args.output)
    return 0
