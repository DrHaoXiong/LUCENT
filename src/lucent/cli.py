"""Command-line entry points for prepared genetic and spatial inputs."""
import argparse
import csv
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
    mapping = commands.add_parser('map', help='Construct a channel provenance ledger')
    for name in ('variants', 'genes', 'catalogue', 'output'):
        mapping.add_argument('--'+name, required=True, type=Path)
    mapping.add_argument('--window-bp', type=int, default=50000)
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
        elif args.command == 'map':
            mapped = core_mapping(read_table(args.variants), read_table(args.genes), window_bp=args.window_bp)
            channels, ledger = nominate_channels(mapped, read_table(args.catalogue))
            write_json(args.output, {'channels': channels, 'genetic_links': ledger, 'mapping_records': mapped})
        else:
            from .report.exporter import main as export
            return export(['--input',str(args.input),'--output',str(args.output)])
    except (ValueError, KeyError, OSError) as exc:
        print(f'lucent: {exc}', file=sys.stderr)
        return 2
    print(args.output)
    return 0
