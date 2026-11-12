# Practical tutorial: from genetic inputs to a spatial report

This walkthrough runs one synthetic ligand–receptor hypothesis through the command-line interface. It then shows how to prepare your own spatial section. Run commands from the repository root in a Unix-like shell; Windows users can use WSL. The example is synthetic and its small simulation budgets are for testing the workflow only.

## 1. Install and create example inputs

Download or clone the repository, enter its root directory, then run:

```bash
python -m pip install .
lucent --help
python examples/synthetic.py --output outputs/example
```

The example generates `variants.tsv`, `genes.tsv`, `catalogue.tsv`, `section.npz`, `hypothesis.json`, `result.json`, `report_bundle.json` and an offline report under `outputs/example/`. The script also runs the Python API; the steps below repeat mapping, testing and reporting with the CLI so you can inspect their hand-offs.

## 2. Map qualified variants to channels

```bash
lucent map --variants outputs/example/variants.tsv \
  --genes outputs/example/genes.tsv \
  --catalogue outputs/example/catalogue.tsv \
  --output outputs/example/cli_provenance.json
```

The JSON contains `mapping_records`, deduplicated `channels`, and endpoint-level `genetic_links`. The example nominates `demo_channel`, connecting `DEMO_L` to `DEMO_R`. Inspect them with:

```bash
python -c "import json; p=json.load(open('outputs/example/cli_provenance.json')); print(json.dumps(p['channels'], indent=2)); print(json.dumps(p['genetic_links'], indent=2))"
```

For your data, GWAS harmonization, QC, LD reference matching and locus qualification must be completed upstream. The mapper accepts the qualifying genome-wide significant variants and the optional PIP values you supply; it does not perform fine-mapping or LD qualification. Gene and endpoint IDs, genome assembly and chromosome naming must agree across tables. See [input formats](inputs.md) for columns and mapping rules. An empty `channels` list warrants checking ID and chromosome conventions before proceeding.

## 3. Test one ordered spatial hypothesis

Inspect `outputs/example/hypothesis.json`. Its ligand and receptor match the nominated channel; `sender_label` and `receiver_label` specify the ordered tissue context. Then run:

```bash
lucent test --input outputs/example/section.npz \
  --config outputs/example/hypothesis.json \
  --output outputs/example/cli_result.json
```

Inspect `hypothesis.centered_score`, `hypothesis.n_edges`, and `calibration.evaluable`, `calibration.p` and `calibration.reason` where present:

```bash
python -c "import json; r=json.load(open('outputs/example/cli_result.json')); print(json.dumps(r, indent=2))"
```

`p` is the calibrated one-sided spatial association probability for this hypothesis, not a PIP. A non-evaluable result has a reason and no spatial p value. The example uses 39 inner and 49 outer draws; for real analyses choose budgets appropriate to the desired Monte Carlo resolution and computational resources (the defaults are 999 and 1,999). Keep a seed for reproducibility.

`test` runs **one** ordered channel–context hypothesis at a time. To analyze more, construct a configuration for each eligible channel and ordered label pair, run `test` for each, and retain unavailable results as well as the chosen family and any subsequent multiple-testing adjustment. The CLI does not automatically loop through the channels in `cli_provenance.json` or calculate cross-sample consistency.

## 4. Connect the CLI outputs to a report

The report exporter accepts an assembled JSON bundle; it does not automatically combine a provenance file with test results. For this one-unit example, `examples/synthetic.py` has already prepared the sections and metadata in `report_bundle.json`. Save the following as `outputs/example/assemble_cli_bundle.py`, then run `python outputs/example/assemble_cli_bundle.py` from the repository root:

```python
import json
from pathlib import Path

root = Path('outputs/example')
bundle = json.loads((root / 'report_bundle.json').read_text(encoding='utf-8'))
provenance = json.loads((root / 'cli_provenance.json').read_text(encoding='utf-8'))
result = json.loads((root / 'cli_result.json').read_text(encoding='utf-8'))

# This demonstration has exactly one channel, one unit and one version.
channel = provenance['channels'][0]
unit = bundle['units'][0]
assert channel['channel_id'] == unit['channel_id']
assert (channel['ligand'], channel['receptor']) == (
    result['hypothesis']['ligand'], result['hypothesis']['receptor'])
assert (unit['sender_state'], unit['receiver_state']) == (
    result['hypothesis']['sender'], result['hypothesis']['receiver'])

cal = result['calibration']
unit['versions']['lucent'].update(
    score=cal['score'], p=cal['p'], q=None,
    status='complete' if cal['evaluable'] else 'not evaluable',
    evaluable=cal['evaluable'])
bundle['genetic_links'] = [
    {**link, 'gwas_p': link['p_value'],
     'pip_status': 'Provided PIP' if link.get('pip') is not None else 'Not provided',
     'build': 'synthetic'}
    for link in provenance['genetic_links']
    if link['channel_id'] == unit['channel_id']
]
(root / 'cli_bundle.json').write_text(json.dumps(bundle, indent=2) + '\n', encoding='utf-8')
```

Export and check the report:

```bash
lucent report --input outputs/example/cli_bundle.json --output outputs/example/cli_report
python -m lucent.report.verify outputs/example/cli_report
```

Open `outputs/example/cli_report/index.html`. The verification command checks the report bundle and output hashes; it does not rerun statistical inference or perform browser visual validation. The bundle's `donor_evidence` is empty in this example. Cross-sample evidence must be computed separately before it is added to a report.

For real data, create an explicit unit for each channel and ordered context, link genetic records to the correct channel and endpoint, and provide section coordinates, labels, edges and expression in the bundle. The assembly snippet above is specific to the single synthetic unit; it is not a general multi-sample report builder. See the `bundle` construction in [the example](../examples/synthetic.py) and the [bundle overview](inputs.md#reports) for the required structure.

## 5. Prepare a spatial section from your own data

Create a `section.npz` containing aligned arrays: raw integer `counts` (observations × genes), `genes`, `coordinates`, `labels`, and directed zero-based `sender`/`receiver` edge indices. Keep all QC-passed section observations, including zero counts, and build the graph before selecting a sender or receiver label. An annotated spot remains a spot; it should not automatically be interpreted as a single cell.

For example, if you already have aligned NumPy arrays in memory:

```python
import numpy as np
from lucent.spatial import directed_knn

# Supply counts, genes, coordinates and labels from your own prepared data.
# Row i of counts, coordinates and labels must describe the same observation;
# column j of counts must match genes[j]. Coordinates have shape (n, 2).
sender, receiver = directed_knn(coordinates, k=8)
np.savez_compressed('section.npz', counts=counts, genes=genes,
                    coordinates=coordinates, labels=labels,
                    sender=sender, receiver=receiver)
```

This is a preparation template, not a Visium/Xenium file importer. Convert platform files, apply QC and choose biologically suitable labels and coordinates upstream; alternatively supply your own directed graph. Save string arrays as ordinary NumPy string arrays, not object arrays, so the file can be loaded without pickle. Match `genes` to the catalogue endpoint IDs. Full requirements and the complete receiver-pool rule are in [inputs.md](inputs.md#spatial-input).

For each selected channel and label pair, create a JSON configuration modeled on the example, set the endpoint IDs and labels, and run `lucent test`. Before scaling up, check one hypothesis and estimate runtime and memory: the spatial reference uses dense receiver-pool operations, so cost depends on pool size as well as the number of hypotheses and simulation draws.
