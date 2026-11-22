# LUCENT

**Locus-traceable and Uncertainty-calibrated Cell–cell communication Evidence in Native Tissue**

LUCENT connects trait-associated loci to ligand–receptor channels and ordered tissue interfaces. It preserves genetic provenance, measures centered spatial coordination, and calibrates association evidence with a fitted count-and-spatial reference.

## Install

Python 3.11 or later:

```bash
python -m pip install .
lucent --help
```

`requirements-tested.txt` records the dependency versions used to verify the frozen numerical example.

## Run the example

```bash
python examples/synthetic.py --output outputs/example
```

Open `outputs/example/report/index.html` to explore the spatial map, association result and variant-to-endpoint links. The example generates its own synthetic inputs and uses small bootstrap budgets for a quick run.

The offline report combines the supplied tissue map and annotations, candidate channels, genetic provenance and cross-sample evidence. It uses a white background, Times New Roman and italic gene symbols.

The spatial toolbar exports **SVG**, **PDF** and **PNG**. SVG retains editable vectors. PDF opens a figure-only view: choose **Save as PDF**, select the browser's PDF destination and disable headers and footers. Exports retain the current zoom, channel, color mode and legends. Points, edges and text remain vectors; an optional tissue image retains its original raster representation.

## Use your data

1. Prepare GWAS variants that pass the chosen QC and locus/LD qualification rules, gene coordinates on the same assembly, and a direct ligand–receptor catalogue.
2. Construct the genetic ledger:

```bash
lucent map --variants qualified_variants.tsv --genes genes.tsv \
  --catalogue channels.tsv --output outputs/provenance.json
```

3. Supply whole-section counts, coordinates, annotations and directed neighbour edges. Calibrate an ordered hypothesis:

```bash
lucent test --input section.npz --config hypothesis.json --output outputs/result.json
```

4. Export an assembled evidence bundle as an offline interactive report:

```bash
lucent report --input report_bundle.json --output outputs/report
```

Follow the [practical tutorial](docs/tutorial.md) for an end-to-end example and spatial input preparation. See [input formats](docs/inputs.md) for the complete input requirements. The Python API exposes the same mapping and calibration routines.

For an annotated H5AD section and an existing whole-section graph, install `.[h5ad]` and use `lucent prepare-spatial` to create `section.npz`. Prepared NPZ inputs do not require this optional dependency.

## Statistical workflow

The coordination score subtracts the product of edge-weighted marginal means from the mean ligand–receptor product. Calibration holds the ligand field and spatial geometry fixed, fits the receiver count distribution and spatial dependence, and refits the reference in outer bootstrap fields. Expression ranks are calculated over the whole section, retaining zero counts and the complete receiver-label pool.

`inner_draws` and `outer_draws` are configurable; defaults are 999 and 1,999. Each result records its seed, score, calibrated probability and simulation counts. Optional fine-mapping PIPs remain annotations on the genetic ledger. Cross-sample evidence can be supplied to the report alongside the association results.

## Repository contents

- `src/lucent/calibration`: count reference, spatial covariance and refit calibration.
- `src/lucent/spatial.py`: full-section input and ordered hypothesis preparation.
- `src/lucent/provenance.py`: positional Core mapping and endpoint provenance.
- `src/lucent/report`: offline interactive report and downloadable tables.
- `examples`: small synthetic workflow.
- `tests`: numerical and data-contract regression checks.

Paper-specific benchmarks, figure scripts and research datasets belong to the accompanying reproducibility archive. This repository contains the reusable analysis software.

## Citation

The LUCENT manuscript is in preparation. A link and full citation for the preprint or published article will be added here when available. `CITATION.cff` currently uses a collective project name until the author list is available.

## License

LUCENT software is released under the [MIT License](LICENSE). Third-party inputs retain their own terms; see [THIRD_PARTY.md](THIRD_PARTY.md).

## Tests

```bash
python -m unittest discover -s tests -v
```

## Grouped two-stage calibration

The optional Python API builds the four optimized RC1 kernels and routines on x86-64 Linux or WSL2. Install GCC with OpenMP and glibc development libraries on an AVX2-capable CPU. The standard CLI uses the fixed-budget reference implementation.

This example takes a whole-section NPZ in the [input format](docs/inputs.md). Replace the example names with positive hypotheses from your section. Open one session per worker and reuse it across receiver groups.

```python
from lucent.spatial import SpatialData
from lucent.calibration.optimized import optimized_session, adjust_family

data = SpatialData.load("section.npz")
with optimized_session() as session:
    rows = session.calibrate_group(
        data, section="section1", receptor="ITGB1", receiver_label="tumour",
        hypotheses={"h1": ("FN1", "stroma"), "h2": ("COL4A1", "stroma")},
        shallow=(999, 1999), deep=(39999, 19999), target=20, cutoff=0.02,
    )
```

Each budget gives the inner and maximum outer draws. Hypotheses sharing a receptor and receiving region reuse simulated fields and outer refits. The native kernels accelerate likelihood-grid evaluation, the 256-node check, discretization and positive-rank scoring. Each row records shallow and, on promotion, independently seeded deep tail counts, draws and sequential p-value. To compute BH and BY across a section, collect rows for every eligible positive hypothesis and call `adjust_family(full_family_ids, eligible_positive_ids, rows)` with the complete original family. Nonpromoted and ineligible entries contribute 1 to adjustment.
