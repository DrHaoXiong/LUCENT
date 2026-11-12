# Inputs and outputs

## Genetic inputs

All tables are tab-separated. Use consistent chromosome names, gene identifiers and genome assembly. Gene intervals are 1-based and closed.

| Table | Required columns | Optional columns |
|---|---|---|
| Qualified variants | `variant_id`, `locus_id`, `chrom`, `position`, `p_value` | `rsid`, `pip`, and other provenance columns |
| Gene coordinates | `gene_id`, `chrom`, `start`, `end` | `gene_symbol` |
| Direct channels | `channel_id`, `ligand`, `receptor` | `evidence` (default `Direct`) |

`variant_id` should encode chromosome, position and both alleles. Variant rows are the GWS variants retained after upstream QC and locus qualification; the `map` command checks GWS significance and performs mapping. LD qualification and raw GWAS harmonization are upstream inputs, rather than operations performed by this command.

Core mapping retains gene-body overlaps and links within the configured distance from a gene interval (50 kb by default). If a locus has no such links, nearest-gene mapping is applied to its eligible variants, retaining ties. Channels enter when at least one endpoint has a link. The output preserves all variant, locus and endpoint records, with a deduplicated channel list. These commands use single-gene ligand and receptor endpoints.

`provenance.json` contains `mapping_records`, `channels` and `genetic_links`. Missing PIP values remain null. Additional input columns are retained in the mapping records.

## Spatial input

`section.npz` contains these arrays, loadable without pickle:

| Array | Shape | Meaning |
|---|---|---|
| `counts` | n × g | Nonnegative integer-valued raw counts |
| `genes` | g | Unique string identifiers matching the channel endpoints |
| `coordinates` | n × 2 | Spatial coordinates in a common unit |
| `labels` | n | Region or cell-type labels |
| `sender`, `receiver` | e each | Zero-based row indices of directed neighbour edges |

Keep all QC-passed spots/cells of the section, including zero counts. Gene columns may be restricted to the channels of interest. Build the graph over the full section before selecting ordered labels. `lucent.spatial.directed_knn` provides a k-nearest-neighbour graph (default k=8); saved graphs can also be supplied directly.

Example `hypothesis.json`:

```json
{
  "ligand": "GENE_L",
  "receptor": "GENE_R",
  "sender_label": "Stroma",
  "receiver_label": "Tumour",
  "seed": 1234,
  "threads": 4,
  "inner_draws": 999,
  "outer_draws": 1999
}
```

The receiver pool includes every spot/cell carrying the receiver label, including those with no selected edge. The geometry requires at least nine receiver observations. Hypotheses with no edges, insufficient receiver observations or a degenerate expression contrast are returned as unavailable with a reason. For evaluable nonpositive scores, the one-sided probability is 1.

`threads` controls probability integration workers; dense linear algebra uses the frozen backend's own thread limits. Covariance storage grows quadratically with receiver-pool size. Results include the raw product, marginal product, centered score, receiver and edge counts, seed, calibration telemetry and input/configuration SHA-256 hashes.

## Python API

```python
from lucent import SpatialData, CountCalibrationConfig, test_hypothesis

data = SpatialData.load("section.npz")
result = test_hypothesis(
    data, "GENE_L", "GENE_R", "Stroma", "Tumour", seed=1234,
    config=CountCalibrationConfig(inner_draws=999, outer_draws=1999),
    threads=4,
)
```

The lower-level `lucent.calibration.calibrate_count_reference` accepts a complete receiver count vector, centered contrast, full-section count categories, outside-pool histogram and count backend. It retains the paper implementation's numerical metadata, including the internal implementation identifier `RC1`.

## Reports

The report bundle has five sections: `metadata`, `units` (channel–context results), `genetic_links`, `donor_evidence` and `sections` (coordinates, labels, edges and expression). `examples/synthetic.py` constructs a complete minimal bundle. The exporter validates cross-references and preserves supplied association and cross-sample results. Exported reports work offline and include downloadable JSON/TSV tables and an output hash manifest.
