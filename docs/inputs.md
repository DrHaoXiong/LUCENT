# Inputs and outputs

## Genetic inputs

All tables are tab-separated. Use consistent chromosome names, gene identifiers and genome assembly. Gene intervals are 1-based and closed.

| Table | Required columns | Optional columns |
|---|---|---|
| Qualified variants | `variant_id`, `locus_id`, `chrom`, `position`, `p_value` | `rsid`, `pip`, `neglog10p`, and other provenance columns |
| Gene coordinates | `gene_id`, `chrom`, `start`, `end` | `gene_symbol` |
| Direct channels | `channel_id`, `ligand`, `receptor` | `evidence` (default `Direct`) |

`variant_id` should encode chromosome, position and both alleles. Variant rows are the GWS variants retained after upstream QC and locus qualification; the `map` command checks GWS significance and performs mapping. LD qualification and raw GWAS harmonization are upstream inputs, rather than operations performed by this command.

When the source reports `p_value` as zero because the probability was not represented, include its original finite `neglog10p` (the source's −log₁₀(p)); it must pass the same genome-wide significance threshold. Mapping records retain `p_value_raw`, `neglog10p` and `p_value_status`, with `p_value` set to null rather than an exact zero. An explicitly supplied positive value below floating-point range is likewise retained as source text. A zero without a qualifying source `neglog10p` is rejected; no probability floor is imputed.

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

`lucent prepare-spatial` converts a whole QC-passed H5AD section and a previously built graph into this NPZ format. The graph NPZ must include string `unit_ids` in exactly the H5AD observation order and directed `sender`/`receiver` arrays (or `sender_index`/`receiver_index`). The converter checks their alignment, retains all observations and edges, and selects only the requested measured gene columns. It requires the optional `h5ad` dependency; ordinary NPZ use does not. See the [tutorial](tutorial.md#5-prepare-a-spatial-section-from-your-own-data).

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

### Tissue context and typography

The map retains all supplied spots or cells and their annotations. An optional registered tissue image is displayed beneath the spatial overlay and can be toggled on or off. Image bounds and observation coordinates must share the same coordinate system.

The interface and exported figures use Times New Roman, with Times and serif fallbacks. Gene symbols are italicized in candidate labels, expression legends and genetic evidence; identifiers, tissue labels and numbers remain upright. Install Times New Roman on the viewing computer for an exact font match.

### Association evidence

The default candidate view shows the coordination score and spatial association p value, with a `Nominal p < 0.05` filter. Supplied multiple-testing adjustments and the corresponding family size appear under **Statistical details**. JSON and TSV downloads retain these values; figure exports show the score and spatial p value.

### Figure export

Select a channel and section, set the color mode and zoom, then use the spatial toolbar:

- **SVG** downloads an editable vector figure with the current view, annotations and legend.
- **PDF** opens the same vector figure in a separate print view. Click **Save as PDF**, select the browser PDF destination, use 100% scale and disable headers and footers. The default figure width is 180 mm; height follows the view and legend. This path is tested in Chromium/Edge.
- **PNG** exports a raster image for quick sharing.

SVG and PDF use the same drawing commands as the on-screen map. Export does not rerun inference. Tissue-image pixels remain raster; points, edges, axes and text are vector elements. Exported SVG metadata retains the selected result, coordinate system and view settings. Existing report folders can be regenerated from their `downloads/bundle.json` with `lucent report` into a new output directory.
