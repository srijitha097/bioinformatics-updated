> **Start here on macOS:** read `../START_HERE_MAC.md` if using the combined ZIP. Keep your original demo folder backed up before changing files.

# Original Research Project 02 — Human DLPFC snRNA-seq × Spatial Transcriptomics

**Portfolio status:** complete scripts and **locally tested synthetic Python demo**. The real LIBD download/conversion script requires R/Bioconductor, internet access, and has **not** been executed in the supplied bundle.

## Original research objective

Can independent human DLPFC single-nucleus cell-type expression signatures predict the location and relative prevalence of neuronal/glial populations across spatial transcriptomic tissue sections? How reproducible are those distributions among brain sections, and how well can the signatures reconstruct **held-out genes**?

This is an original computational question on existing open data, not a reproduction of the original LIBD publication figures or conclusions.

## Data — **directly relevant to LIBD staff scientist role**

Official project: https://github.com/LieberInstitute/spatialDLPFC

Data access: https://research.libd.org/spatialLIBD/reference/fetch_data.html

- A processed annotated `SingleCellExperiment` contains roughly **77,604 nuclei** across **19 snRNA-seq samples**.
- Spatial transcriptomics: `spatialDLPFC_Visium_example_subset` provides a **three-section Visium subset**, around **12,107 spots** before our memory-limited sampling. The larger 30-sample Visium dataset is also available.
- We **derive our own per-cell-type reference signatures from raw single-nucleus UMI counts**, align gene identifiers to Visium, and run a simple, interpretable non-negative least squares (NNLS) deconvolution baseline.
- These measurements are from **neurotypical** donors. Do NOT describe them as schizophrenia case-control or claim direct disease differential expression.

## 1 — Quick run with synthetic data (Python only)

```bash
cd 02_dlpfc_spatial_integration
conda env create -f environment.yml
conda activate libd-spatial
bash run_demo.sh
python -m pytest -q tests
```

This simulates 972 spatial spots and five cell populations. It is purely a **software smoke test**; the preview figures and metrics are not biological findings.

## 2 — Real analysis with LIBD human DLPFC data

First install Python dependencies (`environment.yml`) and R/Bioconductor dependencies:

```r
install.packages('BiocManager')
BiocManager::install(c(
 'spatialLIBD', 'HDF5Array','SpatialExperiment',
 'SummarizedExperiment','SingleCellExperiment'
), ask=FALSE, update=FALSE)
```

Then:

```bash
bash run_real.sh
```

The R downloader uses `spatialLIBD::fetch_data('spatialDLPFC_snRNAseq')` and `spatialLIBD::fetch_data('spatialDLPFC_Visium_example_subset')`. It finds the HDF5-backed nucleus object, chooses common gene symbols, exports cell-type reference signatures and a compact raw Visium count matrix. Full source datasets remain in `data/raw/` or Bioconductor cache (ignored by Git). The R script will show metadata columns and stop clearly if LIBD changes its annotation schema.

> Resource tip: the full snRNA-seq data are HDF5-backed; the example Visium subset and downsampling keep memory lower, but downloads can still be sizable. Prefer 16–32 GB RAM if possible; expand via HPC when ready. Data transfer requires internet access.

## Script breakdown

| Stage | Script | Purpose |
|---|---|---|
| Synthetic demo | `00_make_demo.py` | Generate 3 simulated sections and 5 cell-type signatures |
| Real import | `01_fetch_libd.R` | Fetch HDF5 snRNAseq + Visium subset, match gene symbols, construct signatures |
| Cell type assignment | `02_deconvolve.py` | Fit nonnegative weighted reference signatures to each Visium spot via NNLS |
| Spatial organization | `03_spatial_stats.py` | Per-section kNN Moran's I, coordinate permutations, FDR, layer summaries if available |
| Robustness | `04_gene_holdout_validation.py` | Withhold 20% of genes from fitting, evaluate held-out expression reconstruction |

## Outputs after run

```text
results/
  spot_celltype_proportions.csv
  moran_by_section.csv
  layer_celltype_estimates.csv       (if layers available)
  heldout_gene_validation.csv
  figures/01_mean_composition.png
  figures/02_spatial_<celltype>.png
  figures/03_layer_celltype_composition.png (if layers available)
```

### Scientific interpretation / caution

- A Visium spot is a mixture, not one cell. NNLS fractions are **estimates**, not observed fractions or proof of cell-cell contact.
- Moran's I quantifies **spatial autocorrelation**, not molecular interaction. Per-section permutation tests break the spatial correspondence, not donor-level confounding.
- Hold-out gene reconstruction measures predictive consistency with reference expression; it is not validation against ground-truth cell counts. Follow up with orthogonal histology/ISH if available.
- If multiple sections come from the same individual, do **not** treat them as independent biological replicates. Use donor-aware/mixed-effect modeling for independent cohort inference.
- Inferred cortical layer patterns are exploratory if annotated layers are present; a layer label must be checked before biological claims.

### GitHub readiness / originality

All scripts, a documented pipeline, versioned environment, falsifiable hypothesis, exploratory quality checks, reproducible random seeds, and report templates are included. Do not upload bulk human donor count matrices. Cite the data producers and official LIBD resource.

The project does NOT claim discovery of new psychiatric mechanisms: it is an independently written analysis of pre-existing public data.

### Preview included

`examples/synthetic_preview/` contains clearly marked simulated figures (spatial maps, estimated cell composition, layer patterns). The original package contains **no human donor data**. Run `bash run_demo.sh` to reproduce the synthetic results, then `bash run_real.sh` after R dependencies are installed to compute results from public measured data.

## Publish to GitHub

See `docs/GITHUB_PUBLICATION.md` for a safety-checked sequence.
