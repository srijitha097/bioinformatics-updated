See [Exploratory Results](docs/FINAL_RESULTS.md) for the completed analysis and limitations.

> **Start here on macOS:** read `../START_HERE_MAC.md` if using the combined ZIP. Keep your original demo folder backed up before changing files.

# Original Research Project 01 — Schizophrenia-associated Cell-Type Programs in Human Brain snRNA-seq

**Portfolio status:** Real-data exploratory analysis completed locally. Full-cohort validation and astrocyte contamination assessment remain outstanding.

## Research question and independently designed analysis

Which neuronal and glial populations in postmortem human orbitofrontal cortex show schizophrenia-associated gene expression after adjustment for donor-level technical/biological covariates? Can those changes be interpreted through molecular gene programs?

This is an **original secondary-analysis design**, not a replication of the dataset authors' published claims. The project uses public existing measurements rather than generating new patient tissue. It does **not** establish causality, clinical utility, or novel target validation.

### Dataset

- NCBI GEO accession **GSE254569**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254569
- Human orbitofrontal cortex, postmortem human single-*nucleus* RNA sequencing.
- The GEO record advertises one **10.7 GB compressed** processed AnnData `.h5ad.gz`, containing annotations; the download expands substantially. Plan significant disk space and possibly HPC resources.
- Sample metadata (diagnosis, age, sex, PMI, RNA integrity where present) originates from GEO, not guessed. Exclude any non-schizophrenia psychiatric diagnoses from this specific contrast.
- Original investigators retain credit for **data generation**. Cite the GEO accession and original data-generation study in a GitHub repo/publication; citing provenance does not mean reproducing their analysis.

## Quick start — **no downloads, no R required for the Python stages**

Tested with Python 3.13, NumPy, pandas, SciPy, matplotlib.

```bash
cd 01_schizophrenia_snrnaseq
conda env create -f environment.yml
conda activate libd-scz
bash run_demo.sh
python -m pytest -q tests
```

The demo generates *synthetic* counts for 24 simulated donors, four simulated cell types, and built-in case/control signal. It tests data handling only. **Never report demo p-values, simulated cell counts, or demo plots as scientific findings.** `run_demo.sh` will run DESeq2 if R is installed, otherwise it stops after donor-level pseudobulk preparation without pretending differential expression was completed.

## Real run — download and run from public data

```bash
# 1) Data download: ~10.7 GB compressed plus uncompressed copy
bash scripts/00_download_real.sh

# 2) Inspect exact columns and whether raw integer count matrix is present
python scripts/00_inspect_h5ad.py

# 3) Run full pipeline (initial, memory-limited per-donor subset)
bash run_real.sh

# Optional explicit metadata mapping when auto-detection does not match the file
python scripts/01_import_real.py --donor-col SAMPLE_COLUMN \
  --celltype-col CELL_TYPE_COLUMN --diagnosis-col DIAGNOSIS_COLUMN \
  --max-nuclei-per-donor 150
python scripts/01_prepare_pseudobulk.py
Rscript scripts/02_deseq2.R
```

If AnnData contains **log-normalized expression but no raw integer counts**, stop. DESeq2 cannot use normalized/log values as raw counts; retrieve original UMI matrices or an appropriate raw-count layer. The importer checks for non-integer values and errors. GEO SOFT is parsed automatically into `data/raw/geo_donor_metadata.csv` to infer diagnosis if it is absent from AnnData `obs` (the donor label must still be mappable). If the schema does not match, inspect columns and provide explicit flags rather than modifying diagnoses by guesswork.

`run_real.sh` accepts additional arguments forwarded to the importer, e.g. `bash run_real.sh --donor-col Donor --celltype-col celltype --diagnosis-col classification`. Initial results use a **subsample of nuclei per donor**, intentionally for manageable local validation. Before claiming comprehensive dataset results, run with appropriately larger per-donor settings and compute resource estimates; do not call a subset a full cohort analysis. The data were originally quality-controlled by the submitters; this project performs post-processing checks, not a raw FASTQ reprocessing workflow.

### R setup if Conda DESeq2 is unavailable

```r
install.packages('BiocManager')
BiocManager::install('DESeq2', ask=FALSE, update=FALSE)
```

## What happens in each step?

| Script | Input | Action | Output |
|---|---|---|---|
| `scripts/00_make_demo.py` | random generator | Synthetic count/metadata test data | `data/processed` |
| `scripts/00_download_real.sh` | GEO | Download AnnData & SOFT metadata | `data/raw` |
| `scripts/00_parse_geo_soft.py` | SOFT | Parse sample metadata | `data/raw/geo_donor_metadata.csv` |
| `scripts/00_inspect_h5ad.py` | `.h5ad` | Expose available `obs`, `var`, count layers | console |
| `scripts/01_import_real.py` | AnnData + GEO metadata | Map conditions & cell types; select raw integer UMIs; donor stratification | `counts.npz`, `cells.csv`, `genes.csv` |
| `scripts/01_prepare_pseudobulk.py` | raw sparse count matrix | Aggregate by **donor × cell type** | `results/pseudobulk/*`, composition plot |
| `scripts/02_deseq2.R` | pseudobulk UMI counts | Negative binomial case/control model; FDR BH | DE tables, volcano plots |
| `scripts/03_enrichment.py` | DE tables + a valid GMT | GO / pathway hypergeometric enrichment in *tested gene universe* | `results/enrichment/*` |

Optional enrichment (after real DE and downloading an appropriate GMT from a recognized gene-set source):

```bash
python scripts/03_enrichment.py --gmt path/to/your_gene_sets.gmt
```

## Statistical design (important for Staff Scientist interviews)

- The independent biological replicate is **the donor**, NOT the nucleus.
- Pseudobulk counts are summed per donor and cell type; a negative-binomial generalized model (`DESeq2`) compares schizophrenia vs control.
- Age, sex, and postmortem interval are added only if present, non-missing, and estimable. Inspect the design summary. Consider RIN/batch and confounding with disease in an expanded preregistered analysis.
- Filter low-count genes; control false discovery rate (Benjamini–Hochberg); interpret effect sizes alongside `padj`.
- Diagnostic groups, number of donors, sex balance, library batch, RNA quality, and postmortem artifacts can confound association. The package does not magically remove these limitations.
- Cell-type abundance plots show **descriptive** sampled nuclei proportions, not validated disease-associated composition tests.

## Results expected **only after actual data analysis**

```text
results/
  donor_celltype_counts.csv
  pseudobulk_manifest.csv
  pseudobulk/<cell-type>/counts.csv, metadata.csv
  figures/01_celltype_composition.png
  deseq2/<cell-type>_deseq2.csv
  deseq2/<cell-type>_volcano.png
  deseq2/<cell-type>_summary.txt
  enrichment/<cell-type>_deseq2_enrichment.csv  [optional]
```

## Responsible publication and data licensing

Do not commit the original H5AD, raw files, pseudobulk matrices containing donor-linked measurements, or identifiable sample metadata. Share scripts, safe aggregated visualizations, and a transparent provenance statement instead. Record GEO accession, original publication, software versions, random seeds, and exact sample selections.

## Resume/portfolio

See `docs/RESUME_AFTER_RUN.md`; **do not claim the analysis was completed until the real-data pipeline has been run and results reviewed**.

### Preview included

`examples/synthetic_preview/` includes one clearly marked simulated composition figure. The raw synthetic input and intermediate result files are generated when you run `bash run_demo.sh`; they are deliberately absent from the initial download and excluded from GitHub by `.gitignore`.

## Publish to GitHub

See `docs/GITHUB_PUBLICATION.md` for a safety-checked sequence.
