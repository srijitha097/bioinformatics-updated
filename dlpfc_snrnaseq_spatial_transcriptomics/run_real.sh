#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
command -v Rscript >/dev/null 2>&1 || { echo 'Rscript missing. Install R and R/Bioconductor packages from README.md' >&2; exit 2; }
Rscript -e 'stopifnot(all(vapply(c("spatialLIBD","HDF5Array","SpatialExperiment","SummarizedExperiment","Matrix"), requireNamespace, logical(1), quietly=TRUE)))'
if [[ -d results && -n "$(ls -A results 2>/dev/null)" ]]; then
  stamp=$(date +%Y%m%d_%H%M%S)
  mv results "results_backup_${stamp}"
  echo "Previous results saved to results_backup_${stamp}"
fi
mkdir -p results
Rscript scripts/01_fetch_libd.R
python scripts/02_deconvolve.py
python scripts/03_spatial_stats.py
python scripts/04_gene_holdout_validation.py
python scripts/05_make_report.py --real
echo 'Real LIBD analysis ran. Review section-level estimates: no disease-case/control comparison is possible.'
