#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
[[ -f data/raw/GSE254569_adata_RNA.h5ad ]] || { echo 'Download real H5AD first: bash scripts/00_download_real.sh' >&2; exit 2; }
[[ -f data/raw/geo_donor_metadata.csv ]] || python scripts/00_parse_geo_soft.py
python scripts/00_inspect_h5ad.py
# Do not delete output without retaining the user's already completed synthetic demo.
if [[ -d results && -n "$(ls -A results 2>/dev/null)" ]]; then
  stamp=$(date +%Y%m%d_%H%M%S)
  mv results "results_backup_${stamp}"
  echo "Previous results saved to results_backup_${stamp}"
fi
mkdir -p results
python scripts/01_import_real.py "$@"
python scripts/01_prepare_pseudobulk.py
Rscript scripts/02_deseq2.R
python scripts/04_make_report.py --real
printf '\nREAL exploratory analysis complete; verify metadata mapping, covariates, and QC before interpreting.\n'
