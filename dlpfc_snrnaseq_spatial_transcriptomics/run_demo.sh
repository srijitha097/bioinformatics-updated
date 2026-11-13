#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ -f data/processed/real_data_provenance.txt ]]; then
  echo 'Real-data intermediates detected; demo would overwrite them. Move data/processed out of this directory before running synthetic demo.'; exit 2
fi
python scripts/00_make_demo.py
python scripts/02_deconvolve.py
python scripts/03_spatial_stats.py
python scripts/04_gene_holdout_validation.py
echo 'Synthetic demo finished. DO NOT describe the synthetic output as experimental results.'
