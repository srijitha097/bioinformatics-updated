#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ -f data/processed/real_data_provenance.txt ]]; then
 echo 'STOP: real data were imported. Running demo would overwrite them. Use a fresh copy of this repository.' >&2; exit 2
fi
python scripts/00_make_demo.py
python scripts/01_prepare_pseudobulk.py
if command -v Rscript >/dev/null 2>&1 && Rscript -e 'quit(status=ifelse(requireNamespace("DESeq2",quietly=TRUE),0,1))'; then
  Rscript scripts/02_deseq2.R
else
  echo 'R/DESeq2 unavailable. The demo stopped after pseudobulk preparation; DE has NOT run.'
fi
printf '\nSYNTHETIC DEMO ONLY — results are NOT biological findings.\n'
