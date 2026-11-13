#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/raw
BASE='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254569'
echo 'GSE254569 GEO: approximately 10.7 GB compressed plus decompressed data.'
if [[ ! -s data/raw/GSE254569_adata_RNA.h5ad.gz ]]; then
  curl --fail --location --retry 5 --retry-delay 5 --continue-at - \
    --output data/raw/GSE254569_adata_RNA.h5ad.gz "$BASE/suppl/GSE254569_adata_RNA.h5ad.gz"
else
  echo 'Compressed H5AD already exists; verify its integrity if a prior download was interrupted.'
fi
if [[ ! -s data/raw/GSE254569_family.soft.gz ]]; then
  curl --fail --location --retry 5 --retry-delay 5 \
    --output data/raw/GSE254569_family.soft.gz "$BASE/soft/GSE254569_family.soft.gz"
fi
if [[ ! -s data/raw/GSE254569_adata_RNA.h5ad ]]; then
  echo 'Checking gzip archive before extraction ...'
  gzip -t data/raw/GSE254569_adata_RNA.h5ad.gz
  gzip -dk data/raw/GSE254569_adata_RNA.h5ad.gz
fi
python scripts/00_parse_geo_soft.py
echo 'Download done. NEXT: python scripts/00_inspect_h5ad.py'
