#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw/GSE147352

URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147352/suppl/GSE147352_RAW.tar"
TAR="data/raw/GSE147352/GSE147352_RAW.tar"

echo "Downloading GSE147352 raw HTSeq count files..."
wget -c "${URL}" -O "${TAR}"

echo "Extracting..."
tar -xf "${TAR}" -C data/raw/GSE147352

echo "Count files:"
find data/raw/GSE147352 -type f -name "GSM*_htseq_count*.txt.gz" | wc -l

echo "Done."
