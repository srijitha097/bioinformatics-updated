#!/usr/bin/env bash
set -euo pipefail

RAW="data/raw"
MATRIX_DIR="${RAW}/filtered_feature_bc_matrix"
SPATIAL_DIR="${RAW}/spatial"
TMP="${RAW}/tmp_spatial"

mkdir -p "${MATRIX_DIR}" "${SPATIAL_DIR}" "${TMP}"

BASE="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6506115&format=file&file="

echo "Downloading barcodes..."
curl -L "${BASE}GSM6506115_SP6_barcodes.tsv.gz" \
  -o "${MATRIX_DIR}/barcodes.tsv.gz"

echo "Downloading features..."
curl -L "${BASE}GSM6506115_SP6_features.tsv.gz" \
  -o "${MATRIX_DIR}/features.tsv.gz"

echo "Downloading matrix..."
curl -L "${BASE}GSM6506115_SP6_matrix.mtx.gz" \
  -o "${MATRIX_DIR}/matrix.mtx.gz"

echo "Downloading spatial files..."
curl -L "${BASE}GSM6506115_SP6_spatial.zip" \
  -o "${RAW}/spatial.zip"

echo "Extracting spatial files..."
rm -rf "${TMP}"
mkdir -p "${TMP}"
unzip -q -o "${RAW}/spatial.zip" -d "${TMP}"

# Copy the standard Visium spatial files regardless of nesting inside the ZIP.
for pattern in \
  "tissue_positions_list.csv" \
  "tissue_positions.csv" \
  "scalefactors_json.json" \
  "tissue_lowres_image.png" \
  "tissue_hires_image.png"
do
  f=$(find "${TMP}" -type f -name "${pattern}" | head -n 1 || true)
  if [[ -n "${f}" ]]; then
    cp "${f}" "${SPATIAL_DIR}/$(basename "${f}")"
  fi
done

echo
echo "Matrix files:"
ls -lh "${MATRIX_DIR}"

echo
echo "Spatial files:"
ls -lh "${SPATIAL_DIR}"

echo
echo "Download complete."
