#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data results
if [[ ! -s data/hapmap1.zip ]]; then
  echo 'Downloading official historical PLINK tutorial data...'
  if ! curl -fL --retry 2 'https://zzz.bwh.harvard.edu/plink/hapmap1.zip' -o data/hapmap1.zip; then
    echo 'Historical PLINK host unavailable; trying Biostar Handbook archive...'
    if ! curl -fL --retry 2 'https://data.biostarhandbook.com/data/hapmap1.zip' -o data/hapmap1.zip; then
      echo 'Both archive hosts failed. Download hapmap1.zip manually and save it at data/hapmap1.zip' >&2
      exit 1
    fi
  fi
fi
unzip -qo data/hapmap1.zip -d data
[[ -s data/hapmap1.ped && -s data/hapmap1.map ]] || {
  echo 'Expected hapmap1.ped and hapmap1.map; inspect archive path:' >&2
  unzip -l data/hapmap1.zip | head -25 >&2
  exit 1
}
printf 'HapMap PED/MAP ready.\n'
