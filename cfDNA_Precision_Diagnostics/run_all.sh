#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python -m pytest -q
python -m molecular_quantification.run
python -m prenatal_genomic_screening.run
python -m single_gene_screening.run
