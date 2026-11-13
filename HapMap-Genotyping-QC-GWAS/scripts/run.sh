#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
command -v plink >/dev/null || { echo 'Please install PLINK 1.9' >&2; exit 1; }
python3 -c 'import numpy,pandas,matplotlib' || { echo 'pip install -r requirements.txt' >&2; exit 1; }
[[ -s data/hapmap1.bed && -s data/hapmap1.bim && -s data/hapmap1.fam ]] || { echo 'Run scripts/get_data.sh first' >&2; exit 1; }
mkdir -p results

plink --bfile data/hapmap1 --make-bed --out results/raw
plink --bfile results/raw --missing --out results/initial
plink --bfile results/raw --freq --out results/initial_freq
plink --bfile results/raw --hardy --out results/initial_hwe

# Report sample/SNP QC from RAW files, then construct a FILTERED analysis set.
plink --bfile results/raw --mind 0.05 --geno 0.05 --maf 0.01 --hwe 1e-6 \
  --make-bed --out results/filtered

# Prune LD before calculating PCA and sample heterozygosity.
plink --bfile results/filtered --indep-pairwise 50 5 0.2 --out results/pruned
plink --bfile results/filtered --extract results/pruned.prune.in --pca 3 --out results/pca
plink --bfile results/filtered --extract results/pruned.prune.in --het --out results/heterozygosity

# Educational phenotype is simulated in historical PED/FAM.
plink --bfile results/filtered --assoc --allow-no-sex --out results/gwas_unadjusted
plink --bfile results/filtered --logistic --covar results/pca.eigenvec \
  --covar-number 1-3 --hide-covar --allow-no-sex --out results/gwas_adjusted

python3 scripts/report.py --results results
printf 'Complete: results/report.html and results/{manhattan,qq,pca}.png\n'
