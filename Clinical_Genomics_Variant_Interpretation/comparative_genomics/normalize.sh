#!/usr/bin/env bash
set -euo pipefail
d=exome_gene_panel_cnv_analysis/data
for f in "$d"/GIAB*.vcf.gz "$d"/converted_ARUP*.vcf.gz; do
  name=$(basename "$f")
  bcftools view -R "$d/verify_wes.bed" -f PASS -m2 -M2 -Ou "$f" |
    bcftools norm -c e -f comparative_genomics/data/b37.fa -Oz -o "comparative_genomics/data/$name" 2> "comparative_genomics/sequencing/$name.norm.log"
done
