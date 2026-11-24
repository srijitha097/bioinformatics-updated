#!/usr/bin/env bash
# Recompute exact genotype counts via htslib region filtering and sorted TSV sets.
set -euo pipefail
d=exome_gene_panel_cnv_analysis/data
o=exome_gene_panel_cnv_analysis/results
bc=${BCFTOOLS:-bcftools}
export LC_ALL=C
truth="$d/GIAB_v3.3.2_NA12878-decomposed-normalized.vcf.gz"
query="$d/converted_ARUP_NA12878_Exome-decomposed-normalized.vcf.gz"
for input in "$truth" "$query"; do
  if [[ ! -f "$input.tbi" && ! -f "$input.csi" ]]; then
    "$bc" index -t "$input"
  fi
done
tmp=$(mktemp -d)
trap 'rm -f "$tmp/confident.bed" "$tmp/targets.bed" "$tmp/truth.tsv" "$tmp/query.tsv"; rmdir "$tmp"' EXIT
autosomal_bed() {
  awk 'BEGIN{OFS="\t"} {sub(/^chr/,"",$1); if ($1 ~ /^([1-9]|1[0-9]|2[0-2])$/) print $1,$2,$3}' "$1"
}
autosomal_bed "$d/GIAB_v3.3.2_NA12878_HC_regions.bed" > "$tmp/confident.bed"
version=$("$bc" --version | sed -n '1p')
printf 'scope\ttool\ttruth\tquery\ttp\tfp\tfn\n'
for scope in wes panel; do
  bed="$d/ARUP_SeqCap_EZ_Exome.bed"
  if [[ "$scope" = panel ]]; then bed="$o/hereditary_cancer_panel.bed"; fi
  autosomal_bed "$bed" > "$tmp/targets.bed"
  for kind in truth query; do
    input="$truth"
    if [[ "$kind" = query ]]; then input="$query"; fi
    "$bc" view -R "$tmp/targets.bed" -T "$tmp/confident.bed" --regions-overlap 0 --targets-overlap 0 \
      -v snps -m2 -M2 -f PASS -Ou "$input" |
      "$bc" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' |
      awk 'BEGIN{OFS="\t"} $3 ~ /^[ACGT]$/ && $4 ~ /^[ACGT]$/ {
        if($5 ~ /^(0[|\/]1|1[|\/]0)$/) $5="0/1";
        else if($5 ~ /^1[|\/]1$/) $5="1/1";
        else next;
        print $1,$2,$3,$4,$5
      }' | sort -u > "$tmp/$kind.tsv"
  done
  t=$(wc -l < "$tmp/truth.tsv")
  q=$(wc -l < "$tmp/query.tsv")
  tp=$(comm -12 "$tmp/truth.tsv" "$tmp/query.tsv" | wc -l)
  printf '%s\t%s\t%d\t%d\t%d\t%d\t%d\n' "$scope" "$version" "$t" "$q" "$tp" "$((q - tp))" "$((t - tp))"
done
