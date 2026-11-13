#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
REGION="${REGION:-chr20:10000000-11000000}"
BAM="${BAM:-data/HG002.region.bam}"
REF="${REF:-data/chr20.fa}"
TRUTH="${TRUTH:-data/HG002.truth.vcf.gz}"
BED="${BED:-data/HG002.high_confidence.bed}"
OUT="${OUT:-results}"
mkdir -p "$OUT"
for path in "$BAM" "$REF" "$TRUTH" "$BED"; do [[ -s "$path" ]] || { echo "Missing: $path" >&2; exit 1; }; done
bash scripts/check_tools.sh

# Region BAM comes from a preprocessed/duplicate-marked alignment dataset.
samtools flagstat "$BAM" > "$OUT/flagstat.txt"
samtools depth -aa -r "$REGION" "$BAM" > "$OUT/depth.tsv"
python3 scripts/coverage.py "$OUT/depth.tsv" "$OUT/coverage_summary.tsv"

# Germline call; this is regional evaluation of a WGS source, not whole genome.
gatk --java-options '-Xmx4g' HaplotypeCaller \
  -R "$REF" -I "$BAM" -L "$REGION" -O "$OUT/variants.raw.vcf.gz"

# Subset then normalize GIAB to the same reference/region.
bcftools view -r "$REGION" -Oz -o "$OUT/truth.region.vcf.gz" "$TRUTH"
bcftools index -f -t "$OUT/truth.region.vcf.gz"
bcftools norm -f "$REF" -m -any -Oz -o "$OUT/truth.norm.vcf.gz" "$OUT/truth.region.vcf.gz"
bcftools index -f -t "$OUT/truth.norm.vcf.gz"

bcftools norm -f "$REF" -m -any -Oz -o "$OUT/variants.norm.vcf.gz" "$OUT/variants.raw.vcf.gz"
bcftools index -f -t "$OUT/variants.norm.vcf.gz"
# Simple exploratory quality flag; NOT a clinical or GATK best-practice filter.
bcftools filter -s LOW_QUAL -e 'QUAL<30' -Oz \
  -o "$OUT/variants.filtered.vcf.gz" "$OUT/variants.norm.vcf.gz"
bcftools index -f -t "$OUT/variants.filtered.vcf.gz"
bcftools stats "$OUT/variants.filtered.vcf.gz" > "$OUT/variant_stats.txt"

python3 scripts/compare_vcfs.py \
  --calls "$OUT/variants.filtered.vcf.gz" --truth "$OUT/truth.norm.vcf.gz" \
  --bed "$BED" --region "$REGION" --output "$OUT/benchmark.json"
python3 scripts/report.py --benchmark "$OUT/benchmark.json" \
  --coverage "$OUT/coverage_summary.tsv" --out "$OUT/report.html"
{
  samtools --version | head -1
  bcftools --version | head -1
  gatk --version 2>&1 | head -2
  python3 --version
} > "$OUT/versions.txt" || true
printf 'DONE. Read %s/report.html and %s/benchmark.json\n' "$OUT" "$OUT"
