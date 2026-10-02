#!/usr/bin/env bash
set -euo pipefail

# Example only.
# Supply matching GRCh38 HISAT2 index and GTF paths.

R1="${1:?R1 FASTQ required}"
R2="${2:?R2 FASTQ required}"
SAMPLE="${3:?sample ID required}"
INDEX_PREFIX="${4:?HISAT2 index prefix required}"
GTF="${5:?GTF required}"

mkdir -p results/bam results/counts

hisat2 \
  -x "${INDEX_PREFIX}" \
  -1 "${R1}" \
  -2 "${R2}" \
  -p 8 \
| samtools sort -@ 8 -o "results/bam/${SAMPLE}.sorted.bam"

samtools index "results/bam/${SAMPLE}.sorted.bam"

htseq-count \
  -f bam \
  -r pos \
  -s reverse \
  -t exon \
  -i gene_id \
  "results/bam/${SAMPLE}.sorted.bam" \
  "${GTF}" \
  > "results/counts/${SAMPLE}.counts.txt"
