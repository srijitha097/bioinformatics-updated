#!/usr/bin/env bash
# Optional: demonstrate read extraction, FastQC, bwa-mem2 and GATK MarkDuplicates.
# The input BAM is already aligned/deduplicated. Extracted reads are NOT original FASTQs.
set -euo pipefail
cd "$(dirname "$0")/.."
for tool in samtools bwa-mem2 gatk; do command -v "$tool" >/dev/null || { echo "Missing $tool" >&2; exit 1; }; done
[[ -s data/HG002.region.bam && -s data/chr20.fa ]] || { echo 'First run scripts/get_data.sh' >&2; exit 1; }
mkdir -p data results
# Discard incomplete read pairs (some mates may lie outside the extracted region).
samtools collate -u -O data/HG002.region.bam | \
  samtools fastq -@ 2 -n -1 data/HG002_R1.fastq.gz -2 data/HG002_R2.fastq.gz \
    -0 /dev/null -s /dev/null -
if command -v fastqc >/dev/null; then
  fastqc -t 2 -o results data/HG002_R1.fastq.gz data/HG002_R2.fastq.gz
fi
bwa-mem2 index data/chr20.fa
bwa-mem2 mem -t 4 -R '@RG\tID:HG002_RG1\tSM:HG002\tPL:ILLUMINA' \
  data/chr20.fa data/HG002_R1.fastq.gz data/HG002_R2.fastq.gz | \
  samtools sort -@ 2 -o data/HG002.realigned.bam -
gatk --java-options '-Xmx4g' MarkDuplicates \
  -I data/HG002.realigned.bam -O data/HG002.realigned.dedup.bam \
  -M results/markduplicates_metrics.txt
samtools index data/HG002.realigned.dedup.bam
if command -v multiqc >/dev/null; then multiqc results -o results/multiqc || true; fi
echo 'Optional realignment complete. For the main analysis, run:'
echo '  BAM=data/HG002.realigned.dedup.bam bash scripts/run.sh'
echo 'Caution: region-extracted read pairs are not a full original WGS FASTQ experiment.'
