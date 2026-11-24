#!/usr/bin/env bash

set -euo pipefail

REF="ngs_pipeline/data/reference/aneuploidy/combined_reference.fa"

DATA="ngs_pipeline/data/aneuploidy"

OUT="results/ngs/aneuploidy"

mkdir -p "$OUT"/{alignment,logs,fastqc,multiqc,benchmark}

TOTAL=0
PROCESSED=0
SKIPPED=0

for R1 in "$DATA"/*_R1.fastq.gz; do

    SAMPLE=$(basename "$R1" _R1.fastq.gz)

    R2="$DATA/${SAMPLE}_R2.fastq.gz"

    BAM="$OUT/alignment/${SAMPLE}.bam"

    TOTAL=$((TOTAL + 1))

    if [[ ! -s "$R2" ]]; then
        echo "Missing paired FASTQ: $R2"
        exit 1
    fi

    if [[ -s "$BAM" && -s "${BAM}.bai" \
        && -s "$OUT/alignment/${SAMPLE}.flagstat.txt" \
        && -s "$OUT/alignment/${SAMPLE}.stats.txt" \
        && -s "$OUT/alignment/${SAMPLE}.coverage.tsv" ]] \
        && samtools quickcheck "$BAM"; then

        echo "Skipping completed sample: $SAMPLE"

        SKIPPED=$((SKIPPED + 1))

        continue

    fi

    echo "========================================"
    echo "Processing sample $TOTAL: $SAMPLE"
    echo "========================================"

    bwa mem -t 4 \
        "$REF" \
        "$R1" \
        "$R2" \
        2> "$OUT/logs/${SAMPLE}.bwa.log" \
        | samtools sort -@ 2 \
          -o "${BAM}.tmp" -

    samtools quickcheck -v "${BAM}.tmp"

    mv "${BAM}.tmp" "$BAM"

    samtools index "$BAM"

    samtools flagstat "$BAM" \
        > "$OUT/alignment/${SAMPLE}.flagstat.txt"

    samtools stats "$BAM" \
        > "$OUT/alignment/${SAMPLE}.stats.txt"

    samtools coverage "$BAM" \
        > "$OUT/alignment/${SAMPLE}.coverage.tsv"

    PROCESSED=$((PROCESSED + 1))

    echo "Completed: $SAMPLE"

done

echo ""
echo "Batch alignment completed"
echo "Total samples: $TOTAL"
echo "Newly processed: $PROCESSED"
echo "Previously completed: $SKIPPED"

