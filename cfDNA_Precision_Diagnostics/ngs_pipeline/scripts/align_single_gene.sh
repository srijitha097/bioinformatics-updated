#!/usr/bin/env bash

set -euo pipefail

REF="ngs_pipeline/data/reference/single_gene/combined_reference.fa"

DATA="ngs_pipeline/data/single_gene"

OUT="results/ngs/single_gene"

THREADS="${THREADS:-4}"

mkdir -p "$OUT"/{alignment,logs,fastqc,multiqc,benchmark}

if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "$REF"
fi

if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

TOTAL=0

for R1 in "$DATA"/*_R1.fastq.gz; do

    SAMPLE=$(basename "$R1" _R1.fastq.gz)

    R2="$DATA/${SAMPLE}_R2.fastq.gz"

    BAM="$OUT/alignment/${SAMPLE}.bam"

    TOTAL=$((TOTAL + 1))

    if [[ ! -s "$R2" ]]; then
        echo "Missing R2: $SAMPLE"
        exit 1
    fi

    if [[ -s "$BAM" && -s "${BAM}.bai" ]] \
        && samtools quickcheck "$BAM"; then

        echo "[$TOTAL/135] Existing BAM: $SAMPLE"

    else

        echo "[$TOTAL/135] Aligning: $SAMPLE"

        bwa mem \
            -t "$THREADS" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
            "$REF" "$R1" "$R2" \
            2> "$OUT/logs/${SAMPLE}.bwa.log" \
        | samtools sort \
            -@ 2 \
            -o "${BAM}.tmp" -

        samtools quickcheck "${BAM}.tmp"

        mv "${BAM}.tmp" "$BAM"

        samtools index -@ 2 "$BAM"

    fi

    if [[ ! -s "$OUT/alignment/${SAMPLE}.flagstat.txt" ]]; then

        samtools flagstat "$BAM" \
            > "$OUT/alignment/${SAMPLE}.flagstat.txt"

    fi

    if [[ ! -s "$OUT/alignment/${SAMPLE}.stats.txt" ]]; then

        samtools stats "$BAM" \
            > "$OUT/alignment/${SAMPLE}.stats.txt"

    fi

    if [[ ! -s "$OUT/alignment/${SAMPLE}.coverage.txt" ]]; then

        samtools coverage "$BAM" \
            > "$OUT/alignment/${SAMPLE}.coverage.txt"

    fi

    echo "Completed: $SAMPLE"

done

echo ""
echo "ALIGNMENT COMPLETE"
echo "Samples processed: $TOTAL"

samtools quickcheck -v "$OUT"/alignment/*.bam

echo "All BAM files passed integrity checks"
