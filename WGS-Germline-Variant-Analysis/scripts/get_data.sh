#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data results
REGION="${REGION:-chr20:10000000-11000000}"
BAM_URL='https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam'
TRUTH_URL='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38'

if [[ ! -s data/chr20.fa ]]; then
  echo 'Downloading GRCh38 chr20 reference (UCSC hg38)...'
  curl -fL --retry 3 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz' | gzip -dc > data/chr20.fa
fi
samtools faidx data/chr20.fa
# GATK HaplotypeCaller requires a matching .dict file
if [[ ! -s data/chr20.dict ]]; then
  gatk CreateSequenceDictionary -R data/chr20.fa -O data/chr20.dict
fi

if [[ ! -s data/HG002.truth.vcf.gz ]]; then
  curl -fL --retry 3 "$TRUTH_URL/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" -o data/HG002.truth.vcf.gz
fi
if [[ ! -s data/HG002.truth.vcf.gz.tbi ]]; then
  curl -fL --retry 3 "$TRUTH_URL/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" -o data/HG002.truth.vcf.gz.tbi
fi
if [[ ! -s data/HG002.high_confidence.bed ]]; then
  curl -fL --retry 3 "$TRUTH_URL/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" -o data/HG002.high_confidence.bed
fi

if [[ ! -s data/HG002.region.bam ]]; then
  if [[ "${FULL_BAM:-0}" == '1' ]]; then
    echo 'Downloading FULL chr20 BAM (multiple GB) and index: this may take a while.'
    curl -fL --retry 3 -C - "$BAM_URL" -o data/HG002.chr20.bam
    curl -fL --retry 3 "$BAM_URL.bai" -o data/HG002.chr20.bam.bai
    samtools view -b -o data/HG002.region.bam data/HG002.chr20.bam "$REGION"
  else
    echo 'Trying a remotely indexed 1 Mb extraction (requires HTTPS range access)...'
    if ! samtools view -b -o data/HG002.region.bam "$BAM_URL" "$REGION"; then
      rm -f data/HG002.region.bam
      echo 'Remote extraction unavailable. Retry: FULL_BAM=1 bash scripts/get_data.sh' >&2
      exit 1
    fi
  fi
fi
samtools quickcheck -v data/HG002.region.bam
samtools index -@ 2 data/HG002.region.bam
printf 'Reference and variant truth data prepared; sample region BAM ready.\n'
