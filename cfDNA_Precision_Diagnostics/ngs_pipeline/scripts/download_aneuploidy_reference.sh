#!/usr/bin/env bash

set -euo pipefail

OUT="ngs_pipeline/data/reference/aneuploidy"

mkdir -p "$OUT"

while read -r chrom accession; do

    echo "Downloading chromosome ${chrom}"

    curl -fL --retry 5 \
      --connect-timeout 30 \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&seq_start=30000001&seq_stop=30100000&rettype=fasta&retmode=text" \
      -o "$OUT/chr${chrom}.fa"

    samtools faidx "$OUT/chr${chrom}.fa"

    sleep 1

done <<'EOF'
1 NC_000001.11
2 NC_000002.12
3 NC_000003.12
13 NC_000013.11
18 NC_000018.10
21 NC_000021.9
EOF

echo "Reference downloads completed"
