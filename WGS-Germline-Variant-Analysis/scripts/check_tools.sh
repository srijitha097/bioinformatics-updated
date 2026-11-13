#!/usr/bin/env bash
set -euo pipefail
for tool in samtools bcftools gatk python3 curl; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: missing dependency $tool" >&2
    exit 1
  fi
done
printf 'All essential WGS tools found.\n'
