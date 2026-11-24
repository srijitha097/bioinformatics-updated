#!/usr/bin/env bash
# Replay the analyses from archived inputs, from the repository root.
set -euo pipefail
for tool in samtools bcftools; do
  version=$("$tool" --version | sed -n '1p')
  if [[ "$version" != "$tool 1.25" ]]; then
    echo "Expected $tool 1.25, found $version" >&2
    exit 1
  fi
done
d=exome_gene_panel_cnv_analysis/data
for pair in ARUP_SeqCap_EZ_Exome:wes GIAB_v3.3.2_NA12878_HC_regions:confident; do
 awk 'BEGIN{OFS="\t"} {sub(/^chr/,"",$1); if ($1 ~ /^([1-9]|1[0-9]|2[0-2])$/) print $1,$2,$3}' "$d/${pair%:*}.bed" > "$d/verify_${pair#*:}.bed"
done
p=rare_disease_variant_curation
python "$p/run.py" prioritize --clinvar "$p/data/variant_summary.txt.gz" --hpo "$p/data/genes_to_phenotype.txt" --cases "$p/cases.tsv" --limit 20 --out "$p/results"
python "$p/review_gjb2.py" --data "$p/data" --out "$p/results" --retrieved 2026-11-24
p=oncology_somatic_variant_interpretation
python "$p/run.py" extract --assertions "$p/data/01-Nov-2026-AcceptedAssertionSummaries.tsv" --accepted-only-export --out "$p/results"
python "$p/run.py" audit --source "$p/results/accepted_assertions.tsv" --review "$p/review_completed.tsv" --evidence "$p/data/01-Nov-2026-AcceptedClinicalEvidenceSummaries.tsv" --out "$p/results"
p=exome_gene_panel_cnv_analysis
o="$p/results"
python "$p/make_panel.py" --genes "$d/Omim_Genes.bed" --capture "$d/ARUP_SeqCap_EZ_Exome.bed" --out "$o/hereditary_cancer_panel.bed"
for label in wes panel; do
bed="$d/ARUP_SeqCap_EZ_Exome.bed"
if [ "$label" = panel ]; then bed="$o/hereditary_cancer_panel.bed"; fi
python "$p/run.py" --truth "$d/GIAB_v3.3.2_NA12878-decomposed-normalized.vcf.gz" --query "$d/converted_ARUP_NA12878_Exome-decomposed-normalized.vcf.gz" --targets "$bed" --confidence "$d/GIAB_v3.3.2_NA12878_HC_regions.bed" --autosomes-only --out "$o/$label.json"
done
bash "$p/crosscheck.sh" > "$o/benchmark_crosscheck.tsv"
python "$p/slice_depth.py" --data "$d" --out "$o"
python "$p/cnv_demo.py" --out "$o/synthetic"
python plot_results.py
python comparative_genomics/ranking.py --select
for app in get_af ranking citation_audit select_oncology reviews benchmark; do python comparative_genomics/$app.py; done
bash comparative_genomics/normalize.sh
for app in check_benchmark downsample array_compare cnv figures final_checks; do python comparative_genomics/$app.py; done
python -m unittest discover -s tests -v
python - <<'PY'
from hashlib import file_digest
from pathlib import Path
for line in Path("comparative_genomics/source_hashes.txt").read_text().splitlines():
    digest, path = line.split("  ", 1)
    with open(path, "rb") as stream:
        if file_digest(stream, "sha256").hexdigest() != digest:
            raise SystemExit(f"Source checksum mismatch: {path}")
PY
