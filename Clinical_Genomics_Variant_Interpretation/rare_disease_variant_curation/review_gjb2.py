"""Reconstruct two published ClinGen assessments and tabulate live gnomAD data."""

import argparse
import json
from pathlib import Path

from run import audit, rows, sha, write

SOURCE = "https://search.clinicalgenome.org/kb/conditions/MONDO%3A0019497/by-gene"
REVIEWS = {
    "17000": (
        "PS4 (Strong), PP1 (Supporting), PP3 (Supporting), PM3 (Very Strong)",
        "Case enrichment in biallelic genotypes (European OR 29), observed in trans with "
        "truncating alleles and segregation supported the expert-panel assessment. PP3 "
        "came from the panel's REVEL assessment. Reduced penetrance and mild hearing "
        "loss limit interpretation of unphenotyped population controls.",
        "https://erepo.clinicalgenome.org/evrepo/ui/classification/e8b1ce82-bc97-4baf-b442-294c2a6849cd",
        "nfe",
    ),
    "17023": (
        "PS4 (Strong), PP1 (Strong), PM3 (Moderate)",
        "Case enrichment in biallelic genotypes (Asian OR 19) and familial segregation "
        "supported the expert-panel assessment. The panel retained moderate PM3 despite "
        "many compound heterozygotes given the high population frequency. Reduced "
        "penetrance makes a high ancestry-specific allele frequency compatible with this assessment.",
        "https://actionability.clinicalgenome.org/evrepo/ui/classification/69193dd4-e172-4874-984a-28596f644cae?version=1.0",
        "eas",
    ),
}


def review(data, out, retrieved):
    queue = out / "candidate_queue.tsv"
    candidates = {row["variation_id"]: row for row in rows(queue)}
    assessments, frequencies = [], []
    for vid, (codes, rationale, expert, ancestry) in REVIEWS.items():
        candidate = candidates[vid]
        path = data / f"gnomad_{vid}.json"
        variant = json.loads(path.read_text())["data"]["variant"]
        frequency_strings = []
        for dataset in ("exome", "genome"):
            counts = variant[dataset]
            population = next(p for p in counts["populations"] if p["id"] == ancestry)
            for label, record in (("all", counts), (ancestry, population)):
                ac, an = record["ac"], record["an"]
                af = ac / an
                frequencies.append(
                    dict(
                        variation_id=vid,
                        variant_id=candidate["variant_id"],
                        dataset=dataset,
                        group=label,
                        ac=ac,
                        an=an,
                        af=af,
                        api_dataset="gnomad_r4",
                        retrieved=retrieved,
                        raw_sha256=sha(path),
                        url=candidate["gnomad_url"],
                    )
                )
                frequency_strings.append(f"{dataset} {label} {ac}/{an} ({af:.3%})")
        assessments.append(
            dict(
                variation_id=vid,
                gene="GJB2",
                condition=candidate["condition"],
                reviewed_condition="DFNB1A (OMIM 220290)",
                inheritance="Autosomal recessive",
                gene_disease_validity_source=SOURCE,
                gnomad_version=f"gnomad_r4 API, retrieved {retrieved}",
                population_af=" | ".join(frequency_strings),
                pmids="31160754",
                acmg_codes_and_strengths=codes,
                research_classification="Pathogenic",
                evidence_rationale=f"{candidate['hgvs']} in DFNB1A (OMIM 220290, HPO HP:0000407). "
                f"Published ClinGen Hearing Loss VCEP v1 ACMG/AMP reconstruction. {rationale} "
                f"The later v2 specifications exclude this allele from the general frequency rules. "
                f"Evidence overlap across codes was not independently re-adjudicated. Expert record {expert} .",
                reviewer="Codex literature reconstruction",
                review_date=retrieved,
            )
        )
    write(out / "gnomad_frequencies.tsv", frequencies, list(frequencies[0]))
    curated = out / "curation_completed.tsv"
    write(curated, assessments, list(assessments[0]))
    audit(queue, curated, out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--retrieved", required=True)
    args = parser.parse_args()
    review(args.data, args.out, args.retrieved)
