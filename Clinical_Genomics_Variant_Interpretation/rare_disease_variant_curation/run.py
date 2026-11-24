"""Build a candidate curation queue. ClinVar assertions are reference labels, not independent classifications."""

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import report, rows, sha, write

GENES = {
    "CFTR",
    "HBB",
    "GJB2",
    "PAH",
    "DMD",
    "BRCA1",
    "BRCA2",
    "MLH1",
    "MSH2",
    "MSH6",
    "PMS2",
    "LDLR",
    "MYBPC3",
    "MYH7",
    "PTEN",
    "TSC1",
    "TSC2",
    "MECP2",
    "SCN1A",
    "FBN1",
}


def prioritize(clinvar, hpo, cases, out, limit):
    if limit < 1:
        raise ValueError("Limit must be positive")
    associations = defaultdict(lambda: defaultdict(set))
    excluded_annotations = 0
    for r in rows(hpo):
        frequency = r.get("frequency", "")
        if frequency == "HP:0040285" or re.fullmatch(
            r"0(?:\.0+)?%|0/[1-9][0-9]*", frequency
        ):
            excluded_annotations += 1
            continue
        associations[r["gene_symbol"]][r.get("disease_id", "")].add(r["hpo_id"])
    case_rows = list(rows(cases))
    candidates = defaultdict(list)
    counts = Counter()
    seen = set()
    for r in rows(clinvar):
        if (
            r.get("Assembly") != "GRCh38"
            or r.get("Type") != "single nucleotide variant"
        ):
            continue
        gene = r.get("GeneSymbol", "")
        vid = r.get("VariationID", "")
        if gene not in GENES or not vid or vid in seen:
            continue
        review = r.get("ReviewStatus", "")
        sig = r.get("ClinicalSignificance", "")
        if (
            "reviewed by expert panel" not in review
            and "practice guideline" not in review
        ):
            continue
        if sig not in {
            "Pathogenic",
            "Likely pathogenic",
            "Benign",
            "Likely benign",
            "Uncertain significance",
        }:
            continue
        ref = r.get("ReferenceAlleleVCF", "")
        alt = r.get("AlternateAlleleVCF", "")
        pos = r.get("PositionVCF", "")
        chrom = r.get("Chromosome", "")
        if (
            ref not in "ACGT"
            or alt not in "ACGT"
            or len(ref) != 1
            or len(alt) != 1
            or not pos.isdigit()
            or not chrom
        ):
            continue
        seen.add(vid)
        counts[sig] += 1
        candidates[gene].append(
            dict(
                variation_id=vid,
                gene=gene,
                hgvs=r.get("Name", ""),
                condition=r.get("PhenotypeList", ""),
                condition_ids=r.get("PhenotypeIDS", ""),
                clinvar_reference=sig,
                review_status=review,
                variant_id=f"{chrom}-{pos}-{ref}-{alt}",
                clinvar_url=f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/",
                gnomad_url=f"https://gnomad.broadinstitute.org/variant/{chrom}-{pos}-{ref}-{alt}?dataset=gnomad_r4",
            )
        )
    ranked = []
    for c in case_rows:
        terms = set(filter(None, c["hpo_ids"].split(";")))
        if not terms or any(not re.fullmatch(r"HP:\d{7}", x) for x in terms):
            raise ValueError("Invalid HPO IDs")
        for gene, vs in candidates.items():
            disease_matches = {
                disease: terms & annotations
                for disease, annotations in associations[gene].items()
            }
            best_count = max(map(len, disease_matches.values()), default=0)
            if not best_count:
                continue
            best_diseases = sorted(
                disease
                for disease, matches in disease_matches.items()
                if len(matches) == best_count
            )
            matches = disease_matches[best_diseases[0]]
            # The score describes the gene-disease phenotype fit and is shared by the gene's variants.
            for v in sorted(vs, key=lambda x: int(x["variation_id"]))[:limit]:
                ranked.append(
                    dict(
                        v,
                        case_id=c["case_id"],
                        hpo_matches=";".join(sorted(matches)),
                        hpo_overlap=best_count / len(terms),
                        hpo_best_disease_ids=";".join(best_diseases),
                    )
                )
    ranked.sort(
        key=lambda r: (
            r["case_id"],
            -r["hpo_overlap"],
            r["gene"],
            int(r["variation_id"]),
        )
    )
    if not ranked:
        raise ValueError(
            "No candidates matched: verify HPO annotations, case terms and data schemas"
        )
    fields = [
        "case_id",
        "gene",
        "variation_id",
        "variant_id",
        "hgvs",
        "condition",
        "condition_ids",
        "hpo_matches",
        "hpo_overlap",
        "hpo_best_disease_ids",
        "clinvar_reference",
        "review_status",
        "clinvar_url",
        "gnomad_url",
    ]
    write(out / "candidate_queue.tsv", ranked, fields)
    write(
        out / "curation_template.tsv",
        [
            dict(
                variation_id=r["variation_id"], gene=r["gene"], condition=r["condition"]
            )
            for r in {x["variation_id"]: x for x in ranked}.values()
        ],
        [
            "variation_id",
            "gene",
            "condition",
            "inheritance",
            "gene_disease_validity_source",
            "gnomad_version",
            "population_af",
            "pmids",
            "acmg_codes_and_strengths",
            "evidence_rationale",
            "research_classification",
            "reviewer",
            "review_date",
        ],
    )
    report(
        out / "summary.json",
        dict(
            research_only=True,
            case_type="synthetic phenotype scenarios, not patient genomes",
            hpo_zero_frequency_rows_excluded=excluded_annotations,
            eligible_expert_reviewed_snvs=sum(counts.values()),
            reference_class_counts=dict(counts),
            candidate_case_pairs=len(ranked),
            unique_candidate_variants=len({r["variation_id"] for r in ranked}),
            input_sha256={str(p): sha(p) for p in (clinvar, hpo, cases)},
            limitation="Best exact HPO overlap within a single annotated disease per gene; not variant-level phenotype association or independent ACMG classification.",
        ),
    )


def audit(queue, curation, out):
    reference = {r["variation_id"]: r for r in rows(queue)}
    completed = []
    seen = set()
    required = [
        "inheritance",
        "gene_disease_validity_source",
        "gnomad_version",
        "population_af",
        "pmids",
        "acmg_codes_and_strengths",
        "evidence_rationale",
        "research_classification",
        "reviewer",
        "review_date",
    ]
    classes = {
        "Pathogenic",
        "Likely pathogenic",
        "Benign",
        "Likely benign",
        "Uncertain significance",
    }
    for r in rows(curation):
        if not r.get("research_classification"):
            continue
        vid = r["variation_id"]
        if vid not in reference or vid in seen:
            raise ValueError(f"Unknown or duplicate variant {vid}")
        if r["research_classification"] not in classes or any(
            not r.get(k) for k in required
        ):
            raise ValueError(f"Incomplete curation {vid}")
        if (
            r["gene"] != reference[vid]["gene"]
            or r["condition"] != reference[vid]["condition"]
        ):
            raise ValueError(f"Gene/condition mismatch {vid}")
        seen.add(vid)
        ref = reference[vid]["clinvar_reference"]
        completed.append(
            dict(
                r,
                reference_class=ref,
                concordant=int(r["research_classification"] == ref),
            )
        )
    if not completed:
        raise ValueError(
            "Complete the manual curation template first; no automated clinical classification is performed."
        )
    write(out / "curation_review.tsv", completed, list(completed[0]))
    report(
        out / "curation_summary.json",
        dict(
            research_only=True,
            reviewed=len(completed),
            limitation="Agreement against ClinVar is descriptive, not independent validation. Reviewers must evaluate disease-specific evidence and document any use of reference assertions.",
        ),
    )
    (out / "variant_assessments.md").write_text(
        "# Research-only variant assessments\n\n"
        + "\n\n".join(
            f"## {r['gene']} — ClinVar {r['variation_id']}\n- **Condition** {r.get('reviewed_condition') or r['condition']}\n- **Inheritance** {r['inheritance']}\n- **Research classification** {r['research_classification']}\n- **Evidence** {r['acmg_codes_and_strengths']}\n- **Rationale** {r['evidence_rationale']}\n- **PMIDs** {r['pmids']}\n- **Population AF** {r['population_af']} ({r['gnomad_version']})\n- **Gene validity** {r['gene_disease_validity_source']}\n- **Reference** {r['reference_class']}\n- **Reviewer/date** {r['reviewer']} / {r['review_date']}"
            for r in completed
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    s = p.add_subparsers(dest="mode", required=True)
    a = s.add_parser("prioritize")
    a.add_argument("--clinvar", type=Path, required=True)
    a.add_argument("--hpo", type=Path, required=True)
    a.add_argument("--cases", type=Path, required=True)
    a.add_argument("--limit", type=int, default=5)
    a.add_argument("--out", type=Path, required=True)
    b = s.add_parser("audit")
    b.add_argument("--queue", type=Path, required=True)
    b.add_argument("--curation", type=Path, required=True)
    b.add_argument("--out", type=Path, required=True)
    x = p.parse_args()
    x.out.mkdir(parents=True, exist_ok=True)
    if x.mode == "prioritize":
        prioritize(x.clinvar, x.hpo, x.cases, x.out, x.limit)
    else:
        audit(x.queue, x.curation, x.out)
