"""Extract accepted CIViC clinical assertions; retain source tier and disease context."""

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import report, rows, sha, write


def norm(s):
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")


def pick(r, *keys):
    return next((r[k] for k in keys if r.get(k)), "")


def extract(assertions, out, accepted_only=False):
    selected = []
    missing = 0
    for raw in rows(assertions):
        r = {norm(k): v for k, v in raw.items()}
        status = pick(r, "assertion_status", "status").lower()
        if status != "accepted" and (status or not accepted_only):
            continue
        tier = pick(
            r,
            "amp_asco_cap_tier",
            "amp_asco_cap_category",
            "amp_asco_cap_classification",
            "amp_category",
            "amp_tier",
        )
        if not tier:
            tier = next(
                (
                    v
                    for k, v in r.items()
                    if "amp" in k and ("tier" in k or "category" in k) and v
                ),
                "",
            )
        if not tier or tier.lower() in {"n/a", "na", "none", "not applicable"}:
            missing += 1
            continue
        aid = pick(r, "assertion_id", "id")
        disease = pick(r, "disease", "disease_name")
        profile = pick(
            r, "molecular_profile", "molecular_profile_name", "variant", "variant_name"
        )
        kind = pick(r, "assertion_type", "evidence_type", "type")
        if not aid or not disease or not profile:
            continue
        selected.append(
            dict(
                assertion_id=aid,
                molecular_profile=profile,
                disease=disease,
                assertion_type=kind,
                source_tier=tier,
                direction=pick(r, "assertion_direction"),
                significance=pick(r, "significance", "clinical_significance"),
                therapies=pick(r, "therapies", "therapy", "drug", "drugs"),
                evidence_ids=pick(
                    r, "evidence_item_ids", "evidence_items", "evidence_ids"
                ),
                source_url=pick(r, "assertion_civic_url", "assertion_url")
                or f"https://civicdb.org/assertions/{aid}/summary",
            )
        )
    if not selected:
        raise ValueError(
            "No accepted tiered assertions found. Inspect TSV header and release documentation before proceeding."
        )
    write(out / "accepted_assertions.tsv", selected, list(selected[0]))
    report(
        out / "summary.json",
        dict(
            research_only=True,
            accepted_tiered_assertions=len(selected),
            accepted_without_tier=missing,
            unique_diseases=len({r["disease"] for r in selected}),
            tiers=dict(Counter(r["source_tier"] for r in selected)),
            source_sha256=sha(assertions),
            accepted_only_export=accepted_only,
            limitation="Imported source tiers, not independently assigned AMP tiers. Molecular profiles may contain multiple variants.",
        ),
    )
    # A separate manual review sheet, never auto-populated with a predicted tier.
    write(
        out / "review_template.tsv",
        [
            {
                k: r[k]
                for k in (
                    "assertion_id",
                    "molecular_profile",
                    "disease",
                    "assertion_type",
                    "source_url",
                )
            }
            for r in selected[:30]
        ],
        [
            "assertion_id",
            "molecular_profile",
            "disease",
            "assertion_type",
            "source_url",
            "reviewed_pmids",
            "evidence_summary",
            "research_tier",
            "reviewer",
            "review_date",
        ],
    )


def audit(source, review, out, evidence=None):
    accepted_evidence = {}
    if evidence is not None:
        accepted_evidence = {
            row["evidence_id"]: row
            for row in rows(evidence)
            if row["evidence_status"].lower() == "accepted"
        }
    ref = {r["assertion_id"]: r for r in rows(source)}
    done = []
    seen = set()
    for r in rows(review):
        if not r.get("research_tier"):
            continue
        k = r["assertion_id"]
        required = ["reviewed_pmids", "evidence_summary", "reviewer", "review_date"]
        if k not in ref or k in seen or any(not r.get(x) for x in required):
            raise ValueError(f"Invalid/incomplete review {k}")
        if any(
            r[x] != ref[k][x]
            for x in ["molecular_profile", "disease", "assertion_type"]
        ):
            raise ValueError("Context mismatch")
        if evidence is not None:
            reviewed_ids = set(
                filter(None, r.get("reviewed_evidence_ids", "").split(","))
            )
            linked_ids = set(ref[k]["evidence_ids"].split(","))
            if (
                not reviewed_ids
                or reviewed_ids - linked_ids
                or reviewed_ids - accepted_evidence.keys()
            ):
                raise ValueError(f"Invalid or unlinked evidence for assertion {k}")
            pmids = {
                accepted_evidence[eid]["citation_id"]
                for eid in reviewed_ids
                if accepted_evidence[eid]["source_type"].lower() == "pubmed"
            }
            if pmids != set(re.split(r"[,; ]+", r["reviewed_pmids"])):
                raise ValueError(
                    f"Reviewed PMIDs do not match evidence for assertion {k}"
                )
        seen.add(k)
        done.append(dict(r, source_tier=ref[k]["source_tier"]))
    if not done:
        raise ValueError("No manually reviewed assertions")
    write(out / "reviewed_assertions.tsv", done, list(done[0]))
    report(
        out / "review_summary.json",
        dict(
            research_only=True,
            reviewed=len(done),
            evidence_sha256=sha(evidence) if evidence is not None else None,
            limitation="Compare source and research tiers manually; source tier notation may include evidence level. No clinical validation.",
        ),
    )
    (out / "oncology_assessments.md").write_text(
        "# Research-only oncology assessments\n\n"
        + "\n\n".join(
            f"## {r['molecular_profile']} — {r['disease']}\n- **Context** {r['assertion_type']}\n- **Research tier** {r['research_tier']}\n- **Source tier** {r['source_tier']}\n- **Evidence summary** {r['evidence_summary']}\n- **PMIDs** {r['reviewed_pmids']}\n- **Evidence links** {r.get('evidence_source', '')}\n- **CIViC assertion** {r['source_url']}"
            for r in done
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    s = p.add_subparsers(dest="mode", required=True)
    a = s.add_parser("extract")
    a.add_argument("--assertions", type=Path, required=True)
    a.add_argument("--accepted-only-export", action="store_true")
    a.add_argument("--out", type=Path, required=True)
    b = s.add_parser("audit")
    b.add_argument("--source", type=Path, required=True)
    b.add_argument("--review", type=Path, required=True)
    b.add_argument("--evidence", type=Path)
    b.add_argument("--out", type=Path, required=True)
    x = p.parse_args()
    x.out.mkdir(parents=True, exist_ok=True)
    if x.mode == "extract":
        extract(x.assertions, x.out, x.accepted_only_export)
    else:
        audit(x.source, x.review, x.out, x.evidence)
