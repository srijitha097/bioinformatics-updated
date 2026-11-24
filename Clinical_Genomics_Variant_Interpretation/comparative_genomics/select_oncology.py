import csv
from pathlib import Path

P = Path("oncology_somatic_variant_interpretation/data")


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


ev = {r["evidence_id"]: r for r in read(next(P.glob("*AcceptedClinical*")))}
seen = set()
chosen = []
for a in read(next(P.glob("*AcceptedAssertion*"))):
    profile = a["molecular_profile"]
    if (
        a["assertion_type"] not in ("Predictive", "Diagnostic", "Prognostic")
        or a["is_flagged"] == "true"
    ):
        continue
    if profile in seen and not (
        profile == "BRAF V600E" and a["assertion_type"] == "Prognostic"
    ):
        continue
    if a["assertion_id"] in {"2", "3", "6", "7", "9", "39"}:
        continue
    candidates = [
        ev[i]
        for i in a["evidence_item_ids"].split(",")
        if i in ev
        and ev[i]["source_type"] == "PubMed"
        and ev[i]["is_flagged"] == "false"
        and ev[i]["evidence_type"] == a["assertion_type"]
    ]
    if not candidates:
        continue
    r = candidates[0]
    seen.add(profile)
    chosen.append(
        [
            a["assertion_id"],
            profile,
            a["disease"],
            a["therapies"],
            a["assertion_type"],
            a["significance"],
            r["citation_id"],
            r["evidence_id"],
        ]
    )
    if len(chosen) == 23:
        break
with open("comparative_genomics/oncology/contexts.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow("aid profile disease therapy context significance pmid eid".split())
    w.writerows(chosen)
