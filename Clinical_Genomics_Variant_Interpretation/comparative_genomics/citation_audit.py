"""Evaluate phenotype overlap using only explicit non-case PMID annotations."""

import csv
import re
import statistics
from collections import defaultdict

from ranking import P, save

cases = list(csv.DictReader(open(P / "rare/cases.tsv"), delimiter="\t"))
genes = {c["gene"] for c in cases}
diseases = defaultdict(set)
for r in csv.DictReader(
    open("rare_disease_variant_curation/data/genes_to_phenotype.txt"), delimiter="\t"
):
    if r["gene_symbol"] in genes:
        diseases[r["disease_id"]].add(r["gene_symbol"])
ann = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
with (P / "data/phenotype.hpoa").open() as f:
    for r in csv.DictReader(
        (line for line in f if not line.startswith("#")), delimiter="\t"
    ):
        if (
            r["qualifier"] == "NOT"
            or r["frequency"] == "HP:0040285"
            or re.fullmatch(r"0(?:\.0+)?%|0/\d+", r["frequency"])
        ):
            continue
        for gene in diseases.get(r["database_id"], ()):
            ann[gene][r["database_id"]][r["hpo_id"]].update(
                re.findall(r"PMID:(\d+)", r["reference"])
            )
result = []
for case in cases:
    source = re.search(r"/PMID_(\d+)", case["packet"]).group(1)
    terms = set(case["terms"].split(";"))
    scores = {}
    for gene in genes:
        scores[gene] = max(
            (
                sum(bool(lookup.get(t, set()) - {source}) for t in terms)
                for lookup in ann[gene].values()
            ),
            default=0,
        )
    value = scores[case["gene"]]
    lo = 1 + sum(v > value for v in scores.values())
    hi = sum(v >= value for v in scores.values())
    result.append(
        dict(
            gene=case["gene"],
            case_pmid=source,
            retained_term_matches=value,
            terms=len(terms),
            best=lo,
            worst=hi,
            top1=(1 / (hi - lo + 1) if lo == 1 else 0),
            rr=statistics.mean(1 / r for r in range(lo, hi + 1)),
        )
    )
save("rare/distinct_publication.tsv", result)
print(
    "Top1",
    statistics.mean(r["top1"] for r in result),
    "MRR",
    statistics.mean(r["rr"] for r in result),
    "no matches",
    sum(r["retained_term_matches"] == 0 for r in result),
)
