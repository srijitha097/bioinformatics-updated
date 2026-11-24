"""Rank a fixed panel of literature alleles for 30 published cases."""

import csv
import json
import math
import random
import statistics
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

P = Path("comparative_genomics")


def save(name, records):
    with (P / name).open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]), delimiter="\t")
        w.writeheader()
        w.writerows(records)


def select():
    genes = {}
    with tarfile.open(P / "data/packets.tgz") as archive:
        for member in archive:
            if "/phenopackets/PMID_" not in member.name or not member.name.endswith(
                ".json"
            ):
                continue
            doc = json.load(archive.extractfile(member))
            terms = sorted(
                {
                    x["type"]["id"]
                    for x in doc.get("phenotypicFeatures", [])
                    if not x.get("excluded") and x["type"]["id"].startswith("HP:")
                }
            )
            calls = [
                g
                for i in doc.get("interpretations", [])
                if i.get("progressStatus") == "SOLVED"
                for g in i.get("diagnosis", {}).get("genomicInterpretations", [])
                if g.get("interpretationStatus") == "CAUSATIVE"
            ]
            if len(calls) != 1 or len(terms) < 2:
                continue
            v = calls[0].get("variantInterpretation", {}).get("variationDescriptor", {})
            rec = v.get("vcfRecord", {})
            gene = v.get("geneContext", {}).get("symbol")
            if (
                not gene
                or rec.get("genomeAssembly") not in ("hg38", "GRCh38")
                or not rec.get("ref")
                or not rec.get("alt")
                or max(len(rec["ref"]), len(rec["alt"])) > 40
            ):
                continue
            allele = "-".join(
                [
                    rec["chrom"].removeprefix("chr"),
                    str(rec["pos"]),
                    rec["ref"],
                    rec["alt"],
                ]
            )
            genes.setdefault(
                gene,
                dict(
                    gene=gene, allele=allele, terms=";".join(terms), packet=member.name
                ),
            )
    names = sorted(genes)
    random.Random(181).shuffle(names)
    preferred = "GJB2 CFTR HBB PAH DMD BRCA1 BRCA2".split()
    names = [n for n in preferred if n in genes] + [
        n for n in names if n not in preferred
    ]
    cases = [genes[n] for n in names[:30]]
    save("rare/cases.tsv", cases)
    return cases


def rank(cases):
    hpo = defaultdict(lambda: defaultdict(set))
    for r in csv.DictReader(
        open("rare_disease_variant_curation/data/genes_to_phenotype.txt"),
        delimiter="\t",
    ):
        if r["frequency"] in ("HP:0040285", "0%") or r["frequency"].startswith("0/"):
            continue
        hpo[r["gene_symbol"]][r["disease_id"]].add(r["hpo_id"])
    af = json.loads((P / "data/af.json").read_text())
    allterms = sorted({t for c in cases for t in c["terms"].split(";")})
    results = []
    for i, c in enumerate(cases):
        terms = set(c["terms"].split(";"))
        rng = random.Random(i + 456)
        for condition, observed in [
            ("original", terms),
            ("half missing", set(rng.sample(sorted(terms), max(1, len(terms) // 2)))),
            ("three added", terms | set(rng.sample(sorted(set(allterms) - terms), 3))),
        ]:
            scores = defaultdict(list)
            for j, v in enumerate(cases):
                phenotype = max(
                    (
                        len(observed & ts) / len(observed)
                        for ts in hpo[v["gene"]].values()
                    ),
                    default=0,
                )
                frequency = af[v["allele"]]
                rarity = (
                    0.5
                    if frequency is None
                    else min(1, -math.log10(max(frequency, 1e-8)) / 8)
                )
                for method, score in [
                    ("phenotype", phenotype),
                    ("frequency", rarity),
                    ("combined", phenotype + 0.25 * rarity),
                ]:
                    scores[method].append(score)
            for method, s in scores.items():
                lo = 1 + sum(x > s[i] for x in s)
                hi = sum(x >= s[i] for x in s)
                ranks = range(lo, hi + 1)
                results.append(
                    dict(
                        gene=c["gene"],
                        condition=condition,
                        method=method,
                        best=lo,
                        worst=hi,
                        **{
                            f"top{k}": statistics.mean(r <= k for r in ranks)
                            for k in (1, 5, 10)
                        },
                        rr=statistics.mean(1 / r for r in ranks),
                    )
                )
    save("rare/ranks.tsv", results)
    summaries = []
    for key in dict.fromkeys((r["condition"], r["method"]) for r in results):
        subset = [r for r in results if (r["condition"], r["method"]) == key]
        summaries.append(
            dict(
                condition=key[0],
                method=key[1],
                n=len(subset),
                **{
                    m: statistics.mean(r[m] for r in subset)
                    for m in ("top1", "top5", "top10", "rr")
                },
            )
        )
    save("rare/recovery.tsv", summaries)


if __name__ == "__main__":
    cases = select()
    if "--select" not in sys.argv:
        rank(cases)
