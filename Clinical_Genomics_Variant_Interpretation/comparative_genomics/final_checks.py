"""Independent audit of fixed candidate denominators, evidence provenance and CNV target calls."""

import csv
import json
import statistics
import tarfile
import xml.etree.ElementTree as ET

import numpy as np
from cnv import cohort, experiment, metrics
from ranking import P, save


def read(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


cases = read(P / "rare/cases.tsv")
assert (
    len(cases)
    == len({r["gene"] for r in cases})
    == len({r["allele"] for r in cases})
    == 30
)
with tarfile.open(P / "data/packets.tgz") as archive:
    for c in cases:
        doc = json.load(archive.extractfile(c["packet"]))
        gs = [
            g
            for i in doc["interpretations"]
            if i["progressStatus"] == "SOLVED"
            for g in i["diagnosis"]["genomicInterpretations"]
            if g["interpretationStatus"] == "CAUSATIVE"
        ]
        assert len(gs) == 1
        desc = gs[0]["variantInterpretation"]["variationDescriptor"]
        assert desc["geneContext"]["symbol"] == c["gene"] and desc["vcfRecord"][
            "genomeAssembly"
        ] in ("GRCh38", "hg38")
        positive = {
            f["type"]["id"]
            for f in doc["phenotypicFeatures"]
            if not f.get("excluded") and f["type"]["id"].startswith("HP:")
        }
        assert positive == set(c["terms"].split(";")) and len(positive) >= 2
ranks = read(P / "rare/ranks.tsv")
for s in read(P / "rare/recovery.tsv"):
    subset = [
        r
        for r in ranks
        if (r["method"], r["condition"]) == (s["method"], s["condition"])
    ]
    assert len(subset) == len({r["gene"] for r in subset}) == 30
    for k in ("top1", "top5", "top10", "rr"):
        assert abs(statistics.mean(float(r[k]) for r in subset) - float(s[k])) < 1e-10
submissions = {}
for source in ("clinvar.xml", "extra.xml"):
    for v in ET.parse(P / "data" / source).getroot():
        for a in v.findall(".//ClinicalAssertion"):
            acc = a.find("ClinVarAccession")
            if acc is not None:
                submissions[acc.get("Accession") + "." + acc.get("Version", "")] = (
                    v,
                    a,
                )
for r in read(P / "rare/evidence_reviews.tsv"):
    v, a = submissions[r["scv"]]
    assert (
        v.get("VariationID"),
        v.get("VariationName"),
        a.findtext("Classification/GermlineClassification"),
    ) == (r["vid"], r["hgvs"], r["reconstruction"])
    assert a.findtext("Classification/ReviewStatus") == "reviewed by expert panel"
# Recompute predictions independently from the two locked parameter sets.
out = []
for seed in (8841, 8842, 8843, 8844):
    both = lost = new = 0
    for gc, refs, _, depth, truth in cohort(seed):
        library_scaled = np.array([r / np.median(r) for r in refs[0] + refs[1]])
        ratio = np.log2((depth / np.median(depth)) / np.median(library_scaled, axis=0))
        correct = []
        for revised in (False, True):
            z = ratio.copy()
            if revised:
                design = np.column_stack([gc, np.ones(len(gc))])
                z -= design @ np.linalg.lstsq(design, z, rcond=None)[0]
            padded = np.pad(z, 2, mode="edge")
            smoothed = np.array([np.median(padded[i : i + 5]) for i in range(len(z))])
            call = np.where(
                smoothed < -0.5,
                1,
                np.where(smoothed > (0.5 if revised else 0.35), 3, 2),
            )
            if revised:
                for group in np.split(
                    np.arange(len(call)), np.nonzero(np.diff(call))[0] + 1
                ):
                    if len(group) < 3:
                        call[group] = 2
            correct.append((call == truth) & (truth != 2))
        a, b = correct
        both += sum(a & b)
        lost += sum(a & ~b)
        new += sum(~a & b)
    out.append(
        dict(seed=seed, shared_true_targets=both, comparator_only=lost, locked_only=new)
    )
save("sequencing/audit_tp_overlap.tsv", out)

seeds = []
for seed in (8841, 8842, 8843, 8844):
    counts = experiment(seed)
    for label, key in [
        ("comparator", (True, False, False, 1, 0.35)),
        ("locked", (True, False, True, 3, 0.5)),
    ]:
        seeds.append(dict(seed=seed, model=label, **metrics(key, counts[key])))
save("sequencing/audit_seeds.tsv", seeds)
