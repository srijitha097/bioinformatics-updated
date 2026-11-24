import csv
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ranking import P


def read(name):
    return list(csv.DictReader(open(P / name), delimiter="\t"))


def finish(ax, path):
    ax.spines[["top", "right"]].set_visible(False)
    ax.figure.tight_layout()
    ax.figure.savefig(P / path, dpi=180)
    ax.figure.savefig((P / path).with_suffix(".pdf"), metadata={"CreationDate": None})
    plt.close(ax.figure)


fig, ax = plt.subplots(figsize=(6, 3.5))
rr = read("rare/recovery.tsv")
for method in ("phenotype", "combined", "frequency"):
    s = [r for r in rr if r["method"] == method]
    ax.plot(
        [r["condition"] for r in s],
        [float(r["rr"]) for r in s],
        marker="o",
        label=method,
    )
ax.set(
    ylabel="Mean reciprocal rank",
    ylim=(0, 1.05),
    title="30 published cases ranked within a 30-gene panel",
)
ax.legend()
finish(ax, "rare/ranking.png")
fig, ax = plt.subplots(figsize=(6, 3.5))
for kind in ("SNV", "indel"):
    subset = [
        r
        for r in read("sequencing/metrics.tsv")
        if r["GQ"] == "0" and r["scope"] == "exome" and r["kind"] == kind
    ]
    ax.plot(
        [int(r["DP"]) for r in subset],
        [float(r["recall"]) for r in subset],
        marker="o",
        label=f"{kind} (n={subset[0]['truth']} truth genotypes)",
    )
ax.set(
    xlabel="Minimum query read depth (×)",
    ylabel="Recall",
    ylim=(0, 1),
    title="NA12878 exact genotype recovery",
)
ax.legend()
finish(ax, "sequencing/depth.png")
fig, ax = plt.subplots(figsize=(6, 3.5))

c = Counter(r["assessment"] for r in read("oncology/comparison.tsv"))
ax.bar(list(c), list(c.values()), color="#236f86")
ax.set(title="Initial literature assessments across 23 contexts", ylabel="Contexts")
finish(ax, "oncology/tiers.png")
