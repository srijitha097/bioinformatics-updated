import csv
from collections import Counter

from ranking import P, save

with open(P / "data/ccle.seg") as f:
    segments = [
        s
        for s in csv.DictReader(f, delimiter="\t")
        if s["ID"] == "HCC1143_BREAST" and s["chrom"] == "20"
    ]
with open("exome_gene_panel_cnv_analysis/results/real_copy_ratios.tsv") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
out = []
for r in rows:
    hits = [
        s
        for s in segments
        if int(s["loc.start"]) <= int(r["start"]) + 1
        and int(s["loc.end"]) >= int(r["end"])
    ]
    if len(hits) == 1:
        out.append(
            dict(
                r,
                array_start=hits[0]["loc.start"],
                array_end=hits[0]["loc.end"],
                array_mean=hits[0]["seg.mean"],
                markers=hits[0]["num.mark"],
            )
        )
save("sequencing/orthogonal_array.tsv", out)


usable_counts = Counter(
    (r["array_start"], r["array_end"], float(r["array_mean"]) > 0)
    for r in out
    if r["centered_log2_ratio"] != ""
)
save(
    "sequencing/array_overlap_summary.tsv",
    [
        dict(
            capture_targets=len(rows),
            matched_targets=len(out),
            usable_targets=sum(usable_counts.values()),
            usable_segments=len(usable_counts),
            positive_segments=sum(key[2] for key in usable_counts),
            usable_targets_per_segment=",".join(
                map(str, sorted(usable_counts.values(), reverse=True))
            ),
        )
    ],
)
