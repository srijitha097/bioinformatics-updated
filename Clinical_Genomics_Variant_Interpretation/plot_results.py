"""Render publication-style research figures from the saved workflow outputs."""

import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import rows, write

ROOT = Path(__file__).resolve().parent
RARE = ROOT / "rare_disease_variant_curation" / "results"
CIVIC = ROOT / "oncology_somatic_variant_interpretation" / "results"
NGS = ROOT / "exome_gene_panel_cnv_analysis" / "results"
BLUE = "#19699c"
RUST = "#b34b25"
GREY = "#657785"


def save(fig, path):
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(
        path.with_suffix(".pdf"), bbox_inches="tight", metadata={"CreationDate": None}
    )
    plt.close(fig)


def rare_plot():
    grouped = {}
    counts = Counter()
    for row in rows(RARE / "candidate_queue.tsv"):
        key = row["case_id"], row["gene"]
        grouped[key] = row
        counts[key] += 1
    records = []
    for key, row in grouped.items():
        records.append(
            dict(
                case_id=key[0],
                gene=key[1],
                hpo_overlap=float(row["hpo_overlap"]),
                disease_ids=row["hpo_best_disease_ids"],
                candidate_variants=counts[key],
            )
        )
    write(RARE / "phenotype_gene_scores.tsv", records, list(records[0]))
    fig, ax = plt.subplots(figsize=(7, 4.4), layout="constrained")
    labels = [f"{r['case_id']} · {r['gene']}" for r in records]
    ax.barh(labels, [r["hpo_overlap"] for r in records], color=BLUE)
    ax.invert_yaxis()
    ax.set(
        xlim=(0, 1.05),
        xlabel="Fraction of scenario HPO terms matched",
        title="Highest single-disease HPO overlap by gene",
    )
    save(fig, RARE / "phenotype_gene_scores.png")


def civic_plot():
    records = list(rows(CIVIC / "accepted_assertions.tsv"))
    categories = sorted({r["assertion_type"] for r in records})
    tiers = ["Tier I - Level A", "Tier II - Level C"]
    counts = Counter((r["assertion_type"], r["source_tier"]) for r in records)
    table = [
        dict(assertion_type=kind, tier=tier, count=count)
        for (kind, tier), count in sorted(counts.items())
    ]
    write(CIVIC / "tier_by_assertion_type.tsv", table, list(table[0]))
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.2), sharey=True, layout="constrained")
    for ax, tier, color in zip(axes, tiers, [BLUE, RUST]):
        values = [counts[category, tier] for category in categories]
        bars = ax.barh(categories, values, color=color)
        ax.bar_label(bars, padding=4)
        ax.set(xlim=(0, max(values) * 1.25), title=tier)
    axes[1].set_xticks(
        range(max(counts[category, tiers[1]] for category in categories) + 1)
    )
    fig.supxlabel("Accepted assertions (scales differ by tier)")
    fig.suptitle("CIViC source tiers · 1 November 2026")
    save(fig, CIVIC / "tier_by_assertion_type.png")


def sequencing_plots():
    summary = []
    for scope in ("wes", "panel"):
        data = json.loads((NGS / f"{scope}.json").read_text())
        counts = data["true_positive"], data["false_positive"], data["false_negative"]
        row = dict(
            scope=scope,
            tp=counts[0],
            fp=counts[1],
            fn=counts[2],
            precision=data["precision"],
            recall=data["recall"],
            f1=data["f1"],
            confident_target_bases=data["confident_target_bases"],
            target_bases=data["target_bases"],
        )
        summary.append(row)
    write(NGS / "benchmark_metrics.tsv", summary, list(summary[0]))
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True, layout="constrained")
    for ax, label, record in zip(axes, ("ARUP exome", "10-gene exome subset"), summary):
        values = [record["precision"], record["recall"]]
        bars = ax.bar(["Precision", "Recall"], values, color=[BLUE, GREY])
        ax.bar_label(bars, [f"{v:.1%}" for v in values], padding=3)
        ax.set(
            ylim=(0, 1.12),
            title=f"{label}\n{record['tp']} TP · {record['fp']} FP · {record['fn']} FN",
        )
    axes[0].set_ylabel("Fraction")
    fig.suptitle("NA12878 exact SNV and genotype benchmark")
    save(fig, NGS / "benchmark_metrics.png")

    records = [
        r for r in rows(NGS / "real_copy_ratios.tsv") if r["centered_log2_ratio"]
    ]
    fig, ax = plt.subplots(figsize=(7, 3.4), layout="constrained")
    ax.scatter(
        [(int(r["start"]) + int(r["end"])) / 2e6 for r in records],
        [float(r["centered_log2_ratio"]) for r in records],
        s=13,
        color=BLUE,
    )
    ax.axhline(0, color=GREY, linewidth=1)
    ax.set(
        xlabel="Target midpoint on chromosome 20 (Mb)",
        ylabel="Centered log₂(tumor / normal)",
        title=f"HCC1143 exome slices · {len(records)} targets with ≥10× mean depth in both",
    )
    save(fig, NGS / "real_copy_ratios.png")

    records = [
        r
        for r in rows(NGS / "synthetic" / "synthetic_targets.tsv")
        if r["sample"] in {"1", "2", "10", "11"}
    ]
    fig, axes = plt.subplots(
        2, 2, figsize=(8, 5.1), sharex=True, sharey=True, layout="constrained"
    )
    for ax, sample in zip(axes.flat, ("1", "2", "10", "11")):
        rs = [r for r in records if r["sample"] == sample]
        ax.plot(
            [int(r["target"]) for r in rs],
            [float(r["normalized_log2_ratio"]) for r in rs],
            color=BLUE,
            linewidth=1.2,
        )
        affected = [int(r["target"]) for r in rs if r["true_cn"] != "2"]
        ax.axvspan(min(affected) - 0.5, max(affected) + 0.5, color=RUST, alpha=0.22)
        for threshold in (-0.5, 0.35):
            ax.axhline(threshold, color=GREY, linestyle=":", linewidth=1)
        ax.set_title(
            f"Sample {sample} · {len(affected)}-target {'loss' if sample in {'1', '10'} else 'gain'}"
        )
    fig.supxlabel("Target index")
    fig.supylabel("Smoothed log₂ depth ratio")
    fig.suptitle("Synthetic CNV screen · shaded intervals show simulated truth")
    save(fig, NGS / "synthetic" / "cnv_examples.png")


if __name__ == "__main__":
    plt.rcParams.update(
        {"axes.spines.top": False, "axes.spines.right": False, "font.size": 10}
    )
    rare_plot()
    civic_plot()
    sequencing_plots()
