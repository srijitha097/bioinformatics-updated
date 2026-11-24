import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from single_gene_screening.run import posterior

DATA = Path("ngs_pipeline/data/single_gene")
OUT = Path("results/ngs/single_gene/benchmark")


def evaluate(group):

    called = group[group.call >= 0]
    affected = group[group.truth == 2]
    nonaffected = group[group.truth != 2]

    true_positive = int(
        (affected.call == 2).sum()
    )

    false_positive = int(
        (nonaffected.call == 2).sum()
    )

    correct = int(
        (group.call == group.truth).sum()
    )

    called_nonaffected = nonaffected[
        nonaffected.call >= 0
    ]

    return {
        "evaluations": int(len(group)),
        "called": int(len(called)),
        "correct": correct,
        "incorrect": int(len(called) - correct),
        "no_calls": int(len(group) - len(called)),
        "accuracy_all": float(
            correct / len(group)
        ),
        "accuracy_called": (
            float(correct / len(called))
            if len(called) else None
        ),
        "no_call_rate": float(
            1 - len(called) / len(group)
        ),
        "affected_samples": int(len(affected)),
        "true_positives": true_positive,
        "false_positives": false_positive,
        "affected_detection_rate": (
            float(true_positive / len(affected))
            if len(affected) else None
        ),
        "false_positive_rate_all": (
            float(false_positive / len(nonaffected))
            if len(nonaffected) else None
        ),
        "specificity_called": (
            float(
                (called_nonaffected.call != 2).mean()
            )
            if len(called_nonaffected) else None
        )
    }


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--confidence",
        type=float,
        default=0.95
    )

    parser.add_argument(
        "--min-molecules",
        type=int,
        default=100
    )

    args = parser.parse_args()

    if not 0.5 < args.confidence <= 1:
        raise ValueError("Invalid posterior threshold")

    if args.min_molecules < 1:
        raise ValueError("Invalid minimum molecular depth")

    counts = pd.read_csv(
        OUT / "molecular_allele_counts.csv"
    )

    manifest = pd.read_csv(
        DATA / "sample_manifest.csv"
    )

    qc = pd.read_csv(
        OUT / "alignment_qc_summary.csv"
    )

    df = counts.merge(
        manifest[[
            "sample",
            "estimated_ff",
            "molecules_per_gene"
        ]],
        on="sample",
        validate="many_to_one"
    )

    df = df.merge(
        qc[["sample", "qc_status"]],
        on="sample",
        validate="many_to_one"
    )

    assert len(df) == len(counts)

    rows = []

    for row in df.itertuples():

        for mode, paternal in [
            ("carrier_known", 1),
            ("paternal_unknown", None)
        ]:

            n = int(row.consensus_depth)
            alt = int(row.consensus_alt)

            if not 0 <= alt <= n:
                raise ValueError("Invalid allele counts")

            if row.qc_status != "PASS":

                post = [np.nan] * 3
                call = -1
                reason = "FAILED_ALIGNMENT_QC"

            elif n < args.min_molecules:

                post = [np.nan] * 3
                call = -1
                reason = "LOW_MOLECULAR_DEPTH"

            elif not 0 < row.estimated_ff <= 0.5:

                post = [np.nan] * 3
                call = -1
                reason = "INVALID_FETAL_FRACTION"

            else:

                post = posterior(
                    alt,
                    n,
                    maternal=1,
                    ff=row.estimated_ff,
                    paternal=paternal,
                    q=0.02,
                    rho=0.0005,
                    ff_sd=0.005
                )

                if not np.isclose(
                    np.sum(post), 1.0
                ):
                    raise ValueError(
                        "Invalid posterior probabilities"
                    )

                if max(post) >= args.confidence:

                    call = int(np.argmax(post))
                    reason = "CALLED"

                else:

                    call = -1
                    reason = "LOW_POSTERIOR"

            rows.append({
                "sample": row.sample,
                "gene": row.gene,
                "clinvar_id": row.clinvar_id,
                "mode": mode,
                "estimated_ff": row.estimated_ff,
                "input_depth": row.molecules_per_gene,
                "consensus_depth": n,
                "consensus_alt": alt,
                "consensus_vaf": row.consensus_vaf,
                "call": call,
                "call_reason": reason,
                "p0": float(post[0]),
                "p1": float(post[1]),
                "p2": float(post[2])
            })

    predictions = pd.DataFrame(rows)

    predictions.to_csv(
        OUT / "genotype_predictions.csv",
        index=False
    )

    # Ground truth is attached only after inference.
    truth = pd.read_csv(
        DATA / "molecular_ground_truth.csv"
    )

    results = predictions.merge(
        truth[[
            "sample",
            "gene",
            "clinvar_id",
            "fetal_genotype",
            "true_ff"
        ]],
        on=["sample", "gene", "clinvar_id"],
        validate="many_to_one"
    ).rename(
        columns={"fetal_genotype": "truth"}
    )

    assert len(results) == len(predictions)
    assert len(results) == 2 * len(truth)

    results.to_csv(
        OUT / "genotype_benchmark.csv",
        index=False
    )

    metrics = {
        "posterior_threshold": args.confidence,
        "minimum_consensus_depth": args.min_molecules,
        "total_evaluations": int(len(results))
    }

    for mode, group in results.groupby("mode"):

        metrics[mode] = evaluate(group)

        confusion = pd.crosstab(
            group.truth,
            group.call
        ).reindex(
            index=[0, 1, 2],
            columns=[-1, 0, 1, 2],
            fill_value=0
        )

        confusion.to_csv(
            OUT / f"confusion_{mode}.csv"
        )

        print(f"\nCONFUSION MATRIX: {mode}")
        print(confusion.to_string())

    with open(
        OUT / "genotype_screening_metrics.json",
        "w"
    ) as handle:

        json.dump(metrics, handle, indent=2)

    performance = []

    for keys, group in results.groupby([
        "mode", "true_ff", "input_depth"
    ]):

        mode, ff, depth = keys

        entry = evaluate(group)

        entry.update({
            "mode": mode,
            "true_ff": ff,
            "input_depth": depth
        })

        performance.append(entry)

    performance = pd.DataFrame(performance)

    performance.to_csv(
        OUT / "genotype_performance_by_ff_depth.csv",
        index=False
    )

    gene_performance = []

    for keys, group in results.groupby([
        "gene", "mode"
    ]):

        gene, mode = keys

        entry = evaluate(group)

        entry.update({
            "gene": gene,
            "mode": mode
        })

        gene_performance.append(entry)

    pd.DataFrame(gene_performance).to_csv(
        OUT / "genotype_performance_by_gene.csv",
        index=False
    )

    reasons = pd.crosstab(
        results["mode"],
        results["call_reason"]
    )

    reasons.to_csv(
        OUT / "genotype_call_reasons.csv"
    )

    plots = [
        (
            "affected_detection_rate",
            "Affected-genotype detection rate (%)",
            "Affected Genotype Detection",
            "affected_detection_performance.png"
        ),
        (
            "no_call_rate",
            "No-call rate (%)",
            "Genotype No-call Performance",
            "genotype_no_call_performance.png"
        ),
        (
            "accuracy_all",
            "Correct calls / all evaluations (%)",
            "Fetal Genotype Inference Accuracy",
            "genotype_accuracy_performance.png"
        )
    ]

    for metric, ylabel, title, filename in plots:

        fig, ax = plt.subplots(figsize=(10, 6))

        for keys, group in performance.groupby([
            "mode", "input_depth"
        ]):

            mode, depth = keys

            group = group.sort_values("true_ff")

            ax.plot(
                group.true_ff * 100,
                group[metric] * 100,
                marker="o",
                label=f"{mode}, {depth:,}"
            )

        ax.set(
            xlabel="Simulated fetal fraction (%)",
            ylabel=ylabel,
            title=title,
            ylim=(0, 105)
        )

        ax.legend(fontsize=8)
        ax.grid(alpha=0.2)

        fig.tight_layout()

        fig.savefig(
            OUT / filename,
            dpi=300
        )

        plt.close(fig)

    print("\nGENOTYPE SCREENING METRICS")
    print(json.dumps(metrics, indent=2))

    print("\nPERFORMANCE BY GENE")

    print(
        pd.DataFrame(
            gene_performance
        ).round(5).to_string(index=False)
    )

    print("\nCALL REASONS")
    print(reasons.to_string())

    print("\nAnalysis completed")
    print("Results saved to:", OUT)


if __name__ == "__main__":
    main()
