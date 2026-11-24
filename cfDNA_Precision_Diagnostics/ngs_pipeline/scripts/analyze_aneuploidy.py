import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("results/ngs/aneuploidy/benchmark")
DATA = Path("ngs_pipeline/data/aneuploidy")

CONTROL_CHROMS = ["chr1", "chr2", "chr3"]
TARGET_CHROMS = ["chr13", "chr18", "chr21"]


def gc_correct(group):

    group = group.copy()

    controls = group[
        group.chromosome.isin(CONTROL_CHROMS)
    ]

    x = controls.gc_content.to_numpy() - 0.45

    y = np.log(
        controls.fragment_count.to_numpy() + 0.5
    )

    coefficients = np.polyfit(x, y, deg=2)

    all_gc = group.gc_content.to_numpy() - 0.45

    fitted = np.polyval(coefficients, all_gc)

    bias = np.exp(
        fitted - np.median(fitted)
    )

    group["estimated_gc_bias"] = bias

    group["corrected_count"] = (
        group.fragment_count / bias
    )

    return group


def classify(positives, qc_status, ff, min_ff):

    if qc_status != "PASS":
        return "NO_CALL_QC"

    if ff < min_ff:
        return "NO_CALL_LOW_FF"

    if len(positives) > 1:
        return "NO_CALL_MULTIPLE"

    if not positives:
        return "euploid"

    return positives[0].replace("chr", "trisomy")


def evaluate(results, column):

    called = results[
        ~results[column].str.startswith("NO_CALL")
    ]

    positive = called[
        called.truth != "euploid"
    ]

    negative = called[
        called.truth == "euploid"
    ]

    correct = int(
        (called[column] == called.truth).sum()
    )

    return {
        "evaluation_samples": int(len(results)),
        "called_samples": int(len(called)),
        "no_calls": int(len(results) - len(called)),
        "no_call_rate": float(
            1 - len(called) / len(results)
        ),
        "correct_calls": correct,
        "accuracy_all": float(
            correct / len(results)
        ),
        "accuracy_called": float(
            correct / len(called)
        ) if len(called) else None,
        "exact_trisomy_sensitivity_called": float(
            (positive[column] == positive.truth).mean()
        ) if len(positive) else None,
        "specificity_called": float(
            (negative[column] == "euploid").mean()
        ) if len(negative) else None
    }


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--z-threshold", type=float, default=3.0
    )

    parser.add_argument(
        "--min-ff", type=float, default=0.04
    )

    args = parser.parse_args()

    counts = pd.read_csv(
        OUT / "window_fragment_counts.csv"
    )

    manifest = pd.read_csv(
        DATA / "sample_manifest.csv"
    )

    alignment_qc = pd.read_csv(
        OUT / "alignment_qc_summary.csv"
    )

    print("Applying sample-specific GC correction")

    corrected = pd.concat([
        gc_correct(group)
        for _, group in counts.groupby("sample")
    ], ignore_index=True)

    corrected.to_csv(
        OUT / "gc_corrected_window_counts.csv",
        index=False
    )

    dosage_rows = []

    for sample, group in corrected.groupby("sample"):

        raw = group.groupby(
            "chromosome"
        ).fragment_count.sum()

        adjusted = group.groupby(
            "chromosome"
        ).corrected_count.sum()

        raw_control = raw[CONTROL_CHROMS].mean()

        adjusted_control = adjusted[
            CONTROL_CHROMS
        ].mean()

        for chrom in TARGET_CHROMS:

            dosage_rows.append({
                "sample": sample,
                "chromosome": chrom,
                "raw_dosage":
                    raw[chrom] / raw_control,
                "corrected_dosage":
                    adjusted[chrom] / adjusted_control
            })

    dosage = pd.DataFrame(dosage_rows)

    dosage = dosage.merge(
        manifest[[
            "sample",
            "cohort",
            "condition",
            "trisomy",
            "true_ff"
        ]],
        on="sample",
        validate="many_to_one"
    )

    dosage = dosage.merge(
        alignment_qc[[
            "sample",
            "qc_status"
        ]],
        on="sample",
        validate="many_to_one"
    )

    calibration = dosage[
        (dosage.cohort == "calibration")
        & (dosage.qc_status == "PASS")
    ]

    baseline = calibration.groupby(
        "chromosome"
    ).agg(
        raw_mean=("raw_dosage", "mean"),
        raw_sd=("raw_dosage", "std"),
        corrected_mean=("corrected_dosage", "mean"),
        corrected_sd=("corrected_dosage", "std"),
        calibration_samples=("sample", "size")
    ).reset_index()

    if (
        len(baseline) != 3
        or (baseline.calibration_samples < 10).any()
        or (baseline.raw_sd <= 0).any()
        or (baseline.corrected_sd <= 0).any()
    ):
        raise RuntimeError(
            "Insufficient calibration data"
        )

    baseline["raw_cv_pct"] = (
        baseline.raw_sd / baseline.raw_mean * 100
    )

    baseline["corrected_cv_pct"] = (
        baseline.corrected_sd
        / baseline.corrected_mean * 100
    )

    baseline.to_csv(
        OUT / "euploid_baseline.csv",
        index=False
    )

    dosage = dosage.merge(
        baseline,
        on="chromosome",
        validate="many_to_one"
    )

    dosage["raw_z"] = (
        dosage.raw_dosage - dosage.raw_mean
    ) / dosage.raw_sd

    dosage["corrected_z"] = (
        dosage.corrected_dosage
        - dosage.corrected_mean
    ) / dosage.corrected_sd

    dosage.to_csv(
        OUT / "chromosome_dosage_zscores.csv",
        index=False
    )

    evaluation = dosage[
        dosage.cohort == "evaluation"
    ]

    calls = []

    for sample, group in evaluation.groupby("sample"):

        row = group.iloc[0]

        raw_positive = group[
            group.raw_z >= args.z_threshold
        ].chromosome.tolist()

        corrected_positive = group[
            group.corrected_z >= args.z_threshold
        ].chromosome.tolist()

        raw_call = classify(
            raw_positive,
            row.qc_status,
            row.true_ff,
            args.min_ff
        )

        corrected_call = classify(
            corrected_positive,
            row.qc_status,
            row.true_ff,
            args.min_ff
        )

        scores = group.set_index(
            "chromosome"
        ).corrected_z

        calls.append({
            "sample": sample,
            "truth": row.condition,
            "true_ff": row.true_ff,
            "raw_call": raw_call,
            "corrected_call": corrected_call,
            "z13": scores["chr13"],
            "z18": scores["chr18"],
            "z21": scores["chr21"]
        })

    results = pd.DataFrame(calls)

    results.to_csv(
        OUT / "screening_calls.csv",
        index=False
    )

    metrics = {
        "z_threshold": args.z_threshold,
        "minimum_simulated_ff": args.min_ff,
        "raw": evaluate(results, "raw_call"),
        "gc_corrected": evaluate(
            results, "corrected_call"
        )
    }

    with open(
        OUT / "screening_metrics.json", "w"
    ) as handle:

        json.dump(metrics, handle, indent=2)

    results["called"] = (
        ~results.corrected_call.str.startswith(
            "NO_CALL"
        )
    )

    results["correct"] = (
        results.corrected_call == results.truth
    )

    performance = results.groupby(
        ["true_ff", "truth"]
    ).agg(
        samples=("sample", "size"),
        called=("called", "sum"),
        correct=("correct", "sum")
    ).reset_index()

    performance["accuracy_all"] = (
        performance.correct
        / performance.samples
    )

    performance["accuracy_called"] = (
        performance.correct
        / performance.called.replace(0, np.nan)
    )

    performance.to_csv(
        OUT / "performance_by_ff.csv",
        index=False
    )

    confusion = pd.crosstab(
        results.truth,
        results.corrected_call
    )

    confusion.to_csv(
        OUT / "screening_confusion_matrix.csv"
    )

    # Screening performance by fetal fraction.

    fig, ax = plt.subplots(figsize=(9, 6))

    for condition, group in performance.groupby("truth"):

        ax.plot(
            group.true_ff * 100,
            group.accuracy_all * 100,
            marker="o",
            label=condition
        )

    ax.set(
        xlabel="Simulated fetal fraction (%)",
        ylabel="Correct calls / all samples (%)",
        title="Aneuploidy Screening Performance",
        ylim=(0, 105)
    )

    ax.legend()
    ax.grid(alpha=0.2)

    fig.tight_layout()

    fig.savefig(
        OUT / "aneuploidy_screening_performance.png",
        dpi=300
    )

    plt.close(fig)

    # Individual chromosome Z-score visualizations.

    for chrom in [13, 18, 21]:

        fig, ax = plt.subplots(figsize=(9, 6))

        for condition, group in results.groupby("truth"):

            ax.scatter(
                group.true_ff * 100,
                group[f"z{chrom}"],
                label=condition,
                alpha=0.8
            )

        ax.axhline(
            args.z_threshold,
            linestyle="--",
            label="Screening threshold"
        )

        ax.set(
            xlabel="Simulated fetal fraction (%)",
            ylabel=f"Chromosome {chrom} Z-score",
            title=f"Chromosome {chrom} Dosage Screening"
        )

        ax.legend()
        ax.grid(alpha=0.2)

        fig.tight_layout()

        fig.savefig(
            OUT / f"chromosome{chrom}_zscores.png",
            dpi=300
        )

        plt.close(fig)

    # GC correction baseline variability comparison.

    fig, ax = plt.subplots(figsize=(8, 5))

    x = np.arange(len(baseline))

    ax.bar(
        x - 0.18,
        baseline.raw_cv_pct,
        width=0.36,
        label="Raw counts"
    )

    ax.bar(
        x + 0.18,
        baseline.corrected_cv_pct,
        width=0.36,
        label="GC-corrected counts"
    )

    ax.set_xticks(x)
    ax.set_xticklabels(baseline.chromosome)

    ax.set(
        ylabel="Coefficient of variation (%)",
        title="Euploid Baseline Dosage Variability"
    )

    ax.legend()
    ax.grid(axis="y", alpha=0.2)

    fig.tight_layout()

    fig.savefig(
        OUT / "gc_correction_comparison.png",
        dpi=300
    )

    plt.close(fig)

    print("\nEUPLOID CALIBRATION BASELINE")
    print(baseline.round(6).to_string(index=False))

    print("\nSCREENING METRICS")
    print(json.dumps(metrics, indent=2))

    print("\nCONFUSION MATRIX")
    print(confusion.to_string())

    print("\nPERFORMANCE BY FETAL FRACTION")
    print(performance.to_string(index=False))

    print("\nAnalysis completed")
    print("Results saved to:", OUT)


if __name__ == "__main__":
    main()
