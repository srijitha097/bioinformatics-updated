from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("results/ngs/aneuploidy/benchmark")

DEST = OUT / "threshold_analysis"
DEST.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(
    OUT / "screening_calls.csv"
)

rows = []

for threshold in np.arange(1.5, 4.01, 0.25):

    calls = []

    for row in df.itertuples():

        scores = {
            "trisomy13": row.z13,
            "trisomy18": row.z18,
            "trisomy21": row.z21
        }

        positive = [
            chrom for chrom, score in scores.items()
            if score >= threshold
        ]

        if row.corrected_call == "NO_CALL_QC":
            call = "NO_CALL_QC"

        elif row.true_ff < 0.04:
            call = "NO_CALL_LOW_FF"

        elif len(positive) > 1:
            call = "NO_CALL_MULTIPLE"

        elif len(positive) == 1:
            call = positive[0]

        else:
            call = "euploid"

        calls.append({
            "truth": row.truth,
            "call": call
        })

    result = pd.DataFrame(calls)

    called = result[
        ~result.call.str.startswith("NO_CALL")
    ]

    positive = called[
        called.truth != "euploid"
    ]

    negative = called[
        called.truth == "euploid"
    ]

    sensitivity = (
        (positive.call == positive.truth).mean()
        if len(positive) else np.nan
    )

    specificity = (
        (negative.call == "euploid").mean()
        if len(negative) else np.nan
    )

    correct = (
        called.call == called.truth
    ).sum()

    rows.append({
        "z_threshold": threshold,
        "called_samples": len(called),
        "correct_calls": correct,
        "incorrect_calls": len(called) - correct,
        "no_calls": len(result) - len(called),
        "sensitivity_pct": sensitivity * 100,
        "specificity_pct": specificity * 100,
        "balanced_accuracy_pct":
            (sensitivity + specificity) / 2 * 100,
        "accuracy_called_pct":
            correct / len(called) * 100
            if len(called) else np.nan,
        "accuracy_all_pct":
            correct / len(result) * 100,
        "no_call_pct":
            (1 - len(called) / len(result)) * 100
    })

summary = pd.DataFrame(rows)

summary.to_csv(
    DEST / "threshold_comparison.csv",
    index=False
)

print(summary.round(4).to_string(index=False))

fig, ax = plt.subplots(figsize=(10, 6))

for column, label in [
    ("sensitivity_pct", "Trisomy sensitivity"),
    ("specificity_pct", "Specificity"),
    ("accuracy_called_pct", "Called accuracy"),
    ("no_call_pct", "No-call rate")
]:

    ax.plot(
        summary.z_threshold,
        summary[column],
        marker="o",
        label=label
    )

ax.axvline(
    3.0,
    linestyle="--",
    label="Original threshold"
)

ax.set(
    xlabel="Z-score threshold",
    ylabel="Percentage (%)",
    title="Aneuploidy Screening Threshold Analysis",
    ylim=(0, 105)
)

ax.legend()
ax.grid(alpha=0.2)

fig.tight_layout()

fig.savefig(
    DEST / "zscore_threshold_analysis.png",
    dpi=300
)

plt.close(fig)

print("\nThreshold analysis completed")
print("Results saved to:", DEST)
