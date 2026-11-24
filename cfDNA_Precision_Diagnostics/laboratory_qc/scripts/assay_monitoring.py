import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

NGS = ROOT / "results/ngs"
DATA = ROOT / "ngs_pipeline/data"
OUT = ROOT / "results/laboratory_qc"

OUT.mkdir(parents=True, exist_ok=True)

with open(
    ROOT / "laboratory_qc/config/thresholds.json"
) as handle:

    config = json.load(handle)

columns = [
    "sample",
    "assay",
    "gene",
    "mode",
    "metric",
    "observed_value",
    "threshold",
    "severity",
    "category"
]

alerts = []


def alert(
    sample,
    assay,
    metric,
    observed,
    threshold,
    severity,
    category,
    gene="",
    mode=""
):

    alerts.append({
        "sample": sample,
        "assay": assay,
        "gene": gene,
        "mode": mode,
        "metric": metric,
        "observed_value": observed,
        "threshold": threshold,
        "severity": severity,
        "category": category
    })


# -----------------------------------
# Molecular quantification monitoring
# -----------------------------------

molecular = pd.read_csv(
    NGS / "single_gene/benchmark/"
    "molecular_truth_comparison.csv"
)

for row in molecular.itertuples():

    if (
        pd.isna(row.family_count_ratio)
        or row.family_count_ratio
        < config["minimum_family_count_ratio"]
    ):

        alert(
            row.sample,
            "Single-gene",
            "family_count_ratio",
            row.family_count_ratio,
            config["minimum_family_count_ratio"],
            "WARN",
            "BENCHMARK_QC",
            gene=row.gene
        )

    if (
        pd.isna(row.consensus_vaf_error)
        or row.consensus_vaf_error
        > config["maximum_vaf_error"]
    ):

        alert(
            row.sample,
            "Single-gene",
            "consensus_vaf_error",
            row.consensus_vaf_error,
            config["maximum_vaf_error"],
            "WARN",
            "BENCHMARK_QC",
            gene=row.gene
        )

    if row.consensus_depth < 100:

        alert(
            row.sample,
            "Single-gene",
            "consensus_depth",
            row.consensus_depth,
            100,
            "FAIL",
            "TECHNICAL_QC",
            gene=row.gene
        )


# -----------------------------------
# Fetal fraction monitoring
# -----------------------------------

manifest = pd.read_csv(
    DATA / "single_gene/sample_manifest.csv"
)

for row in manifest.itertuples():

    if (
        pd.isna(row.estimated_ff)
        or not 0 < row.estimated_ff <= 0.5
    ):

        alert(
            row.sample,
            "Single-gene",
            "estimated_fetal_fraction",
            row.estimated_ff,
            "0 < FF <= 0.5",
            "FAIL",
            "TECHNICAL_QC"
        )

    elif (
        row.estimated_ff
        < config["minimum_fetal_fraction"]
    ):

        alert(
            row.sample,
            "Single-gene",
            "estimated_fetal_fraction",
            row.estimated_ff,
            config["minimum_fetal_fraction"],
            "WARN",
            "LOW_FETAL_FRACTION"
        )


# -----------------------------------
# Bayesian posterior monitoring
# -----------------------------------

predictions = pd.read_csv(
    NGS / "single_gene/benchmark/"
    "genotype_predictions.csv"
)

for row in predictions.itertuples():

    if row.call_reason == "LOW_POSTERIOR":

        confidence = max(
            row.p0,
            row.p1,
            row.p2
        )

        alert(
            row.sample,
            "Single-gene",
            "posterior_confidence",
            confidence,
            0.95,
            "INFO",
            "NO_CALL",
            gene=row.gene,
            mode=row.mode
        )

    elif row.call_reason != "CALLED":

        alert(
            row.sample,
            "Single-gene",
            "genotype_call_reason",
            row.call_reason,
            "CALLED",
            "FAIL",
            "TECHNICAL_QC",
            gene=row.gene,
            mode=row.mode
        )


# -----------------------------------
# Aneuploidy no-call monitoring
# -----------------------------------

aneuploidy = pd.read_csv(
    NGS / "aneuploidy/benchmark/"
    "screening_calls.csv"
)

for row in aneuploidy.itertuples():

    if str(row.corrected_call).startswith("NO_CALL"):

        alert(
            row.sample,
            "Aneuploidy",
            "screening_call",
            row.corrected_call,
            "Reportable result",
            (
                "FAIL"
                if row.corrected_call == "NO_CALL_QC"
                else "INFO"
            ),
            (
                "TECHNICAL_QC"
                if row.corrected_call == "NO_CALL_QC"
                else "NO_CALL"
            )
        )


# -----------------------------------
# Generate monitoring reports
# -----------------------------------

report = pd.DataFrame(
    alerts,
    columns=columns
)

report.to_csv(
    OUT / "assay_monitoring_alerts.csv",
    index=False
)

summary = report.groupby(
    ["assay", "category", "severity"]
).size().reset_index(
    name="events"
)

summary.to_csv(
    OUT / "assay_monitoring_summary.csv",
    index=False
)

metrics = {
    "molecular_experiments": int(len(molecular)),
    "single_gene_samples": int(len(manifest)),
    "genotype_evaluations": int(len(predictions)),
    "aneuploidy_evaluations": int(len(aneuploidy)),
    "monitoring_events": int(len(report)),
    "technical_failures": int(
        (report.severity == "FAIL").sum()
    ),
    "warnings": int(
        (report.severity == "WARN").sum()
    ),
    "informational_events": int(
        (report.severity == "INFO").sum()
    ),
    "single_gene_no_calls": int(
        (predictions.call == -1).sum()
    ),
    "aneuploidy_no_calls": int(
        aneuploidy.corrected_call
        .str.startswith("NO_CALL").sum()
    )
}

with open(
    OUT / "assay_monitoring_metrics.json",
    "w"
) as handle:

    json.dump(metrics, handle, indent=2)

print("\nASSAY MONITORING SUMMARY")
print(summary.to_string(index=False))

print("\nMONITORING METRICS")
print(json.dumps(metrics, indent=2))

print("\nAssay monitoring completed")
