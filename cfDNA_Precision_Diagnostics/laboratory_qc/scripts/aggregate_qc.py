import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

RESULTS = ROOT / "results/ngs"
OUT = ROOT / "results/laboratory_qc"

CONFIG = ROOT / "laboratory_qc/config/thresholds.json"

OUT.mkdir(parents=True, exist_ok=True)

with open(CONFIG) as handle:
    thresholds = json.load(handle)

sources = {
    "Aneuploidy": (
        RESULTS / "aneuploidy/benchmark/"
        "alignment_qc_summary.csv"
    ),
    "Single-gene": (
        RESULTS / "single_gene/benchmark/"
        "alignment_qc_summary.csv"
    )
}

frames = []

for assay, path in sources.items():

    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(path)

    df["assay"] = assay

    frames.append(df)

qc = pd.concat(
    frames,
    ignore_index=True
)

alerts = []


def add_alert(row, metric, observed, threshold):

    alerts.append({
        "sample": row.sample,
        "assay": row.assay,
        "metric": metric,
        "observed_value": observed,
        "threshold": threshold,
        "severity": "FAIL"
    })


for row in qc.itertuples():

    checks = [
        (
            "mapping_pct",
            row.mapping_pct,
            thresholds["minimum_mapping_pct"]
        ),
        (
            "properly_paired_pct",
            row.properly_paired_pct,
            thresholds["minimum_pairing_pct"]
        ),
        (
            "high_mapq_pct",
            row.high_mapq_pct,
            thresholds["minimum_high_mapq_pct"]
        ),
        (
            "mean_mapq",
            row.mean_mapq,
            thresholds["minimum_mean_mapq"]
        ),
        (
            "mean_fragment_length",
            row.mean_fragment_length,
            thresholds["minimum_mean_fragment_length"]
        )
    ]

    for metric, observed, threshold in checks:

        if pd.isna(observed) or observed < threshold:

            add_alert(
                row, metric, observed, threshold
            )

    if (
        row.mean_fragment_length >
        thresholds["maximum_mean_fragment_length"]
    ):

        add_alert(
            row,
            "mean_fragment_length",
            row.mean_fragment_length,
            thresholds["maximum_mean_fragment_length"]
        )

    if row.expected_reads != row.total_primary_reads:

        add_alert(
            row,
            "primary_read_count",
            row.total_primary_reads,
            row.expected_reads
        )

    if row.qc_status != "PASS":

        add_alert(
            row,
            "upstream_qc_status",
            row.qc_status,
            "PASS"
        )

alert_df = pd.DataFrame(
    alerts,
    columns=[
        "sample",
        "assay",
        "metric",
        "observed_value",
        "threshold",
        "severity"
    ]
)

failed = set(
    zip(alert_df.assay, alert_df["sample"])
)

qc["monitoring_status"] = [
    "FAIL" if (a, s) in failed else "PASS"
    for a, s in zip(qc.assay, qc["sample"])
]

qc.to_csv(
    OUT / "centralized_alignment_qc.csv",
    index=False
)

alert_df.to_csv(
    OUT / "qc_alerts.csv",
    index=False
)

summary = qc.groupby("assay").agg(
    samples=("sample", "size"),
    total_reads=("total_primary_reads", "sum"),
    mean_mapping_pct=("mapping_pct", "mean"),
    mean_pairing_pct=("properly_paired_pct", "mean"),
    mean_mapq=("mean_mapq", "mean"),
    mean_fragment_length=("mean_fragment_length", "mean"),
    qc_passed=(
        "monitoring_status",
        lambda x: (x == "PASS").sum()
    ),
    qc_failed=(
        "monitoring_status",
        lambda x: (x == "FAIL").sum()
    )
).reset_index()

summary.to_csv(
    OUT / "assay_qc_summary.csv",
    index=False
)

metrics = {
    "assays": int(qc.assay.nunique()),
    "samples": int(len(qc)),
    "total_reads": int(
        qc.total_primary_reads.sum()
    ),
    "qc_passed": int(
        (qc.monitoring_status == "PASS").sum()
    ),
    "qc_failed": int(
        (qc.monitoring_status == "FAIL").sum()
    ),
    "total_alerts": int(len(alert_df))
}

with open(
    OUT / "laboratory_qc_metrics.json", "w"
) as handle:

    json.dump(metrics, handle, indent=2)

print("\nASSAY QC SUMMARY")
print(summary.to_string(index=False))

print("\nLABORATORY QC METRICS")
print(json.dumps(metrics, indent=2))

print("\nResults saved to:", OUT)
