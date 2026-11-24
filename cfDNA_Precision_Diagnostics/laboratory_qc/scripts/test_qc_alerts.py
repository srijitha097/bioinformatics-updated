import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

DEST = ROOT / "results/laboratory_qc/validation"

DEST.mkdir(parents=True, exist_ok=True)

SCENARIOS = [
    ("LOW_MAPPING", "mapping_pct", 85.0),
    ("LOW_PAIRING", "properly_paired_pct", 80.0),
    ("LOW_HIGH_MAPQ", "high_mapq_pct", 75.0),
    ("LOW_MEAN_MAPQ", "mean_mapq", 20.0),
    ("SHORT_FRAGMENTS", "mean_fragment_length", 110.0),
    ("LONG_FRAGMENTS", "mean_fragment_length", 220.0),
    ("READ_COUNT_MISMATCH", "total_primary_reads", 100),
    ("UPSTREAM_FAILURE", "qc_status", "FAIL"),
    ("MISSING_MAPPING", "mapping_pct", float("nan"))
]

validation = []

with tempfile.TemporaryDirectory(
    prefix="cfdna_qc_test_"
) as temp:

    sandbox = Path(temp)

    script = (
        sandbox / "laboratory_qc/scripts/aggregate_qc.py"
    )

    script.parent.mkdir(parents=True)

    shutil.copy2(
        ROOT / "laboratory_qc/scripts/aggregate_qc.py",
        script
    )

    config = (
        sandbox / "laboratory_qc/config"
    )

    config.mkdir(parents=True)

    shutil.copy2(
        ROOT / "laboratory_qc/config/thresholds.json",
        config / "thresholds.json"
    )

    for assay in ["aneuploidy", "single_gene"]:

        relative = Path(
            f"results/ngs/{assay}/benchmark/"
            "alignment_qc_summary.csv"
        )

        target = sandbox / relative

        target.parent.mkdir(
            parents=True,
            exist_ok=True
        )

        shutil.copy2(
            ROOT / relative,
            target
        )

    source = (
        sandbox / "results/ngs/aneuploidy/benchmark/"
        "alignment_qc_summary.csv"
    )

    df = pd.read_csv(source)

    injected = {}

    for i, (scenario, metric, value) in enumerate(
        SCENARIOS
    ):

        sample = f"QC_TEST_{scenario}"

        df.loc[i, "sample"] = sample

        df.loc[i, metric] = value

        expected_metric = {
            "total_primary_reads": "primary_read_count",
            "qc_status": "upstream_qc_status"
        }.get(metric, metric)

        injected[sample] = {
            "scenario": scenario,
            "expected_metric": expected_metric
        }

    df.to_csv(source, index=False)

    subprocess.run(
        [sys.executable, str(script)],
        cwd=sandbox,
        check=True
    )

    output = sandbox / "results/laboratory_qc"

    alerts = pd.read_csv(
        output / "qc_alerts.csv"
    )

    samples = pd.read_csv(
        output / "centralized_alignment_qc.csv"
    )

    for sample, scenario in injected.items():

        matching = alerts[
            (alerts["sample"] == sample)
            & (
                alerts.metric
                == scenario["expected_metric"]
            )
        ]

        status = samples.loc[
            samples["sample"] == sample,
            "monitoring_status"
        ].iloc[0]

        detected = (
            not matching.empty
            and status == "FAIL"
        )

        validation.append({
            "scenario": scenario["scenario"],
            "sample": sample,
            "expected_metric":
                scenario["expected_metric"],
            "detected": bool(detected),
            "monitoring_status": status,
            "test_status":
                "PASS" if detected else "FAIL"
        })

    normal = samples[
        ~samples["sample"].isin(injected)
    ]

    false_alert_samples = int(
        (normal.monitoring_status == "FAIL").sum()
    )

    alerts.to_csv(
        DEST / "simulated_qc_alerts.csv",
        index=False
    )

    samples.to_csv(
        DEST / "simulated_sample_qc.csv",
        index=False
    )

report = pd.DataFrame(validation)

report.to_csv(
    DEST / "failure_detection_validation.csv",
    index=False
)

detected = int(report.detected.sum())

metrics = {
    "failure_scenarios": len(SCENARIOS),
    "detected_failures": detected,
    "missed_failures": len(SCENARIOS) - detected,
    "failure_detection_rate":
        detected / len(SCENARIOS),
    "normal_samples": len(normal),
    "false_alert_samples": false_alert_samples,
    "original_results_modified": False
}

with open(
    DEST / "validation_metrics.json", "w"
) as handle:

    json.dump(metrics, handle, indent=2)

print("\nQC FAILURE DETECTION VALIDATION")
print(report.to_string(index=False))

print("\nVALIDATION METRICS")
print(json.dumps(metrics, indent=2))

assert detected == len(SCENARIOS)
assert false_alert_samples == 0

print("\nALL QC FAILURE SCENARIOS PASSED")
