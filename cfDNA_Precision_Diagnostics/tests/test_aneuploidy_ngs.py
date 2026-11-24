import numpy as np
import pandas as pd

from ngs_pipeline.scripts.analyze_aneuploidy import (
    gc_correct,
    classify,
    evaluate
)


def test_euploid_classification():
    assert classify(
        [], "PASS", 0.12, 0.04
    ) == "euploid"


def test_trisomy_classification():
    for chrom in [13, 18, 21]:
        assert classify(
            [f"chr{chrom}"],
            "PASS",
            0.12,
            0.04
        ) == f"trisomy{chrom}"


def test_low_fetal_fraction_no_call():
    assert classify(
        ["chr21"], "PASS", 0.02, 0.04
    ) == "NO_CALL_LOW_FF"


def test_failed_qc_no_call():
    assert classify(
        [], "FAIL", 0.12, 0.04
    ) == "NO_CALL_QC"


def test_multiple_chromosome_no_call():
    assert classify(
        ["chr13", "chr21"],
        "PASS",
        0.12,
        0.04
    ) == "NO_CALL_MULTIPLE"


def test_gc_correction_reduces_variability():
    gc = np.linspace(0.30, 0.60, 30)

    counts = 500 * np.exp(
        -10 * (gc - 0.45)**2
    )

    df = pd.DataFrame({
        "chromosome": ["chr1"] * 30,
        "gc_content": gc,
        "fragment_count": counts
    })

    corrected = gc_correct(df)

    raw_cv = (
        df.fragment_count.std()
        / df.fragment_count.mean()
    )

    corrected_cv = (
        corrected.corrected_count.std()
        / corrected.corrected_count.mean()
    )

    assert corrected_cv < raw_cv


def test_screening_metrics():
    df = pd.DataFrame({
        "truth": [
            "euploid",
            "trisomy13",
            "trisomy21",
            "trisomy18"
        ],
        "call": [
            "euploid",
            "trisomy13",
            "euploid",
            "NO_CALL_LOW_FF"
        ]
    })

    metrics = evaluate(df, "call")

    assert metrics["called_samples"] == 3
    assert metrics["no_calls"] == 1
    assert metrics["correct_calls"] == 2
    assert metrics["specificity_called"] == 1.0
    assert metrics[
        "exact_trisomy_sensitivity_called"
    ] == 0.5

    assert np.isclose(
        metrics["accuracy_called"], 2/3
    )
