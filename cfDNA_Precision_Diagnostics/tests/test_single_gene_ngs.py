import numpy as np
import pandas as pd

from collections import Counter

from ngs_pipeline.scripts.quantify_single_gene import (
    umi_neighbors,
    correct_umis,
    calculate_vaf
)

from ngs_pipeline.scripts.infer_single_gene_genotypes import (
    evaluate
)

from single_gene_screening.run import (
    prior,
    posterior,
    probabilities
)


def test_umi_neighbor_generation():
    umi = "AAAAAAAAAAAA"
    neighbors = set(umi_neighbors(umi))
    assert len(neighbors) == 36
    assert umi not in neighbors

    for neighbor in neighbors:
        assert sum(
            a != b for a, b in zip(umi, neighbor)
        ) == 1


def test_umi_error_correction():
    families = {
        "AAAAAAAAAAAA": Counter({"A": 5}),
        "AAAAAAAAAAAT": Counter({"A": 1})
    }

    corrected = correct_umis(families)
    assert len(corrected) == 1
    assert corrected["AAAAAAAAAAAA"]["A"] == 6


def test_equal_size_umis_not_merged():
    families = {
        "AAAAAAAAAAAA": Counter({"A": 3}),
        "AAAAAAAAAAAT": Counter({"T": 3})
    }

    assert len(correct_umis(families)) == 2


def test_ambiguous_umi_parent():
    families = {
        "AAAAAAAAAAAA": Counter({"A": 1}),
        "CAAAAAAAAAAA": Counter({"A": 5}),
        "GAAAAAAAAAAA": Counter({"G": 5})
    }

    assert len(correct_umis(families)) == 3


def test_vaf_calculation():
    counts = Counter({"A": 60, "T": 40, "C": 5})

    assert np.isclose(
        calculate_vaf(counts, "A", "T"), 0.4
    )


def test_empty_allele_counts():
    assert np.isnan(
        calculate_vaf(Counter(), "A", "T")
    )


def test_mendelian_carrier_prior():
    assert np.allclose(
        prior(1, 1),
        [0.25, 0.50, 0.25]
    )


def test_unknown_paternal_prior():
    assert np.allclose(
        prior(1, None, q=0.02),
        [0.49, 0.50, 0.01]
    )


def test_expected_allele_frequencies():
    assert np.allclose(
        probabilities(1, 0.20, error=0),
        [0.40, 0.50, 0.60]
    )


def test_posterior_normalization():
    result = posterior(
        6000, 10000,
        maternal=1,
        ff=0.20,
        paternal=1,
        ff_sd=0.005
    )

    assert np.all(np.isfinite(result))
    assert np.isclose(result.sum(), 1.0)
    assert np.all(result >= 0)


def test_affected_genotype_prediction():
    result = posterior(
        6000, 10000,
        maternal=1,
        ff=0.20,
        paternal=1,
        ff_sd=0.005
    )

    assert np.argmax(result) == 2
    assert result[2] >= 0.95


def test_no_call_aware_metrics():
    df = pd.DataFrame({
        "truth": [0, 1, 2, 2, 0],
        "call": [0, 1, 2, -1, 2]
    })

    metrics = evaluate(df)

    assert metrics["evaluations"] == 5
    assert metrics["called"] == 4
    assert metrics["correct"] == 3
    assert metrics["incorrect"] == 1
    assert metrics["no_calls"] == 1
    assert metrics["true_positives"] == 1
    assert metrics["false_positives"] == 1

    assert np.isclose(
        metrics["accuracy_called"], 0.75
    )

    assert np.isclose(
        metrics["accuracy_all"], 0.60
    )

    assert np.isclose(
        metrics["affected_detection_rate"], 0.50
    )
