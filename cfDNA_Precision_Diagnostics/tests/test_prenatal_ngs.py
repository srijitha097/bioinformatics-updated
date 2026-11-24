import numpy as np

from prenatal_genomic_screening.run import estimate_ff


def test_negative_control():
    depth = np.full(30, 500)
    alt = np.zeros(30)

    ff, lower, upper = estimate_ff(
        alt, depth, error=0.0003
    )

    assert ff < 0.001
    assert lower == 0
    assert upper > ff


def test_fetal_fraction_estimation():
    depth = np.full(30, 500)
    alt = np.full(30, 25)

    ff, lower, upper = estimate_ff(
        alt, depth, error=0.0003
    )

    assert abs(ff - 0.10) < 0.002
    assert lower < ff < upper


def test_confidence_interval_contains_truth():
    depth = np.full(30, 500)
    alt = np.full(30, 20)

    ff, lower, upper = estimate_ff(
        alt, depth, error=0.0003
    )

    assert lower <= 0.08 <= upper


def test_depth_improves_precision():
    depth_low = np.full(30, 100)
    alt_low = np.full(30, 5)

    depth_high = np.full(30, 1000)
    alt_high = np.full(30, 50)

    _, lo1, hi1 = estimate_ff(
        alt_low, depth_low, error=0.0003
    )

    _, lo2, hi2 = estimate_ff(
        alt_high, depth_high, error=0.0003
    )

    assert (hi2 - lo2) < (hi1 - lo1)
