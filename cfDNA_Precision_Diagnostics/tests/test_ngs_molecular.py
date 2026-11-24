from ngs_pipeline.scripts.quantify_molecules import (
    hamming_distance,
    correct_umis,
    consensus_call
)


def test_hamming_distance():
    assert hamming_distance("AAAA", "AAAT") == 1
    assert hamming_distance("AAAA", "AAAA") == 0


def test_directional_umi_correction():
    families = {
        "AAAAAAAAAAAA": ["T"] * 5,
        "AAAAAAAAAAAT": ["T"]
    }

    corrected = correct_umis(families)

    assert len(corrected) == 1
    assert len(corrected["AAAAAAAAAAAA"]) == 6


def test_unrelated_umis_preserved():
    families = {
        "AAAAAAAAAAAA": ["T"] * 5,
        "CCCCCCCCCCCC": ["G"] * 3
    }

    assert len(correct_umis(families)) == 2


def test_consensus_error_correction():
    assert consensus_call(
        ["T", "T", "T", "G"], 2, 0.75
    ) == "T"


def test_singleton_rejected():
    assert consensus_call(
        ["T"], 2, 0.75
    ) is None


def test_ambiguous_consensus_rejected():
    assert consensus_call(
        ["T", "G"], 2, 0.75
    ) is None
