"""Phase-two algorithm and saved-output integration checks."""

import csv
import hashlib
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "comparative_genomics"))
from check_benchmark import mask
from reviews import ONCOLOGY


def table(name):
    with (ROOT / "comparative_genomics" / name).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


class PhaseTwo(unittest.TestCase):
    def test_interval_endpoints(self):
        with TemporaryDirectory() as tmp:
            a, b = Path(tmp) / "a", Path(tmp) / "b"
            a.write_text("chr1\t0\t5\nchr1\t1\t3\nchr1\t5\t7\n")
            b.write_text("1\t4\t6\n")
            self.assertEqual(mask([a, b])["1"], [4, 6])

    def test_fixed_panel_frequency_null(self):
        cases = table("rare/cases.tsv")
        n = len(cases)
        self.assertEqual(len({r["gene"] for r in cases}), n)
        for r in table("rare/recovery.tsv"):
            if r["method"] == "frequency":
                for k in (1, 5, 10):
                    self.assertAlmostEqual(float(r[f"top{k}"]), k / n)
                self.assertAlmostEqual(
                    float(r["rr"]), sum(1 / i for i in range(1, n + 1)) / n
                )

    def test_tied_ranks(self):
        for r in table("rare/ranks.tsv"):
            l, h = int(r["best"]), int(r["worst"])
            self.assertAlmostEqual(
                float(r["rr"]), sum(1 / i for i in range(l, h + 1)) / (h - l + 1)
            )
            for k in (1, 5, 10):
                self.assertAlmostEqual(
                    float(r[f"top{k}"]), max(0, min(k, h) - l + 1) / (h - l + 1)
                )

    def test_reference_join(self):
        initial = table("oncology/assessments.tsv")
        joined = table("oncology/comparison.tsv")
        self.assertEqual(len(ONCOLOGY), len(initial))
        self.assertEqual(
            initial, [{k: v for k, v in r.items() if k != "reference"} for r in joined]
        )
        self.assertTrue(all(r["pmid"].isdigit() and r["publication"] for r in initial))

    def test_normalized_counts(self):
        cross = {
            (r["scope"], r["kind"]): r
            for r in table("sequencing/normalized_crosscheck.tsv")
        }
        for r in table("sequencing/metrics.tsv"):
            if r["DP"] == r["GQ"] == "0":
                for k in ("tp", "fp", "fn"):
                    self.assertEqual(r[k], cross[r["scope"], r["kind"]][k])
            self.assertEqual(int(r["tp"]) + int(r["fn"]), int(r["truth"]))
            self.assertEqual(int(r["tp"]) + int(r["fp"]), int(r["query"]))

    def test_cnv_lock(self):
        dev = table("sequencing/cnv_development.tsv")
        best = max(dev, key=lambda r: float(r["f1"]))
        selected = table("sequencing/cnv_holdout.tsv")[1]
        for k in (
            "normalize",
            "matched_controls",
            "gc_correction",
            "minimum_targets",
            "gain_threshold",
        ):
            self.assertEqual(best[k], selected[k])
        for r in dev + [selected]:
            tp, fp, fn = [int(r[k]) for k in ("tp", "fp", "fn")]
            self.assertAlmostEqual(float(r["f1"]), 2 * tp / (2 * tp + fp + fn))

    def test_baseline_files(self):
        for line in (
            (ROOT / "comparative_genomics/baseline_hashes.txt").read_text().splitlines()
        ):
            digest, path = line.split("  ", 1)
            self.assertEqual(
                hashlib.sha256((ROOT / path).read_bytes()).hexdigest(), digest
            )
