"""Summarize samtools depth -aa restricted to target intervals; missing positions are zero."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run import contains, contig, intervals

from common import op, report


def summarize(bed, depth, out):
    targets = intervals(bed)
    n = sum(e - s for starts, ends in targets.values() for s, e in zip(starts, ends))
    total = covered20 = covered30 = observed = 0
    seen = set()
    with op(depth) as f:
        for line in f:
            c, p, d = line.split()
            p, d = int(p), int(d)
            if d < 0:
                raise ValueError("Negative depth")
            if not contains(targets, c, p - 1):
                continue
            c = contig(c)
            if (c, p) in seen:
                raise ValueError("Duplicate position")
            seen.add((c, p))
            observed += 1
            total += d
            covered20 += d >= 20
            covered30 += d >= 30
    report(
        out,
        dict(
            target_bases=n,
            observed_positions=observed,
            mean_depth=total / n,
            fraction_20x=covered20 / n,
            fraction_30x=covered30 / n,
            missing_positions_treated_as_zero=n - observed,
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--bed", required=True)
    p.add_argument("--depth", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    summarize(a.bed, a.depth, a.out)
