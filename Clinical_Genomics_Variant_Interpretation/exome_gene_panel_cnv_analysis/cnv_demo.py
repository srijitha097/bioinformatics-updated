"""Benchmark a simple depth-based CNV screen on separately simulated controls and cases."""

import argparse
import math
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import report, write

TARGETS = 100
CONTROLS = 16
CASES = 12


def screen(depths, controls):
    normalized = [normalize(c) for c in controls]
    baseline = [statistics.median(column) for column in zip(*normalized)]
    ratios = [
        math.log2(depth / expected)
        for depth, expected in zip(normalize(depths), baseline)
    ]
    return [
        statistics.median(ratios[max(0, i - 2) : i + 3]) for i in range(len(ratios))
    ]


def normalize(depths):
    median = statistics.median(depths)
    return [depth / median for depth in depths]


def segments(states):
    start = 0
    for end in range(1, len(states) + 1):
        if end == len(states) or states[end] != states[start]:
            if states[start] != 2:
                yield start, end, states[start]
            start = end


def demo(out):
    rng = random.Random(42)
    efficiencies = [rng.uniform(0.5, 1.5) for _ in range(TARGETS)]

    def simulate(copy_numbers):
        library = rng.uniform(65, 140)
        return [
            max(1, round(library * ratio * cn / 2 * rng.lognormvariate(0, 0.20)))
            for ratio, cn in zip(efficiencies, copy_numbers)
        ]

    controls = [simulate([2] * TARGETS) for _ in range(CONTROLS)]
    records = []
    events = []
    for sample in range(CASES):
        truth = [2] * TARGETS
        if sample % 3:
            width = (2, 4, 7, 10)[sample // 3]
            start = rng.randrange(6, TARGETS - width - 5)
            truth[start : start + width] = [1 if sample % 3 == 1 else 3] * width
        depths = simulate(truth)
        ratios = screen(depths, controls)
        calls = [1 if ratio < -0.5 else 3 if ratio > 0.35 else 2 for ratio in ratios]
        for i, (depth, cn, call, ratio) in enumerate(zip(depths, truth, calls, ratios)):
            records.append(
                dict(
                    sample=sample,
                    target=i,
                    depth=depth,
                    true_cn=cn,
                    called_cn=call,
                    normalized_log2_ratio=round(ratio, 5),
                )
            )
        for channel, states in (("truth", truth), ("call", calls)):
            for start, end, cn in segments(states):
                events.append(
                    dict(
                        sample=sample,
                        source=channel,
                        start_target=start,
                        end_target=end,
                        direction="loss" if cn == 1 else "gain",
                    )
                )
    tp = sum(r["true_cn"] != 2 and r["called_cn"] == r["true_cn"] for r in records)
    fp = sum(r["called_cn"] != 2 and r["called_cn"] != r["true_cn"] for r in records)
    fn = sum(r["true_cn"] != 2 and r["called_cn"] != r["true_cn"] for r in records)
    write(out / "synthetic_targets.tsv", records, list(records[0]))
    write(out / "synthetic_segments.tsv", events, list(events[0]))
    report(
        out / "cnv_summary.json",
        dict(
            synthetic_only=True,
            seed=42,
            target_level=True,
            controls=CONTROLS,
            samples=CASES,
            targets_per_sample=TARGETS,
            tp=tp,
            fp=fp,
            fn=fn,
            precision=tp / (tp + fp) if tp + fp else None,
            recall=tp / (tp + fn) if tp + fn else None,
            limitation="Independent simulated copy-neutral control samples estimate per-target baselines. "
            "Shared capture efficiencies and independent log-normal noise simplify real assays. "
            "Five-target smoothing can miss short gains or losses. This is a target-level exercise.",
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    demo(parser.parse_args().out)
