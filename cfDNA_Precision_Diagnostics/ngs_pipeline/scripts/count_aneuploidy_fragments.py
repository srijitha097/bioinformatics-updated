import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd
import pysam

DATA = Path("ngs_pipeline/data/aneuploidy")

REFERENCE = Path(
    "ngs_pipeline/data/reference/aneuploidy"
)

ALIGNMENT = Path(
    "results/ngs/aneuploidy/alignment"
)

OUT = Path(
    "results/ngs/aneuploidy/benchmark"
)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--min-mapq", type=int, default=30
    )

    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)

    windows = pd.read_csv(
        REFERENCE / "selected_windows.tsv",
        sep="\t"
    )

    manifest = pd.read_csv(
        DATA / "sample_manifest.csv"
    )

    window_lookup = {
        (
            row.chromosome,
            int(row.start_0based)
        ): int(row.end_0based)
        for row in windows.itertuples()
    }

    results = []
    sample_qc = []

    for sample in manifest.itertuples():

        counts = Counter()
        filtered = Counter()

        bam_path = (
            ALIGNMENT / f"{sample.sample}.bam"
        )

        with pysam.AlignmentFile(
            bam_path, "rb"
        ) as bam:

            for read in bam.fetch(until_eof=True):

                if (
                    read.is_secondary
                    or read.is_supplementary
                    or not read.is_read1
                ):
                    continue

                if (
                    read.is_unmapped
                    or read.mate_is_unmapped
                    or read.is_qcfail
                ):
                    filtered["alignment_flag"] += 1
                    continue

                if not read.is_proper_pair:
                    filtered["improper_pair"] += 1
                    continue

                if read.mapping_quality < args.min_mapq:
                    filtered["low_mapq"] += 1
                    continue

                if (
                    read.reference_id
                    != read.next_reference_id
                ):
                    filtered["different_chromosome"] += 1
                    continue

                length = abs(
                    read.template_length
                )

                if not 140 <= length <= 190:
                    filtered["fragment_length"] += 1
                    continue

                start = min(
                    read.reference_start,
                    read.next_reference_start
                )

                end = start + length

                chrom = read.reference_name

                window_start = (
                    start // 1000
                ) * 1000

                key = (chrom, window_start)

                if key not in window_lookup:
                    filtered["outside_windows"] += 1
                    continue

                if end > window_lookup[key]:
                    filtered["cross_window_boundary"] += 1
                    continue

                counts[key] += 1

        recovered = sum(counts.values())

        for row in windows.itertuples():

            key = (
                row.chromosome,
                int(row.start_0based)
            )

            results.append({
                "sample": sample.sample,
                "chromosome": row.chromosome,
                "start_0based": row.start_0based,
                "end_0based": row.end_0based,
                "gc_content": row.gc_content,
                "fragment_count": counts[key]
            })

        sample_qc.append({
            "sample": sample.sample,
            "expected_fragments": sample.fragments,
            "recovered_fragments": recovered,
            "recovery_rate":
                recovered / sample.fragments,
            "filtered_fragments":
                sum(filtered.values()),
            "filter_details": json.dumps(
                dict(filtered)
            )
        })

        print(
            sample.sample,
            "| Recovered:", recovered,
            "| Recovery:",
            round(
                recovered / sample.fragments * 100,
                4
            ),
            "%",
            flush=True
        )

    df = pd.DataFrame(results)

    df.to_csv(
        OUT / "window_fragment_counts.csv",
        index=False
    )

    chrom = (
        df.groupby(
            ["sample", "chromosome"],
            as_index=False
        )["fragment_count"].sum()
    )

    chrom.to_csv(
        OUT / "chromosome_fragment_counts.csv",
        index=False
    )

    qc = pd.DataFrame(sample_qc)

    qc.to_csv(
        OUT / "fragment_recovery_qc.csv",
        index=False
    )

    # Ground truth is used only for evaluation.
    truth = pd.read_csv(
        DATA / "window_ground_truth.csv"
    )

    comparison = df.merge(
        truth[[
            "sample",
            "chromosome",
            "start_0based",
            "end_0based",
            "simulated_fragments"
        ]],
        on=[
            "sample",
            "chromosome",
            "start_0based",
            "end_0based"
        ],
        validate="one_to_one"
    )

    assert len(comparison) == len(df) == len(truth)

    comparison["count_error"] = (
        comparison["fragment_count"]
        - comparison["simulated_fragments"]
    )

    comparison.to_csv(
        OUT / "fragment_truth_comparison.csv",
        index=False
    )

    metrics = {
        "samples": int(len(qc)),
        "window_records": int(len(df)),
        "chromosome_records": int(len(chrom)),
        "expected_fragments": int(
            qc.expected_fragments.sum()
        ),
        "recovered_fragments": int(
            qc.recovered_fragments.sum()
        ),
        "mean_recovery_pct": float(
            qc.recovery_rate.mean() * 100
        ),
        "exact_window_count_matches": int(
            (comparison.count_error == 0).sum()
        ),
        "window_count_mae": float(
            comparison.count_error.abs().mean()
        )
    }

    with open(
        OUT / "fragment_count_metrics.json", "w"
    ) as handle:
        json.dump(metrics, handle, indent=2)

    print("\nFRAGMENT COUNT METRICS")

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
