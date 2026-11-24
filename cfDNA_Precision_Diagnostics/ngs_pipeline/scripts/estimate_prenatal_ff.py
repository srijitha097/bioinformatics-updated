import argparse
import sys
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from prenatal_genomic_screening.run import estimate_ff


def count_marker(bam, marker, min_mapq=30, min_baseq=20):

    position = int(marker.local_position_1based) - 1

    observations = defaultdict(set)
    filtered = Counter()

    for read in bam.fetch(
        marker.contig,
        position,
        position + 1
    ):

        if (
            read.is_unmapped
            or read.is_secondary
            or read.is_supplementary
            or read.is_qcfail
        ):
            filtered["alignment_flag"] += 1
            continue

        if not read.is_proper_pair:
            filtered["improper_pair"] += 1
            continue

        if read.mapping_quality < min_mapq:
            filtered["low_mapq"] += 1
            continue

        for qpos, rpos in read.get_aligned_pairs():

            if rpos != position or qpos is None:
                continue

            if read.query_qualities[qpos] < min_baseq:
                filtered["low_baseq"] += 1
                continue

            base = read.query_sequence[qpos].upper()

            observations[read.query_name].add(base)

    counts = Counter()

    for name, bases in observations.items():

        if len(bases) != 1:
            filtered["discordant_mates"] += 1
            continue

        counts[next(iter(bases))] += 1

    ref_count = counts[marker.ref]
    alt_count = counts[marker.alt]

    depth = ref_count + alt_count

    return {
        "rsid": marker.rsid,
        "chromosome": marker.chromosome,
        "genomic_position_1based":
            marker.genomic_position_1based,
        "local_position_1based":
            marker.local_position_1based,
        "ref": marker.ref,
        "alt": marker.alt,
        "ref_count": ref_count,
        "alt_count": alt_count,
        "other_count": sum(counts.values()) - depth,
        "depth": depth,
        "observed_vaf":
            alt_count / depth if depth else np.nan,
        "filtered_observations": sum(filtered.values())
    }


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--min-mapq", type=int, default=30
    )

    parser.add_argument(
        "--min-baseq", type=int, default=20
    )

    parser.add_argument(
        "--min-marker-depth", type=int, default=100
    )

    args = parser.parse_args()

    snps = pd.read_csv(
        ROOT / "ngs_pipeline/data/snps/informative_snps.tsv",
        sep="\t"
    )

    manifest = pd.read_csv(
        ROOT / "ngs_pipeline/data/prenatal/sample_manifest.csv"
    )

    alignment = ROOT / "results/ngs/prenatal/alignment"

    out = ROOT / "results/ngs/fetal_fraction"
    out.mkdir(parents=True, exist_ok=True)

    all_counts = []

    for sample in manifest.itertuples():

        bam_path = alignment / f"{sample.sample}.bam"

        print("Counting SNP alleles:", sample.sample)

        with pysam.AlignmentFile(bam_path, "rb") as bam:

            for marker in snps.itertuples():

                result = count_marker(
                    bam,
                    marker,
                    args.min_mapq,
                    args.min_baseq
                )

                result["sample"] = sample.sample

                all_counts.append(result)

    counts = pd.DataFrame(all_counts)

    counts["passed_depth_qc"] = (
        counts["depth"] >= args.min_marker_depth
    )

    counts.to_csv(
        out / "snp_allele_counts.csv",
        index=False
    )

    control = counts[
        (counts["sample"] == "prenatal_ff00")
        & counts["passed_depth_qc"]
    ]

    if control["depth"].sum() == 0:
        raise RuntimeError(
            "Negative control has no usable depth"
        )

    # Estimate effective alternate allele error
    # using the negative control.
    error = (
        control["alt_count"].sum() + 0.5
    ) / (
        control["depth"].sum() + 1
    )

    print(
        "Calibrated alternate allele error:",
        round(error, 7)
    )

    summary = []

    for sample in manifest.itertuples():

        group = counts[
            (counts["sample"] == sample.sample)
            & counts["passed_depth_qc"]
        ]

        if group.empty:
            raise RuntimeError(
                f"No usable SNP markers: {sample.sample}"
            )

        alt = group["alt_count"].to_numpy()
        depth = group["depth"].to_numpy()

        ff, lower, upper = estimate_ff(
            alt,
            depth,
            error=error
        )

        observed_vaf = alt.sum() / depth.sum()

        qc_status = (
            "PASS"
            if len(group) >= 20
            and depth.sum() >= 10000
            else "FAIL"
        )

        if qc_status == "FAIL":
            screening_status = "NO_CALL_QC"

        elif sample.sample == "prenatal_ff00":
            screening_status = "NEGATIVE_CONTROL"

        elif lower < 0.04:
            screening_status = "NO_CALL_LOW_FF"

        else:
            screening_status = "FF_THRESHOLD_PASS"

        summary.append({
            "sample": sample.sample,
            "target_ff": sample.target_ff,
            "realized_ff": sample.realized_ff,
            "estimated_ff": ff,
            "ff_lower_95": lower,
            "ff_upper_95": upper,
            "absolute_error":
                abs(ff - sample.realized_ff),
            "markers_passed": len(group),
            "total_allele_depth": int(depth.sum()),
            "mean_marker_depth": depth.mean(),
            "reference_count": int(
                group["ref_count"].sum()
            ),
            "alternate_count": int(alt.sum()),
            "other_count": int(
                group["other_count"].sum()
            ),
            "observed_vaf": observed_vaf,
            "true_vaf": sample.true_vaf,
            "qc_status": qc_status,
            "screening_status": screening_status,
            "ci_contains_target_ff":
                lower <= sample.target_ff <= upper
        })

    results = pd.DataFrame(summary)

    results.to_csv(
        out / "fetal_fraction_summary.csv",
        index=False
    )

    evaluation = results[
        results["sample"] != "prenatal_ff00"
    ]

    metrics = {
        "samples": int(len(results)),
        "evaluation_samples": int(len(evaluation)),
        "annotated_markers": int(len(snps)),
        "calibrated_error_rate": float(error),
        "total_allele_observations": int(
            results["total_allele_depth"].sum()
        ),
        "mean_absolute_error_percentage_points": float(
            evaluation["absolute_error"].mean() * 100
        ),
        "rmse_percentage_points": float(
            np.sqrt(
                np.mean(evaluation["absolute_error"] ** 2)
            ) * 100
        ),
        "target_ff_ci_coverage": float(
            evaluation["ci_contains_target_ff"].mean()
        ),
        "qc_passed_samples": int(
            (results["qc_status"] == "PASS").sum()
        ),
        "ff_threshold_passed": int(
            (
                results["screening_status"]
                == "FF_THRESHOLD_PASS"
            ).sum()
        )
    }

    with open(
        out / "benchmark_metrics.json", "w"
    ) as handle:

        json.dump(metrics, handle, indent=2)

    x = results["realized_ff"].to_numpy() * 100
    y = results["estimated_ff"].to_numpy() * 100

    lower = results["ff_lower_95"].to_numpy() * 100
    upper = results["ff_upper_95"].to_numpy() * 100

    yerr = np.vstack([
        np.maximum(y - lower, 0),
        np.maximum(upper - y, 0)
    ])

    fig, ax = plt.subplots(figsize=(7, 6))

    ax.errorbar(
        x,
        y,
        yerr=yerr,
        fmt="o",
        capsize=5,
        label="BAM-derived fetal fraction"
    )

    ax.plot(
        [0, 22],
        [0, 22],
        linestyle="--",
        label="Perfect agreement"
    )

    ax.set(
        xlabel="Realized fetal fraction (%)",
        ylabel="Estimated fetal fraction (%)",
        title="Reference-Based cfDNA Fetal Fraction Estimation"
    )

    ax.legend()
    ax.grid(alpha=0.2)

    fig.tight_layout()

    fig.savefig(
        out / "fetal_fraction_benchmark.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print("\nFetal fraction estimates:")

    print(
        results[[
            "sample",
            "realized_ff",
            "estimated_ff",
            "ff_lower_95",
            "ff_upper_95",
            "markers_passed",
            "total_allele_depth",
            "screening_status"
        ]].round(5).to_string(index=False)
    )

    print("\nBenchmark metrics:")
    print(json.dumps(metrics, indent=2))

    print("\nResults saved to:", out)


if __name__ == "__main__":
    main()
