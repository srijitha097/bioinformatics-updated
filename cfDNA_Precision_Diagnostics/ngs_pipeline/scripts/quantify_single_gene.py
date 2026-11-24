import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

DATA = Path("ngs_pipeline/data/single_gene")
OUT = Path("results/ngs/single_gene/benchmark")
ALIGNMENT = Path("results/ngs/single_gene/alignment")

TARGETS = Path(
    "ngs_pipeline/data/clinvar/selected_clinvar_targets.tsv"
)

BASES = "ACGT"


def umi_neighbors(umi):

    for i, base in enumerate(umi):

        for replacement in BASES:

            if replacement != base:

                yield (
                    umi[:i]
                    + replacement
                    + umi[i+1:]
                )


def correct_umis(families):

    sizes = {
        umi: sum(counts.values())
        for umi, counts in families.items()
    }

    parents = {}

    for umi, size in sizes.items():

        candidates = [
            neighbor
            for neighbor in umi_neighbors(umi)
            if neighbor in sizes
            and sizes[neighbor] > size
            and sizes[neighbor] >= 2 * size - 1
        ]

        if candidates:

            best_size = max(
                sizes[x] for x in candidates
            )

            best = [
                x for x in candidates
                if sizes[x] == best_size
            ]

            if len(best) == 1:
                parents[umi] = best[0]

    corrected = defaultdict(Counter)

    for umi, counts in families.items():

        root = umi

        while root in parents:
            root = parents[root]

        corrected[root].update(counts)

    return corrected


def get_target_base(read, position, min_baseq):

    for query_pos, ref_pos in read.get_aligned_pairs(
        matches_only=True
    ):

        if ref_pos == position:

            if read.query_qualities is None:
                return None

            if (
                read.query_qualities[query_pos]
                < min_baseq
            ):
                return None

            return read.query_sequence[
                query_pos
            ].upper()

    return None


def calculate_vaf(counts, ref, alt):

    depth = counts[ref] + counts[alt]

    return (
        counts[alt] / depth
        if depth else np.nan
    )


def quantify_gene(
    bam,
    target,
    min_mapq,
    min_baseq,
    min_family,
    consensus_threshold
):

    position = int(
        target.local_position_1based
    ) - 1

    families = defaultdict(
        lambda: defaultdict(Counter)
    )

    raw_counts = Counter()
    filters = Counter()

    for read in bam.fetch(
        target.gene,
        position,
        position + 1
    ):

        # The simulator places the target in R1.
        # Count only R1 to avoid double counting mates.
        if (
            not read.is_read1
            or read.is_secondary
            or read.is_supplementary
        ):
            continue

        if (
            read.is_unmapped
            or read.mate_is_unmapped
            or not read.is_proper_pair
            or read.is_qcfail
        ):
            filters["alignment_flag"] += 1
            continue

        if read.mapping_quality < min_mapq:
            filters["low_mapq"] += 1
            continue

        if (
            read.reference_id
            != read.next_reference_id
        ):
            filters["different_contig"] += 1
            continue

        length = abs(read.template_length)

        if not 140 <= length <= 190:
            filters["fragment_length"] += 1
            continue

        if ":UMI_" not in read.query_name:
            filters["missing_umi"] += 1
            continue

        umi = read.query_name.rsplit(
            ":UMI_", 1
        )[-1]

        if (
            len(umi) != 12
            or not set(umi) <= set(BASES)
        ):
            filters["invalid_umi"] += 1
            continue

        base = get_target_base(
            read,
            position,
            min_baseq
        )

        if base is None or base not in BASES:
            filters["invalid_base"] += 1
            continue

        start = min(
            read.reference_start,
            read.next_reference_start
        )

        end = start + length

        key = (start, end)

        families[key][umi][base] += 1

        raw_counts[base] += 1

    uncorrected_families = sum(
        len(group)
        for group in families.values()
    )

    consensus_counts = Counter()
    corrected_families = 0
    rejected_small = 0
    rejected_ambiguous = 0
    accepted_sizes = []

    for group in families.values():

        corrected = correct_umis(group)

        corrected_families += len(corrected)

        for counts in corrected.values():

            size = sum(counts.values())

            if size < min_family:
                rejected_small += 1
                continue

            ranking = counts.most_common()

            base, support = ranking[0]

            if (
                len(ranking) > 1
                and ranking[1][1] == support
            ):
                rejected_ambiguous += 1
                continue

            if support / size < consensus_threshold:

                rejected_ambiguous += 1
                continue

            consensus_counts[base] += 1
            accepted_sizes.append(size)

    ref = target.ref
    alt = target.alt

    raw_depth = (
        raw_counts[ref]
        + raw_counts[alt]
    )

    consensus_depth = (
        consensus_counts[ref]
        + consensus_counts[alt]
    )

    return {
        "raw_ref": raw_counts[ref],
        "raw_alt": raw_counts[alt],
        "raw_other": sum(
            n for b, n in raw_counts.items()
            if b not in {ref, alt}
        ),
        "raw_depth": raw_depth,
        "raw_vaf": calculate_vaf(
            raw_counts, ref, alt
        ),
        "uncorrected_umi_families":
            uncorrected_families,
        "corrected_umi_families":
            corrected_families,
        "merged_umi_families":
            uncorrected_families
            - corrected_families,
        "rejected_small_families":
            rejected_small,
        "rejected_ambiguous_families":
            rejected_ambiguous,
        "consensus_ref":
            consensus_counts[ref],
        "consensus_alt":
            consensus_counts[alt],
        "consensus_other": sum(
            n for b, n in consensus_counts.items()
            if b not in {ref, alt}
        ),
        "consensus_depth": consensus_depth,
        "consensus_vaf": calculate_vaf(
            consensus_counts, ref, alt
        ),
        "mean_accepted_family_size": (
            float(np.mean(accepted_sizes))
            if accepted_sizes else 0
        ),
        "filtered_reads": sum(filters.values()),
        "filter_details": json.dumps(dict(filters))
    }


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--min-mapq",
        type=int,
        default=30
    )

    parser.add_argument(
        "--min-baseq",
        type=int,
        default=20
    )

    parser.add_argument(
        "--min-family",
        type=int,
        default=2
    )

    parser.add_argument(
        "--consensus",
        type=float,
        default=0.75
    )

    args = parser.parse_args()

    if args.min_family < 1:
        raise ValueError("Invalid minimum family size")

    if not 0.5 < args.consensus <= 1:
        raise ValueError("Invalid consensus threshold")

    OUT.mkdir(parents=True, exist_ok=True)

    targets = pd.read_csv(
        TARGETS, sep="\t"
    )

    manifest = pd.read_csv(
        DATA / "sample_manifest.csv"
    )

    rows = []

    for sample in manifest.itertuples():

        path = ALIGNMENT / f"{sample.sample}.bam"

        with pysam.AlignmentFile(
            str(path), "rb"
        ) as bam:

            for target in targets.itertuples():

                result = quantify_gene(
                    bam,
                    target,
                    args.min_mapq,
                    args.min_baseq,
                    args.min_family,
                    args.consensus
                )

                result.update({
                    "sample": sample.sample,
                    "gene": target.gene,
                    "clinvar_id": target.clinvar_id
                })

                rows.append(result)

                print(
                    sample.sample,
                    target.gene,
                    "| Raw depth:",
                    result["raw_depth"],
                    "| Consensus depth:",
                    result["consensus_depth"],
                    flush=True
                )

    df = pd.DataFrame(rows)

    df.to_csv(
        OUT / "molecular_allele_counts.csv",
        index=False
    )

    # Simulation truth is used only for benchmarking.
    truth = pd.read_csv(
        DATA / "molecular_ground_truth.csv"
    )

    comparison = df.merge(
        truth,
        on=["sample", "gene", "clinvar_id"],
        validate="one_to_one"
    )

    assert len(comparison) == len(df) == len(truth)

    comparison["raw_vaf_error"] = (
        comparison.raw_vaf
        - comparison.realized_vaf
    ).abs()

    comparison["consensus_vaf_error"] = (
        comparison.consensus_vaf
        - comparison.realized_vaf
    ).abs()

    comparison["family_count_ratio"] = (
        comparison.consensus_depth
        / comparison.molecules
    )

    comparison.to_csv(
        OUT / "molecular_truth_comparison.csv",
        index=False
    )

    raw_mae = float(
        comparison.raw_vaf_error.mean()
    )

    consensus_mae = float(
        comparison.consensus_vaf_error.mean()
    )

    metrics = {
        "sequencing_samples":
            int(df["sample"].nunique()),
        "gene_experiments":
            int(len(df)),
        "original_molecules":
            int(comparison.molecules.sum()),
        "raw_allele_observations":
            int(df.raw_depth.sum()),
        "uncorrected_umi_families":
            int(df.uncorrected_umi_families.sum()),
        "corrected_umi_families":
            int(df.corrected_umi_families.sum()),
        "merged_umi_families":
            int(df.merged_umi_families.sum()),
        "rejected_small_families":
            int(df.rejected_small_families.sum()),
        "rejected_ambiguous_families":
            int(df.rejected_ambiguous_families.sum()),
        "accepted_ref_alt_families":
            int(df.consensus_depth.sum()),
        "mean_family_count_ratio_pct":
            float(
                comparison.family_count_ratio.mean()
                * 100
            ),
        "raw_vaf_mae": raw_mae,
        "consensus_vaf_mae": consensus_mae,
        "vaf_error_reduction_pct": (
            float(
                (1 - consensus_mae / raw_mae)
                * 100
            )
            if raw_mae > 0 else None
        ),
        "min_family_size": args.min_family,
        "consensus_threshold": args.consensus
    }

    with open(
        OUT / "molecular_quantification_metrics.json",
        "w"
    ) as handle:

        json.dump(metrics, handle, indent=2)

    summary = comparison.groupby(
        ["gene", "molecules"]
    ).agg(
        experiments=("sample", "size"),
        mean_consensus_depth=("consensus_depth", "mean"),
        mean_family_count_ratio=("family_count_ratio", "mean"),
        raw_vaf_mae=("raw_vaf_error", "mean"),
        consensus_vaf_mae=("consensus_vaf_error", "mean")
    ).reset_index()

    summary.to_csv(
        OUT / "molecular_performance_by_depth.csv",
        index=False
    )

    print("\nMOLECULAR QUANTIFICATION METRICS")
    print(json.dumps(metrics, indent=2))

    print("\nGENE AND DEPTH SUMMARY")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
