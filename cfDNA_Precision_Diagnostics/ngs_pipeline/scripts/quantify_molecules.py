import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pysam


def hamming_distance(a, b):
    if len(a) != len(b):
        return float("inf")
    return sum(x != y for x, y in zip(a, b))


def correct_umis(families):
    counts = {
        umi: len(bases)
        for umi, bases in families.items()
    }

    ordered = sorted(
        counts,
        key=lambda u: (-counts[u], u)
    )

    parent = {umi: umi for umi in ordered}

    for major in ordered:
        if parent[major] != major:
            continue

        for minor in ordered:
            if minor == major:
                continue

            if parent[minor] != minor:
                continue

            if counts[major] <= counts[minor]:
                continue

            if counts[major] < 2 * counts[minor] - 1:
                continue

            if hamming_distance(major, minor) == 1:
                parent[minor] = major

    corrected = defaultdict(list)

    for umi, bases in families.items():
        corrected[parent[umi]].extend(bases)

    return corrected


def consensus_call(bases, min_size, threshold):
    if len(bases) < min_size:
        return None

    counts = Counter(bases)
    ranked = counts.most_common()

    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        return None

    base, count = ranked[0]

    if count / len(bases) < threshold:
        return None

    return base


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--sample", default="sample_02")
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--min-family-size", type=int, default=2)
    parser.add_argument("--consensus", type=float, default=0.75)

    args = parser.parse_args()

    sample = args.sample

    bam_path = Path(
        f"results/ngs/alignment/{sample}.sorted.bam"
    )

    variant_path = Path(
        f"ngs_pipeline/data/annotations/{sample}_variant.tsv"
    )

    truth_path = Path(
        f"ngs_pipeline/data/annotations/{sample}_truth.csv"
    )

    output = Path("results/ngs/molecular")
    output.mkdir(parents=True, exist_ok=True)

    with open(variant_path) as handle:
        variant = next(
            csv.DictReader(handle, delimiter="\t")
        )

    with open(truth_path) as handle:
        truth = list(csv.DictReader(handle))

    contig = variant["contig"]
    position = int(variant["local_position_1based"]) - 1

    ref = variant["ref"]
    alt = variant["alt"]
    true_vaf = float(variant["true_vaf"])

    observations = defaultdict(list)
    filtered = Counter()

    with pysam.AlignmentFile(bam_path, "rb") as bam:

        for read in bam.fetch(
            contig, position, position + 1
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

            if read.mapping_quality < args.min_mapq:
                filtered["low_mapq"] += 1
                continue

            if not 140 <= abs(read.template_length) <= 190:
                filtered["fragment_length"] += 1
                continue

            match = re.search(
                r"UMI_([ACGT]{12})$",
                read.query_name
            )

            if not match:
                filtered["missing_umi"] += 1
                continue

            umi = match.group(1)

            for qpos, rpos in read.get_aligned_pairs():

                if rpos != position or qpos is None:
                    continue

                if read.query_qualities[qpos] < args.min_baseq:
                    filtered["low_baseq"] += 1
                    continue

                start = min(
                    read.reference_start,
                    read.next_reference_start
                )

                end = start + abs(read.template_length)

                base = read.query_sequence[qpos].upper()

                observations[read.query_name].append(
                    (start, end, umi, base)
                )

    families = defaultdict(lambda: defaultdict(list))
    raw_counts = Counter()

    for name, records in observations.items():

        if len(set(records)) != 1:
            filtered["discordant_mates"] += 1
            continue

        start, end, umi, base = records[0]

        raw_counts[base] += 1
        families[(start, end)][umi].append(base)

    uncorrected_count = sum(
        len(group) for group in families.values()
    )

    consensus_counts = Counter()
    family_sizes = Counter()

    corrected_count = 0
    rejected_count = 0

    truth_lookup = {
        (
            int(row["start"]),
            int(row["end"]),
            row["umi"]
        ): row["allele"]
        for row in truth
    }

    recovered_truth = set()
    correct_truth_calls = 0

    with open(
        output / f"{sample}_families.csv", "w"
    ) as handle:

        writer = csv.writer(handle)

        writer.writerow([
            "start",
            "end",
            "umi",
            "family_size",
            "consensus",
            "truth_allele"
        ])

        for (start, end), group in families.items():

            corrected = correct_umis(group)
            corrected_count += len(corrected)

            for umi, bases in corrected.items():

                call = consensus_call(
                    bases,
                    args.min_family_size,
                    args.consensus
                )

                family_sizes[len(bases)] += 1

                key = (start, end, umi)
                expected = truth_lookup.get(key)

                if call is None:
                    rejected_count += 1
                    call = "NO_CALL"

                else:
                    consensus_counts[call] += 1

                    if expected is not None:
                        recovered_truth.add(key)

                        if call == expected:
                            correct_truth_calls += 1

                writer.writerow([
                    start,
                    end,
                    umi,
                    len(bases),
                    call,
                    expected or "UNMATCHED"
                ])

    raw_total = raw_counts[ref] + raw_counts[alt]

    consensus_total = (
        consensus_counts[ref] +
        consensus_counts[alt]
    )

    raw_vaf = (
        raw_counts[alt] / raw_total
        if raw_total else None
    )

    consensus_vaf = (
        consensus_counts[alt] / consensus_total
        if consensus_total else None
    )

    summary = {
        "sample": sample,
        "input_molecules": len(truth),
        "reference_allele": ref,
        "alternate_allele": alt,
        "true_vaf": true_vaf,
        "raw_ref_count": raw_counts[ref],
        "raw_alt_count": raw_counts[alt],
        "raw_other_count": sum(raw_counts.values()) - raw_total,
        "raw_vaf": raw_vaf,
        "uncorrected_umi_families": uncorrected_count,
        "corrected_umi_families": corrected_count,
        "rejected_families": rejected_count,
        "consensus_ref_molecules": consensus_counts[ref],
        "consensus_alt_molecules": consensus_counts[alt],
        "consensus_other_molecules": (
            sum(consensus_counts.values()) - consensus_total
        ),
        "consensus_vaf": consensus_vaf,
        "raw_absolute_error": (
            abs(raw_vaf - true_vaf)
            if raw_vaf is not None else None
        ),
        "consensus_absolute_error": (
            abs(consensus_vaf - true_vaf)
            if consensus_vaf is not None else None
        ),
        "exact_truth_matched_molecules": len(recovered_truth),
        "exact_truth_recovery_rate": (
            len(recovered_truth) / len(truth)
            if truth else None
        ),
        "matched_consensus_accuracy": (
            correct_truth_calls / len(recovered_truth)
            if recovered_truth else None
        ),
        "filter_counts": dict(filtered),
        "family_size_distribution": dict(family_sizes)
    }

    with open(
        output / f"{sample}_summary.json", "w"
    ) as handle:
        json.dump(summary, handle, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
