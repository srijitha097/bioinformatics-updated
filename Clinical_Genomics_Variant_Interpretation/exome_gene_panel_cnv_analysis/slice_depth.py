"""Summarize public HCC1143 exome slices and exploratory tumor/normal depth ratios."""

import argparse
import math
import statistics
import subprocess
from pathlib import Path

from coverage import summarize
from run import report, sha

from common import write

NORMAL = "SM-74NEG"
TUMOR = "SM-74P4M"


def analyze(data, out, samtools):
    out.mkdir(parents=True, exist_ok=True)
    version = subprocess.check_output([samtools, "--version"], text=True).splitlines()[
        0
    ]
    source = data / "ice_targets_sample-chr20.interval_list"
    spans = []
    for line in source.read_text().splitlines():
        if not line or line.startswith("@"):
            continue
        chrom, first, last = line.split()[:3]
        start, end = int(first) - 1, int(last)
        if chrom != "20" or not 0 <= start < end <= 1_000_000:
            raise ValueError(f"Unexpected slice interval {line}")
        if spans and start < spans[-1][1]:
            raise ValueError("Capture intervals overlap or are out of order")
        spans.append((start, end))
    bed = out / "ice_chr20.bed"
    bed.write_text("".join(f"20\t{start}\t{end}\n" for start, end in spans))
    input_hashes = {source.name: sha(source)}
    depth_by_sample = {}
    for sample in (NORMAL, TUMOR):
        bam = data / f"{sample}-v1-chr20-downsampled.deduplicated.bam"
        subprocess.run([samtools, "quickcheck", str(bam)], check=True)
        input_hashes[bam.name] = sha(bam)
        depth_file = data / f"{sample}.depth.tsv"
        with depth_file.open("w") as stream:
            subprocess.run(
                [
                    samtools,
                    "depth",
                    "-a",
                    "-b",
                    str(bed),
                    "-s",
                    "-q",
                    "20",
                    "-Q",
                    "20",
                    "-G",
                    "SUPPLEMENTARY",
                    str(bam),
                ],
                stdout=stream,
                check=True,
            )
        summarize(bed, depth_file, out / f"{sample}_coverage.json")
        depth_by_sample[sample] = {
            int(position) - 1: int(depth)
            for _, position, depth in (
                line.split() for line in depth_file.read_text().splitlines()
            )
        }
    records = []
    for start, end in spans:
        normal, tumor = [
            sum(depth_by_sample[sample].get(i, 0) for i in range(start, end))
            / (end - start)
            for sample in (NORMAL, TUMOR)
        ]
        log_ratio = math.log2(tumor / normal) if min(normal, tumor) >= 10 else None
        records.append(
            dict(
                chrom="20",
                start=start,
                end=end,
                normal_depth=normal,
                tumor_depth=tumor,
                raw_log2_ratio=log_ratio,
            )
        )
    retained = [r["raw_log2_ratio"] for r in records if r["raw_log2_ratio"] is not None]
    if not retained:
        raise ValueError("No targets pass the two-sample mean-depth floor")
    center = statistics.median(retained)
    for record in records:
        value = record["raw_log2_ratio"]
        record["centered_log2_ratio"] = value - center if value is not None else None
    write(out / "real_copy_ratios.tsv", records, list(records[0]))
    report(
        out / "real_copy_ratio_summary.json",
        dict(
            data_type="real public downsampled HCC1143 exome alignment slices",
            samtools=version,
            minimum_base_quality=20,
            minimum_mapping_quality=20,
            target_mean_depth_floor_both_samples=10,
            targets=len(records),
            retained_targets=len(retained),
            median_raw_log2_ratio=center,
            range_centered_log2_ratio=[min(retained) - center, max(retained) - center],
            input_sha256=input_hashes,
            limitations="Relative tumor/normal ratios only. Centering can remove regional shifts. "
            "No segmentation, purity/ploidy adjustment, or independent CNV truth for the slice.",
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args()
    analyze(args.data, args.out, args.samtools)
