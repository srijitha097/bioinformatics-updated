import argparse
import csv
import gzip
import random
from pathlib import Path

import pysam

BASES = "ACGT"


def reverse_complement(sequence):
    return sequence.translate(
        str.maketrans("ACGTN", "TGCAN")
    )[::-1]


def introduce_errors(sequence, rate, rng):
    result = []

    for base in sequence:
        if base in BASES and rng.random() < rate:
            base = rng.choice(
                [b for b in BASES if b != base]
            )
        result.append(base)

    return "".join(result)


def write_fastq(handle, name, sequence):
    quality = "?" * len(sequence)
    handle.write(
        f"@{name}\n{sequence}\n+\n{quality}\n"
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--reference",
        default="ngs_pipeline/data/reference/hbb_region.fa"
    )

    parser.add_argument("--molecules", type=int, default=5000)
    parser.add_argument("--vaf", type=float, default=0.02)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sample", default="sample_01")

    args = parser.parse_args()

    rng = random.Random(args.seed)

    output = Path("ngs_pipeline/data/fastq")
    output.mkdir(parents=True, exist_ok=True)

    metadata = Path("ngs_pipeline/data/annotations")
    metadata.mkdir(parents=True, exist_ok=True)

    with pysam.FastaFile(args.reference) as fasta:
        chromosome = fasta.references[0]
        reference = fasta.fetch(chromosome).upper()

    position = int(
        Path(
            "ngs_pipeline/data/annotations/selected_target.txt"
        ).read_text().strip()
    )

    ref = reference[position]
    alt = rng.choice([b for b in BASES if b != ref])

    read_length = 75
    total_reads = 0
    alternate_molecules = 0
    used_umis = set()

    r1_path = output / f"{args.sample}_R1.fastq.gz"
    r2_path = output / f"{args.sample}_R2.fastq.gz"

    truth_path = metadata / f"{args.sample}_truth.csv"

    with (
        gzip.open(r1_path, "wt") as r1,
        gzip.open(r2_path, "wt") as r2,
        open(truth_path, "w") as truth
    ):
        writer = csv.writer(truth)

        writer.writerow([
            "molecule_id",
            "umi",
            "start",
            "end",
            "allele",
            "pcr_copies"
        ])

        for molecule in range(args.molecules):

            while True:
                umi = "".join(
                    rng.choices(BASES, k=12)
                )
                if umi not in used_umis:
                    used_umis.add(umi)
                    break

            fragment_length = max(
                140,
                min(190, int(rng.gauss(166, 12)))
            )

            start = position - rng.randint(20, 65)
            end = start + fragment_length

            fragment = list(reference[start:end])

            is_alt = rng.random() < args.vaf

            if is_alt:
                fragment[position - start] = alt
                alternate_molecules += 1

            fragment = "".join(fragment)

            copies = max(
                1,
                min(15, int(rng.expovariate(1/4)) + 1)
            )

            writer.writerow([
                molecule,
                umi,
                start,
                end,
                alt if is_alt else ref,
                copies
            ])

            for copy in range(copies):

                observed_umi = introduce_errors(
                    umi, 0.002, rng
                )

                seq1 = introduce_errors(
                    fragment[:read_length],
                    0.005,
                    rng
                )

                seq2 = introduce_errors(
                    reverse_complement(
                        fragment[-read_length:]
                    ),
                    0.005,
                    rng
                )

                name = (
                    f"{args.sample}:"
                    f"{total_reads}:"
                    f"UMI_{observed_umi}"
                )

                write_fastq(r1, name, seq1)
                write_fastq(r2, name, seq2)

                total_reads += 1

    variant_path = metadata / f"{args.sample}_variant.tsv"

    with open(variant_path, "w") as handle:
        handle.write(
            "contig\tlocal_position_1based\t"
            "genomic_position_1based\tref\talt\t"
            "true_vaf\n"
        )

        handle.write(
            f"{chromosome}\t{position+1}\t"
            f"{5200001+position}\t{ref}\t{alt}\t"
            f"{alternate_molecules/args.molecules:.6f}\n"
        )

    print("Sample:", args.sample)
    print("Original molecules:", args.molecules)
    print("Alternate molecules:", alternate_molecules)
    print("True VAF:", alternate_molecules/args.molecules)
    print("Paired-end read pairs:", total_reads)
    print("Reference contig:", chromosome)
    print("Local variant position:", position + 1)
    print("Reference allele:", ref)
    print("Alternate allele:", alt)
    print("FASTQ R1:", r1_path)
    print("FASTQ R2:", r2_path)
    print("Ground truth:", truth_path)


if __name__ == "__main__":
    main()
