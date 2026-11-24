import argparse
import gzip
import random
from pathlib import Path

import pandas as pd
import pysam

BASES = "ACGT"


def reverse_complement(seq):
    return seq.translate(
        str.maketrans("ACGTN", "TGCAN")
    )[::-1]


def introduce_errors(seq, rate, rng):
    return "".join(
        rng.choice([b for b in BASES if b != base])
        if base in BASES and rng.random() < rate
        else base
        for base in seq
    )


def write_fastq(handle, name, seq):
    quality = "?" * len(seq)
    handle.write(
        f"@{name}\n{seq}\n+\n{quality}\n"
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--molecules", type=int, default=500)

    args = parser.parse_args()

    reference_path = (
        "ngs_pipeline/data/reference/hbb_region.fa"
    )

    snp_path = (
        "ngs_pipeline/data/snps/informative_snps.tsv"
    )

    out = Path("ngs_pipeline/data/prenatal")
    out.mkdir(parents=True, exist_ok=True)

    snps = pd.read_csv(snp_path, sep="\t")

    with pysam.FastaFile(reference_path) as fasta:
        reference = fasta.fetch(
            fasta.references[0]
        ).upper()

    rng = random.Random(args.seed)

    manifest = []
    all_truth = []

    for ff in [0.00, 0.02, 0.04, 0.08, 0.12, 0.20]:

        sample = f"prenatal_ff{int(ff*100):02d}"

        total_molecules = 0
        fetal_molecules = 0
        alt_molecules = 0

        r1_path = out / f"{sample}_R1.fastq.gz"
        r2_path = out / f"{sample}_R2.fastq.gz"

        with (
            gzip.open(r1_path, "wt") as r1,
            gzip.open(r2_path, "wt") as r2
        ):

            for row in snps.itertuples():

                position = (
                    int(row.local_position_1based) - 1
                )

                ref = row.ref
                alt = row.alt

                if reference[position] != ref:
                    raise ValueError(
                        f"Reference mismatch: {row.rsid}"
                    )

                marker_alt = 0
                marker_fetal = 0

                for molecule in range(args.molecules):

                    is_fetal = rng.random() < ff

                    is_alt = (
                        is_fetal
                        and rng.random() < 0.5
                    )

                    if is_fetal:
                        marker_fetal += 1
                        fetal_molecules += 1

                    if is_alt:
                        marker_alt += 1
                        alt_molecules += 1

                    fragment_length = max(
                        140,
                        min(
                            190,
                            int(rng.gauss(166, 12))
                        )
                    )

                    start = position - rng.randint(20, 65)
                    end = start + fragment_length

                    fragment = list(
                        reference[start:end]
                    )

                    if is_alt:
                        fragment[position-start] = alt

                    fragment = "".join(fragment)

                    seq1 = introduce_errors(
                        fragment[:75], 0.001, rng
                    )

                    seq2 = introduce_errors(
                        reverse_complement(fragment[-75:]),
                        0.001,
                        rng
                    )

                    name = (
                        f"{sample}:{row.rsid}:{molecule}"
                    )

                    write_fastq(r1, name, seq1)
                    write_fastq(r2, name, seq2)

                    total_molecules += 1

                all_truth.append({
                    "sample": sample,
                    "rsid": row.rsid,
                    "target_ff": ff,
                    "molecules": args.molecules,
                    "fetal_molecules": marker_fetal,
                    "alternate_molecules": marker_alt
                })

        manifest.append({
            "sample": sample,
            "target_ff": ff,
            "realized_ff": (
                fetal_molecules / total_molecules
            ),
            "molecules": total_molecules,
            "fetal_molecules": fetal_molecules,
            "alternate_molecules": alt_molecules,
            "true_vaf": (
                alt_molecules / total_molecules
            ),
            "r1": str(r1_path),
            "r2": str(r2_path)
        })

        print(
            sample,
            "Molecules:", total_molecules,
            "Target FF:", round(ff, 4),
            "Realized FF:",
            round(fetal_molecules/total_molecules, 4),
            "True VAF:",
            round(alt_molecules/total_molecules, 4)
        )

    pd.DataFrame(manifest).to_csv(
        out / "sample_manifest.csv",
        index=False
    )

    pd.DataFrame(all_truth).to_csv(
        out / "marker_ground_truth.csv",
        index=False
    )

    print("\nPrenatal sequencing simulation completed")


if __name__ == "__main__":
    main()
