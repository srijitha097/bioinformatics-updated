import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

REFERENCE = Path(
    "ngs_pipeline/data/reference/single_gene/combined_reference.fa"
)

TARGETS = Path(
    "ngs_pipeline/data/clinvar/selected_clinvar_targets.tsv"
)

OUT = Path("ngs_pipeline/data/single_gene")

BASES = "ACGT"

TRANSLATION = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq):
    return seq.translate(TRANSLATION)[::-1]


def mutate(seq, rate, rng):

    n_errors = rng.binomial(len(seq), rate)

    if n_errors == 0:
        return seq

    seq = list(seq)

    positions = rng.choice(
        len(seq),
        size=n_errors,
        replace=False
    )

    for i in positions:

        options = [
            b for b in BASES
            if b != seq[i]
        ]

        seq[i] = options[rng.integers(0, 3)]

    return "".join(seq)


def generate_umi(rng, used):

    while True:

        value = int(
            rng.integers(0, 4**12)
        )

        if value not in used:

            used.add(value)

            break

    bases = []

    for _ in range(12):

        bases.append(BASES[value % 4])
        value //= 4

    return "".join(bases)


def write_fastq(handle, name, seq):

    handle.write(
        f"@{name}\n{seq}\n+\n"
        + "?" * len(seq)
        + "\n"
    )


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--seed", type=int, default=42
    )

    parser.add_argument(
        "--replicates", type=int, default=3
    )

    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)

    targets = pd.read_csv(
        TARGETS, sep="\t"
    )

    with pysam.FastaFile(str(REFERENCE)) as fasta:

        sequences = {
            gene: fasta.fetch(gene).upper()
            for gene in targets.gene
        }

    rng = np.random.default_rng(args.seed)

    manifest = []
    truth = []

    for ff in [0.02, 0.04, 0.08, 0.12, 0.20]:

        for depth in [500, 2000, 10000]:

            for genotype in [0, 1, 2]:

                for rep in range(
                    1, args.replicates + 1
                ):

                    sample = (
                        f"sg_ff{round(ff*100):02d}"
                        f"_n{depth}_g{genotype}_r{rep}"
                    )

                    r1_path = (
                        OUT / f"{sample}_R1.fastq.gz"
                    )

                    r2_path = (
                        OUT / f"{sample}_R2.fastq.gz"
                    )

                    pair_id = 0
                    total_molecules = 0

                    # Simulated independent FF measurement.
                    estimated_ff = float(
                        np.clip(
                            rng.normal(ff, 0.005),
                            0.001,
                            0.5
                        )
                    )

                    with (
                        gzip.open(
                            r1_path, "wt",
                            compresslevel=3
                        ) as r1,
                        gzip.open(
                            r2_path, "wt",
                            compresslevel=3
                        ) as r2
                    ):

                        for row in targets.itertuples():

                            gene = row.gene

                            seq = sequences[gene]

                            position = int(
                                row.local_position_1based
                            ) - 1

                            if seq[position] != row.ref:
                                raise ValueError(
                                    f"REF mismatch: {gene}"
                                )

                            # Heterozygous maternal genotype.
                            expected_vaf = (
                                (1 - ff) * 0.5
                                + ff * genotype / 2
                            )

                            # Match the overdispersion
                            # assumption of the existing
                            # Bayesian model.
                            rho = 0.0005
                            k = 1 / rho - 1

                            sample_vaf = rng.beta(
                                expected_vaf * k,
                                (1 - expected_vaf) * k
                            )

                            alt_flags = (
                                rng.random(depth)
                                < sample_vaf
                            )

                            used_umis = set()

                            alt_count = int(
                                alt_flags.sum()
                            )

                            gene_pairs = 0

                            for molecule in range(depth):

                                umi = generate_umi(
                                    rng, used_umis
                                )

                                length = int(
                                    np.clip(
                                        rng.normal(166, 12),
                                        140,
                                        190
                                    )
                                )

                                start = (
                                    position
                                    - int(
                                        rng.integers(20, 66)
                                    )
                                )

                                end = start + length

                                fragment = list(
                                    seq[start:end]
                                )

                                if alt_flags[molecule]:

                                    fragment[
                                        position - start
                                    ] = row.alt

                                fragment = "".join(fragment)

                                copies = min(
                                    6,
                                    1 + rng.poisson(2)
                                )

                                for _ in range(copies):

                                    observed_umi = mutate(
                                        umi, 0.002, rng
                                    )

                                    seq1 = mutate(
                                        fragment[:75],
                                        0.001,
                                        rng
                                    )

                                    seq2 = mutate(
                                        reverse_complement(
                                            fragment[-75:]
                                        ),
                                        0.001,
                                        rng
                                    )

                                    name = (
                                        f"{sample}:{pair_id}:"
                                        f"UMI_{observed_umi}"
                                    )

                                    write_fastq(
                                        r1, name, seq1
                                    )

                                    write_fastq(
                                        r2, name, seq2
                                    )

                                    pair_id += 1
                                    gene_pairs += 1

                            total_molecules += depth

                            truth.append({
                                "sample": sample,
                                "gene": gene,
                                "clinvar_id": row.clinvar_id,
                                "true_ff": ff,
                                "estimated_ff": estimated_ff,
                                "maternal_genotype": 1,
                                "fetal_genotype": genotype,
                                "molecules": depth,
                                "alternate_molecules": alt_count,
                                "expected_vaf": expected_vaf,
                                "simulated_vaf": sample_vaf,
                                "realized_vaf": alt_count / depth,
                                "read_pairs": gene_pairs
                            })

                    manifest.append({
                        "sample": sample,
                        "true_ff": ff,
                        "estimated_ff": estimated_ff,
                        "fetal_genotype": genotype,
                        "molecules_per_gene": depth,
                        "total_molecules": total_molecules,
                        "replicate": rep,
                        "read_pairs": pair_id,
                        "reads": pair_id * 2,
                        "r1": str(r1_path),
                        "r2": str(r2_path)
                    })

                    print(
                        sample,
                        "| Molecules:", total_molecules,
                        "| Read pairs:", pair_id,
                        flush=True
                    )

    pd.DataFrame(manifest).to_csv(
        OUT / "sample_manifest.csv",
        index=False
    )

    pd.DataFrame(truth).to_csv(
        OUT / "molecular_ground_truth.csv",
        index=False
    )

    print("\nSIMULATION COMPLETED")
    print("Sequencing samples:", len(manifest))
    print("Sample-gene experiments:", len(truth))

    print(
        "Original molecules:",
        sum(
            row["total_molecules"]
            for row in manifest
        )
    )

    print(
        "Sequencing reads:",
        sum(
            row["reads"]
            for row in manifest
        )
    )


if __name__ == "__main__":
    main()
