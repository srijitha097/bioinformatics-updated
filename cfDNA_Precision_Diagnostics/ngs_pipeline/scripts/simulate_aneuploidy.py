import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

REFERENCE = Path(
    "ngs_pipeline/data/reference/aneuploidy/combined_reference.fa"
)

WINDOWS = REFERENCE.parent / "selected_windows.tsv"

OUTPUT = Path("ngs_pipeline/data/aneuploidy")

BASES = "ACGT"

CHROMOSOMES = [
    "chr1", "chr2", "chr3",
    "chr13", "chr18", "chr21"
]

TRANSLATION = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq):
    return seq.translate(TRANSLATION)[::-1]


def introduce_errors(seq, rng, rate=0.001):

    n_errors = rng.binomial(len(seq), rate)

    if n_errors == 0:
        return seq

    bases = list(seq)

    positions = rng.choice(
        len(seq),
        size=n_errors,
        replace=False
    )

    for position in positions:

        original = bases[position]

        alternatives = [
            b for b in BASES
            if b != original
        ]

        bases[position] = alternatives[
            rng.integers(0, 3)
        ]

    return "".join(bases)


def write_fastq(handle, name, sequence):

    quality = "?" * len(sequence)

    handle.write(
        f"@{name}\n{sequence}\n+\n{quality}\n"
    )


def create_design():

    design = []

    for rep in range(1, 21):

        design.append({
            "sample": f"calibration_{rep:02d}",
            "cohort": "calibration",
            "condition": "euploid",
            "trisomy": "none",
            "true_ff": 0.0,
            "replicate": rep
        })

    for condition in [
        "euploid",
        "trisomy13",
        "trisomy18",
        "trisomy21"
    ]:

        trisomy = (
            "none"
            if condition == "euploid"
            else "chr" + condition.replace("trisomy", "")
        )

        for ff in [0.02, 0.04, 0.08, 0.12, 0.20]:

            for rep in range(1, 4):

                design.append({
                    "sample": (
                        f"{condition}_ff"
                        f"{round(ff*100):02d}_r{rep}"
                    ),
                    "cohort": "evaluation",
                    "condition": condition,
                    "trisomy": trisomy,
                    "true_ff": ff,
                    "replicate": rep
                })

    return design


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--seed", type=int, default=42
    )

    parser.add_argument(
        "--fragments", type=int, default=60000
    )

    args = parser.parse_args()

    OUTPUT.mkdir(parents=True, exist_ok=True)

    windows = pd.read_csv(
        WINDOWS, sep="\t"
    )

    with pysam.FastaFile(str(REFERENCE)) as fasta:

        sequences = {
            chrom: fasta.fetch(chrom).upper()
            for chrom in CHROMOSOMES
        }

    gc = windows["gc_content"].to_numpy()

    chromosomes = windows[
        "chromosome"
    ].to_numpy()

    design = create_design()

    seeds = np.random.SeedSequence(
        args.seed
    ).spawn(len(design))

    manifest = []
    chromosome_truth = []
    window_truth = []

    for sample, seed in zip(design, seeds):

        rng = np.random.default_rng(seed)

        name = sample["sample"]

        ff = sample["true_ff"]

        # Sample-specific GC-dependent bias.
        gc_strength = rng.uniform(6, 12)

        gc_center = rng.normal(
            0.45, 0.015
        )

        gc_bias = np.exp(
            -gc_strength * (gc - gc_center)**2
        )

        # Small random technical variation between windows.
        technical_noise = rng.lognormal(
            mean=0,
            sigma=0.02,
            size=len(windows)
        )

        dosage = np.ones(len(windows))

        if sample["trisomy"] != "none":

            dosage[
                chromosomes == sample["trisomy"]
            ] = 1 + ff / 2

        weights = (
            gc_bias
            * dosage
            * technical_noise
        )

        probabilities = weights / weights.sum()

        fragment_counts = rng.multinomial(
            args.fragments,
            probabilities
        )

        r1_path = OUTPUT / f"{name}_R1.fastq.gz"
        r2_path = OUTPUT / f"{name}_R2.fastq.gz"

        fragment_id = 0

        with (
            gzip.open(
                r1_path, "wt", compresslevel=3
            ) as r1,
            gzip.open(
                r2_path, "wt", compresslevel=3
            ) as r2
        ):

            for i, row in enumerate(
                windows.itertuples()
            ):

                count = int(fragment_counts[i])

                sequence = sequences[row.chromosome]

                start_window = int(
                    row.start_0based
                )

                end_window = int(
                    row.end_0based
                )

                lengths = np.clip(
                    rng.normal(166, 12, count),
                    140,
                    190
                ).astype(int)

                for length in lengths:

                    start = int(
                        rng.integers(
                            start_window,
                            end_window - length + 1
                        )
                    )

                    fragment = sequence[
                        start:start+length
                    ]

                    seq1 = introduce_errors(
                        fragment[:75], rng
                    )

                    seq2 = introduce_errors(
                        reverse_complement(
                            fragment[-75:]
                        ),
                        rng
                    )

                    read_name = (
                        f"{name}:{fragment_id}"
                    )

                    write_fastq(
                        r1, read_name, seq1
                    )

                    write_fastq(
                        r2, read_name, seq2
                    )

                    fragment_id += 1

                window_truth.append({
                    "sample": name,
                    "chromosome": row.chromosome,
                    "start_0based": start_window,
                    "end_0based": end_window,
                    "gc_content": row.gc_content,
                    "simulated_fragments": count,
                    "expected_relative_dosage":
                        float(dosage[i])
                })

        for chrom in CHROMOSOMES:

            chrom_count = int(
                fragment_counts[
                    chromosomes == chrom
                ].sum()
            )

            chromosome_truth.append({
                "sample": name,
                "chromosome": chrom,
                "simulated_fragments": chrom_count,
                "fragment_fraction":
                    chrom_count / args.fragments
            })

        manifest.append({
            **sample,
            "fragments": args.fragments,
            "reads": args.fragments * 2,
            "gc_strength": gc_strength,
            "gc_center": gc_center,
            "r1": str(r1_path),
            "r2": str(r2_path)
        })

        print(
            name,
            "| Condition:", sample["condition"],
            "| FF:", ff,
            "| Fragments:", fragment_id,
            flush=True
        )

    pd.DataFrame(manifest).to_csv(
        OUTPUT / "sample_manifest.csv",
        index=False
    )

    pd.DataFrame(chromosome_truth).to_csv(
        OUTPUT / "chromosome_ground_truth.csv",
        index=False
    )

    pd.DataFrame(window_truth).to_csv(
        OUTPUT / "window_ground_truth.csv",
        index=False
    )

    print("\nSimulation completed")
    print("Samples:", len(manifest))

    print(
        "Total fragments:",
        len(manifest) * args.fragments
    )

    print(
        "Total sequencing reads:",
        len(manifest) * args.fragments * 2
    )


if __name__ == "__main__":
    main()
