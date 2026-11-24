from collections import Counter
from pathlib import Path

import pandas as pd
import pysam

REFERENCE = Path(
    "ngs_pipeline/data/reference/aneuploidy/combined_reference.fa"
)

OUTPUT = REFERENCE.parent / "selected_windows.tsv"

WINDOW_SIZE = 1000
K = 75

BASES = set("ACGT")


def reverse_complement(seq):
    return seq.translate(
        str.maketrans("ACGT", "TGCA")
    )[::-1]


def canonical(seq):
    return min(seq, reverse_complement(seq))


with pysam.FastaFile(str(REFERENCE)) as fasta:

    sequences = {
        chrom: fasta.fetch(chrom).upper()
        for chrom in fasta.references
    }

print("Calculating combined-reference k-mer uniqueness")

kmers = Counter()

for chrom, seq in sequences.items():

    for i in range(len(seq) - K + 1):

        kmer = seq[i:i+K]

        if set(kmer) <= BASES:
            kmers[canonical(kmer)] += 1

candidates = []

for chrom, seq in sequences.items():

    for start in range(
        0,
        len(seq) - WINDOW_SIZE + 1,
        WINDOW_SIZE
    ):

        end = start + WINDOW_SIZE

        region = seq[start:end]

        if not set(region) <= BASES:
            continue

        gc = (
            region.count("G")
            + region.count("C")
        ) / WINDOW_SIZE

        if not 0.25 <= gc <= 0.65:
            continue

        local_kmers = [
            canonical(region[i:i+K])
            for i in range(len(region)-K+1)
        ]

        uniqueness = sum(
            kmers[k] == 1
            for k in local_kmers
        ) / len(local_kmers)

        if uniqueness < 0.99:
            continue

        candidates.append({
            "chromosome": chrom,
            "start_0based": start,
            "end_0based": end,
            "length": WINDOW_SIZE,
            "gc_content": gc,
            "local_uniqueness": uniqueness
        })

df = pd.DataFrame(candidates)

if df.empty:
    raise RuntimeError("No suitable windows found")

chromosomes = [
    "chr1", "chr2", "chr3",
    "chr13", "chr18", "chr21"
]

available = (
    df.groupby("chromosome")
    .size()
    .reindex(chromosomes, fill_value=0)
)

windows_per_chromosome = min(
    30,
    int(available.min())
)

if windows_per_chromosome < 10:
    raise RuntimeError(
        "Insufficient high-quality chromosome windows"
    )

selected = []

for chrom, group in df.groupby("chromosome"):

    group = group.sort_values(
        ["local_uniqueness", "start_0based"],
        ascending=[False, True]
    )

    selected.append(
        group.head(windows_per_chromosome)
    )

final = pd.concat(selected, ignore_index=True)

final.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)

summary = final.groupby("chromosome").agg(
    windows=("chromosome", "size"),
    selected_bases=("length", "sum"),
    mean_gc=("gc_content", "mean"),
    mean_uniqueness=("local_uniqueness", "mean")
)

summary.to_csv(
    REFERENCE.parent / "window_qc_summary.csv"
)

print("\nAvailable candidate windows:")
print(available.to_string())

print("\nSelected chromosome windows:")
print(summary.round(5).to_string())

print("\nTotal windows:", len(final))
print("Total selected bases:", final.length.sum())
print("Output:", OUTPUT)
