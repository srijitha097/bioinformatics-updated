from collections import Counter
from pathlib import Path

import pandas as pd
import pysam

OUT = Path(
    "ngs_pipeline/data/reference/single_gene"
)

targets = pd.read_csv(
    "ngs_pipeline/data/clinvar/selected_clinvar_targets.tsv",
    sep="\t"
)

sequences = {}
rows = []

K = 75

BASES = set("ACGT")

TRANSLATION = str.maketrans("ACGT", "TGCA")


def canonical(seq):

    reverse = seq.translate(TRANSLATION)[::-1]

    return min(seq, reverse)


for row in targets.itertuples():

    path = OUT / f"{row.gene}.fa"

    with pysam.FastaFile(str(path)) as fasta:

        seq = fasta.fetch(
            fasta.references[0]
        ).upper()

    if len(seq) != 10001:

        raise ValueError(
            f"Unexpected sequence length: {row.gene}"
        )

    position = int(row.local_position_1based) - 1

    actual_ref = seq[position]

    if actual_ref != row.ref:

        raise ValueError(
            f"{row.gene}: ClinVar REF={row.ref}, "
            f"GRCh38 REF={actual_ref}"
        )

    sequences[row.gene] = seq

    print(
        row.gene,
        "| Length:", len(seq),
        "| Reference allele:", actual_ref,
        "| ClinVar allele:", row.ref,
        "| MATCH"
    )

kmers = Counter()

for seq in sequences.values():

    for i in range(len(seq) - K + 1):

        kmer = seq[i:i+K]

        if set(kmer) <= BASES:

            kmers[canonical(kmer)] += 1

for row in targets.itertuples():

    seq = sequences[row.gene]

    position = int(row.local_position_1based) - 1

    region = seq[position-200:position+200]

    if not set(region) <= BASES:

        raise ValueError(
            f"Ambiguous bases near {row.gene} target"
        )

    local_kmers = [
        canonical(region[i:i+K])
        for i in range(len(region)-K+1)
    ]

    uniqueness = sum(
        kmers[k] == 1 for k in local_kmers
    ) / len(local_kmers)

    gc = (
        region.count("G")
        + region.count("C")
    ) / len(region)

    rows.append({
        "gene": row.gene,
        "clinvar_id": row.clinvar_id,
        "chromosome": row.chromosome,
        "genomic_position_1based": row.position,
        "local_position_1based":
            row.local_position_1based,
        "reference_length": len(seq),
        "ref": row.ref,
        "alt": row.alt,
        "reference_match": True,
        "local_gc_content": gc,
        "local_kmer_uniqueness": uniqueness,
        "qc_status": (
            "PASS" if uniqueness >= 0.99
            else "REVIEW"
        )
    })

with open(
    OUT / "combined_reference.fa", "w"
) as handle:

    for gene, seq in sequences.items():

        handle.write(f">{gene}\n")

        for i in range(0, len(seq), 70):

            handle.write(
                seq[i:i+70] + "\n"
            )

qc = pd.DataFrame(rows)

qc.to_csv(
    OUT / "reference_qc.csv",
    index=False
)

print("\nREFERENCE QC SUMMARY")
print(qc.to_string(index=False))

print(
    "\nCombined reference bases:",
    sum(map(len, sequences.values()))
)
