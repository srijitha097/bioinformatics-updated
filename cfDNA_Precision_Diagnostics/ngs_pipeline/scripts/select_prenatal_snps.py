import json
from collections import Counter
from pathlib import Path

import pandas as pd
import pysam

REFERENCE = "ngs_pipeline/data/reference/hbb_region.fa"
VARIANTS = "ngs_pipeline/data/snps/chr11_variants.json"
OUTPUT = "ngs_pipeline/data/snps/informative_snps.tsv"

REGION_START = 5200001
K = 75
BASES = set("ACGT")


def reverse_complement(seq):
    return seq.translate(
        str.maketrans("ACGT", "TGCA")
    )[::-1]


def canonical(seq):
    return min(seq, reverse_complement(seq))


with pysam.FastaFile(REFERENCE) as fasta:
    contig = fasta.references[0]
    sequence = fasta.fetch(contig).upper()

with open(VARIANTS) as handle:
    variants = json.load(handle)

kmers = Counter(
    canonical(sequence[i:i+K])
    for i in range(len(sequence)-K+1)
    if set(sequence[i:i+K]) <= BASES
)

candidates = []
excluded = Counter()

for variant in variants:

    if variant.get("assembly_name", "GRCh38") != "GRCh38":
        excluded["assembly"] += 1
        continue

    if str(variant.get("seq_region_name")) != "11":
        excluded["chromosome"] += 1
        continue

    if variant.get("strand", 1) != 1:
        excluded["strand"] += 1
        continue

    rsid = str(variant.get("id", ""))

    if not rsid.startswith("rs"):
        excluded["identifier"] += 1
        continue

    alleles = variant.get("alleles", [])

    if (
        len(alleles) != 2
        or any(a not in BASES for a in alleles)
    ):
        excluded["not_biallelic_snp"] += 1
        continue

    genomic = int(variant["start"])

    if genomic != int(variant["end"]):
        excluded["not_single_base"] += 1
        continue

    position = genomic - REGION_START

    if position < 250 or position >= len(sequence)-250:
        excluded["reference_boundary"] += 1
        continue

    ref = sequence[position]

    if alleles[0] != ref:
        excluded["reference_mismatch"] += 1
        continue

    alt = alleles[1]

    if ref == alt:
        excluded["identical_alleles"] += 1
        continue

    region = sequence[position-200:position+200]

    if not set(region) <= BASES:
        excluded["ambiguous_reference"] += 1
        continue

    gc = (
        region.count("G") + region.count("C")
    ) / len(region)

    if not 0.35 <= gc <= 0.65:
        excluded["gc_content"] += 1
        continue

    local_kmers = [
        canonical(region[i:i+K])
        for i in range(len(region)-K+1)
    ]

    uniqueness = sum(
        kmers[k] == 1 for k in local_kmers
    ) / len(local_kmers)

    if uniqueness != 1.0:
        excluded["mapping_uniqueness"] += 1
        continue

    candidates.append({
        "rsid": rsid,
        "chromosome": "11",
        "contig": contig,
        "genomic_position_1based": genomic,
        "local_position_1based": position+1,
        "ref": ref,
        "alt": alt,
        "gc_content": gc,
        "local_uniqueness": uniqueness,
        "maternal_genotype": "AA",
        "fetal_genotype": "AB"
    })

candidates.sort(
    key=lambda row: (
        abs(row["gc_content"] - 0.5),
        row["local_position_1based"]
    )
)

selected = []

for candidate in candidates:

    position = candidate["local_position_1based"]

    if any(
        abs(position - row["local_position_1based"]) < 500
        for row in selected
    ):
        continue

    selected.append(candidate)

    if len(selected) == 30:
        break

if not selected:
    raise RuntimeError("No suitable SNP markers found")

df = pd.DataFrame(selected).sort_values(
    "local_position_1based"
)

Path(OUTPUT).parent.mkdir(
    parents=True, exist_ok=True
)

df.to_csv(OUTPUT, sep="\t", index=False)

print("Total annotated variants:", len(variants))
print("Candidate SNPs:", len(candidates))
print("Selected SNP markers:", len(df))

print("\nExclusion summary:")
for reason, count in excluded.most_common():
    print(reason, count)

print("\nSelected SNPs:")
print(df.to_string(index=False))

print("\nOutput:", OUTPUT)
