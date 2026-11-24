from collections import Counter
import pysam

REFERENCE = "ngs_pipeline/data/reference/hbb_region.fa"
K = 75

def reverse_complement(seq):
    return seq.translate(
        str.maketrans("ACGT", "TGCA")
    )[::-1]

def canonical(seq):
    return min(seq, reverse_complement(seq))

with pysam.FastaFile(REFERENCE) as fasta:
    contig = fasta.references[0]
    sequence = fasta.fetch(contig).upper()

kmers = Counter(
    canonical(sequence[i:i+K])
    for i in range(len(sequence)-K+1)
    if set(sequence[i:i+K]) <= set("ACGT")
)

candidates = []

for position in range(500, len(sequence)-500, 50):
    region = sequence[position-200:position+200]

    if not set(region) <= set("ACGT"):
        continue

    gc = (region.count("G") + region.count("C")) / len(region)

    if not 0.35 <= gc <= 0.65:
        continue

    local_kmers = [
        canonical(region[i:i+K])
        for i in range(len(region)-K+1)
    ]

    uniqueness = sum(
        kmers[k] == 1 for k in local_kmers
    ) / len(local_kmers)

    candidates.append(
        (uniqueness, position, gc)
    )

candidates.sort(
    key=lambda x: (-x[0], abs(x[2]-0.5))
)

if not candidates:
    raise RuntimeError("No suitable target found")

print("Top candidate targets:")

for uniqueness, position, gc in candidates[:10]:
    print(
        f"Local position: {position+1}, "
        f"Genomic position: {5200001+position}, "
        f"Uniqueness: {uniqueness:.2%}, "
        f"GC: {gc:.2%}"
    )

best = candidates[0]

with open(
    "ngs_pipeline/data/annotations/selected_target.txt", "w"
) as handle:
    handle.write(str(best[1]) + "\n")

print("\nSelected local position:", best[1]+1)
print("Selected genomic position:", 5200001+best[1])
print("Uniqueness:", best[0])
