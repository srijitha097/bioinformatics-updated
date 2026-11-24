"""Target-restricted, exact SNV/genotype comparison on single-sample normalized VCFs."""

import argparse
import bisect
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import op, report, sha

AUTOSOMES = {str(number) for number in range(1, 23)}


def contig(name):
    bare = name.removeprefix("chr")
    return bare if bare in AUTOSOMES | {"X", "Y"} else name


def intervals(path):
    d = defaultdict(list)
    with op(path) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            c, s, e, *_ = line.split()
            s, e = int(s), int(e)
            if s < 0 or e <= s:
                raise ValueError("Invalid BED interval")
            d[contig(c)].append((s, e))
    result = {}
    for c, arr in d.items():
        merged = []
        for s, e in sorted(arr):
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(e, merged[-1][1]))
            else:
                merged.append((s, e))
        result[c] = ([s for s, e in merged], [e for s, e in merged])
    if not result:
        raise ValueError("Empty BED")
    return result


def contains(d, c, p):
    c = contig(c)
    if c not in d:
        return False
    s, e = d[c]
    i = bisect.bisect_right(s, p) - 1
    return i >= 0 and p < e[i]


def vcf(path, targets, confident):
    variants = {}
    skipped = 0
    header = False
    with op(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                if len(line.rstrip().split("\t")) != 10:
                    raise ValueError("Expected single-sample VCF")
                header = True
                continue
            if line.startswith("#"):
                continue
            c, pos, _, ref, alt, _, flt, _, fmt, sample = line.rstrip().split("\t")
            pos = int(pos)
            c = contig(c)
            if not (contains(targets, c, pos - 1) and contains(confident, c, pos - 1)):
                continue
            if (
                flt != "PASS"
                or len(ref) != 1
                or len(alt) != 1
                or ref not in "ACGT"
                or alt not in "ACGT"
            ):
                skipped += 1
                continue
            gt = dict(zip(fmt.split(":"), sample.split(":"))).get("GT", ".")
            alleles = gt.replace("|", "/").split("/")
            if len(alleles) != 2 or any(a not in ("0", "1") for a in alleles):
                skipped += 1
                continue
            if alleles == ["0", "0"]:
                continue
            key = (c, pos, ref, alt)
            if key in variants:
                raise ValueError(f"Duplicate SNV {key}")
            variants[key] = "/".join(sorted(alleles))
    if not header:
        raise ValueError("Missing VCF header")
    return variants, skipped


def overlap(left, right):
    for chrom in sorted(left.keys() & right.keys()):
        starts, ends = left[chrom]
        other_starts, other_ends = right[chrom]
        i = j = 0
        while i < len(starts) and j < len(other_starts):
            start = max(starts[i], other_starts[j])
            end = min(ends[i], other_ends[j])
            if start < end:
                yield chrom, start, end
            if ends[i] <= other_ends[j]:
                i += 1
            else:
                j += 1


def compare(truth, query, targets, confidence, out, autosomes_only=False):
    t, c = intervals(targets), intervals(confidence)
    if autosomes_only:
        t = {chrom: spans for chrom, spans in t.items() if chrom in AUTOSOMES}
    target_bases = sum(
        end - start for starts, ends in t.values() for start, end in zip(starts, ends)
    )
    confident_target_bases = sum(end - start for _, start, end in overlap(t, c))
    if not confident_target_bases:
        raise ValueError("No target bases overlap benchmark confident regions")
    a, sa = vcf(truth, t, c)
    b, sb = vcf(query, t, c)
    if not a:
        raise ValueError(
            "No truth SNVs in intersected regions; check genome build and chromosome naming"
        )
    shared = a.keys() & b.keys()
    tp = sum(a[k] == b[k] for k in shared)
    fp = len(b) - tp
    fn = len(a) - tp
    report(
        out,
        dict(
            research_only=True,
            method="Exact biallelic PASS diploid SNV and genotype match; no indel or haplotype equivalence",
            target_bases=target_bases,
            confident_target_bases=confident_target_bases,
            autosomes_only=autosomes_only,
            truth_snvs=len(a),
            query_snvs=len(b),
            true_positive=tp,
            false_positive=fp,
            false_negative=fn,
            shared_allele_genotype_disagreements=len(shared) - tp,
            precision=tp / (tp + fp) if tp + fp else None,
            recall=tp / (tp + fn),
            f1=2 * tp / (2 * tp + fp + fn),
            skipped_truth_records=sa,
            skipped_query_records=sb,
            input_sha256={str(p): sha(p) for p in (truth, query, targets, confidence)},
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for k in ("truth", "query", "targets", "confidence", "out"):
        p.add_argument("--" + k, type=Path, required=True)
    p.add_argument("--autosomes-only", action="store_true")
    a = p.parse_args()
    compare(a.truth, a.query, a.targets, a.confidence, a.out, a.autosomes_only)
