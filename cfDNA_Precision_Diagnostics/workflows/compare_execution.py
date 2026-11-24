import json
from pathlib import Path

import pandas as pd
import pysam
import hashlib

ROOT = Path("results/nextflow")

native = ROOT / "smoke_test/alignment"
docker = ROOT / "docker_test/alignment"


def summarize(path):

    primary = 0
    mapped = 0
    paired = 0
    mapq_sum = 0

    digest = hashlib.sha256()

    with pysam.AlignmentFile(str(path), "rb") as bam:

        for read in bam.fetch(until_eof=True):

            if read.is_secondary or read.is_supplementary:
                continue

            primary += 1

            if not read.is_unmapped:
                mapped += 1
                mapq_sum += read.mapping_quality

            if read.is_proper_pair:
                paired += 1

            fields = (
                read.query_name,
                read.flag,
                read.reference_name,
                read.reference_start,
                read.cigarstring,
                read.mapping_quality,
                read.template_length
            )

            digest.update(
                (str(fields) + "\n").encode()
            )

    return {
        "primary_reads": primary,
        "mapped_reads": mapped,
        "properly_paired_reads": paired,
        "mean_mapq": mapq_sum / mapped if mapped else 0,
        "ordered_alignment_signature": digest.hexdigest()
    }


native_files = sorted(native.glob("*.bam"))
docker_files = sorted(docker.glob("*.bam"))

assert len(native_files) == len(docker_files) == 3

assert {p.name for p in native_files} == {
    p.name for p in docker_files
}

rows = []

for path in native_files:

    a = summarize(path)
    b = summarize(docker / path.name)

    for metric in a:

        rows.append({
            "sample": path.stem,
            "metric": metric,
            "native_value": a[metric],
            "docker_value": b[metric],
            "exact_match": a[metric] == b[metric]
        })

df = pd.DataFrame(rows)

df.to_csv(
    ROOT / "execution_reproducibility.csv",
    index=False
)

numeric = df[
    df.metric != "ordered_alignment_signature"
]

signatures = df[
    df.metric == "ordered_alignment_signature"
]

metrics = {
    "samples_compared": len(native_files),
    "numeric_metric_comparisons": len(numeric),
    "numeric_exact_matches": int(
        numeric.exact_match.sum()
    ),
    "numeric_agreement_pct": float(
        numeric.exact_match.mean() * 100
    ),
    "alignment_signatures_compared": len(signatures),
    "matching_alignment_signatures": int(
        signatures.exact_match.sum()
    )
}

with open(
    ROOT / "execution_reproducibility.json", "w"
) as handle:

    json.dump(metrics, handle, indent=2)

print("\nALIGNMENT METRIC COMPARISON")
print(numeric.to_string(index=False))

print("\nALIGNMENT SIGNATURE COMPARISON")
print(
    signatures[["sample", "exact_match"]]
    .to_string(index=False)
)

print("\nREPRODUCIBILITY METRICS")
print(json.dumps(metrics, indent=2))
