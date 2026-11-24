import json
from pathlib import Path

import pandas as pd
import pysam

DATA = Path("ngs_pipeline/data/single_gene")
OUT = Path("results/ngs/single_gene")
DEST = OUT / "benchmark"

DEST.mkdir(parents=True, exist_ok=True)

manifest = pd.read_csv(
    DATA / "sample_manifest.csv"
)

rows = []

for sample in manifest.itertuples():

    total = 0
    mapped = 0
    paired = 0
    high_mapq = 0
    mapq_sum = 0
    fragments = 0
    fragment_length_sum = 0

    bam_path = (
        OUT / "alignment" / f"{sample.sample}.bam"
    )

    with pysam.AlignmentFile(
        str(bam_path), "rb"
    ) as bam:

        for read in bam.fetch(until_eof=True):

            if (
                read.is_secondary
                or read.is_supplementary
            ):
                continue

            total += 1

            if read.is_unmapped:
                continue

            mapped += 1

            mapq_sum += read.mapping_quality

            if read.mapping_quality >= 30:
                high_mapq += 1

            if read.is_proper_pair:

                paired += 1

                if read.is_read1:

                    fragments += 1

                    fragment_length_sum += abs(
                        read.template_length
                    )

    mapping_pct = (
        mapped / total * 100 if total else 0
    )

    pairing_pct = (
        paired / total * 100 if total else 0
    )

    high_mapq_pct = (
        high_mapq / mapped * 100 if mapped else 0
    )

    mean_mapq = (
        mapq_sum / mapped if mapped else 0
    )

    mean_fragment_length = (
        fragment_length_sum / fragments
        if fragments else 0
    )

    status = (
        "PASS"
        if (
            mapping_pct >= 95
            and pairing_pct >= 95
            and high_mapq_pct >= 95
            and total == sample.reads
        )
        else "FAIL"
    )

    rows.append({
        "sample": sample.sample,
        "true_ff": sample.true_ff,
        "fetal_genotype": sample.fetal_genotype,
        "molecules_per_gene":
            sample.molecules_per_gene,
        "expected_reads": sample.reads,
        "total_primary_reads": total,
        "mapped_reads": mapped,
        "mapping_pct": mapping_pct,
        "properly_paired_pct": pairing_pct,
        "high_mapq_pct": high_mapq_pct,
        "mean_mapq": mean_mapq,
        "proper_read_pairs": fragments,
        "mean_fragment_length":
            mean_fragment_length,
        "qc_status": status
    })

    print(
        sample.sample,
        "| Mapping:", round(mapping_pct, 3),
        "| Pairing:", round(pairing_pct, 3),
        "| QC:", status,
        flush=True
    )

df = pd.DataFrame(rows)

df.to_csv(
    DEST / "alignment_qc_summary.csv",
    index=False
)

metrics = {
    "samples": int(len(df)),
    "expected_reads": int(
        df.expected_reads.sum()
    ),
    "total_primary_reads": int(
        df.total_primary_reads.sum()
    ),
    "total_mapped_reads": int(
        df.mapped_reads.sum()
    ),
    "total_proper_read_pairs": int(
        df.proper_read_pairs.sum()
    ),
    "mean_mapping_pct": float(
        df.mapping_pct.mean()
    ),
    "mean_properly_paired_pct": float(
        df.properly_paired_pct.mean()
    ),
    "mean_high_mapq_pct": float(
        df.high_mapq_pct.mean()
    ),
    "mean_mapq": float(
        df.mean_mapq.mean()
    ),
    "mean_fragment_length": float(
        df.mean_fragment_length.mean()
    ),
    "qc_passed": int(
        (df.qc_status == "PASS").sum()
    ),
    "qc_failed": int(
        (df.qc_status == "FAIL").sum()
    )
}

with open(
    DEST / "alignment_qc_metrics.json", "w"
) as handle:

    json.dump(metrics, handle, indent=2)

print("\nALIGNMENT QC METRICS")
print(json.dumps(metrics, indent=2))

print("\nQC SUMMARY BY MOLECULAR DEPTH")

print(
    df.groupby("molecules_per_gene").agg(
        samples=("sample", "size"),
        mean_mapping_pct=("mapping_pct", "mean"),
        mean_pairing_pct=("properly_paired_pct", "mean"),
        mean_mapq=("mean_mapq", "mean"),
        qc_passed=("qc_status", lambda x: (x == "PASS").sum())
    ).to_string()
)

print("\nQC summary saved to:", DEST)
