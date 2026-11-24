import json
from pathlib import Path

import pandas as pd
import pysam

ROOT = Path("results/ngs/aneuploidy")
ALIGNMENT = ROOT / "alignment"
OUT = ROOT / "benchmark"

OUT.mkdir(parents=True, exist_ok=True)

manifest = pd.read_csv(
    "ngs_pipeline/data/aneuploidy/sample_manifest.csv"
)

rows = []

for sample in manifest.itertuples():

    bam_path = ALIGNMENT / f"{sample.sample}.bam"

    if not bam_path.exists():
        raise FileNotFoundError(bam_path)

    total = 0
    mapped = 0
    paired = 0
    high_mapq = 0
    mapq_sum = 0
    fragments = 0
    fragment_length_sum = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:

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

    mapping_rate = mapped / total if total else 0

    pairing_rate = paired / total if total else 0

    high_mapq_rate = (
        high_mapq / mapped if mapped else 0
    )

    status = (
        "PASS"
        if total == sample.reads
        and mapping_rate >= 0.95
        and pairing_rate >= 0.95
        and high_mapq_rate >= 0.95
        else "FAIL"
    )

    rows.append({
        "sample": sample.sample,
        "cohort": sample.cohort,
        "condition": sample.condition,
        "true_ff": sample.true_ff,
        "expected_reads": sample.reads,
        "total_primary_reads": total,
        "mapped_reads": mapped,
        "mapping_pct": mapping_rate * 100,
        "properly_paired_pct": pairing_rate * 100,
        "high_mapq_pct": high_mapq_rate * 100,
        "mean_mapq": (
            mapq_sum / mapped if mapped else 0
        ),
        "proper_fragments": fragments,
        "mean_fragment_length": (
            fragment_length_sum / fragments
            if fragments else 0
        ),
        "qc_status": status
    })

    print(
        sample.sample,
        "| Reads:", total,
        "| Mapping:", round(mapping_rate * 100, 3),
        "| Pairing:", round(pairing_rate * 100, 3),
        "| QC:", status,
        flush=True
    )

df = pd.DataFrame(rows)

df.to_csv(
    OUT / "alignment_qc_summary.csv",
    index=False
)

metrics = {
    "samples": int(len(df)),
    "total_primary_reads": int(
        df.total_primary_reads.sum()
    ),
    "total_mapped_reads": int(
        df.mapped_reads.sum()
    ),
    "total_proper_fragments": int(
        df.proper_fragments.sum()
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
    OUT / "alignment_qc_metrics.json", "w"
) as handle:

    json.dump(metrics, handle, indent=2)

print("\nBATCH QC METRICS")
print(json.dumps(metrics, indent=2))

print("\nQC summary saved to:", OUT)
