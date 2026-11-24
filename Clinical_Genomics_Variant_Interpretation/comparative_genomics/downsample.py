"""Measure depth sensitivity using seeded alignment subsampling on fixed public slices."""

import os
import subprocess
import tempfile
from pathlib import Path

from ranking import save

sg = Path("exome_gene_panel_cnv_analysis")
sam = os.environ.get("SAMTOOLS", "samtools")
rows = []
with tempfile.TemporaryDirectory() as tmp:
    for bam in sorted((sg / "data").glob("SM-*-v1-chr20-downsampled.deduplicated.bam")):
        for fraction in (0.25, 0.5, 0.75, 1):
            for seed in [17, 83, 211] if fraction < 1 else [17]:
                run = bam
                if fraction < 1:
                    run = Path(tmp) / "slice.bam"
                    subprocess.run(
                        [
                            sam,
                            "view",
                            "--subsample",
                            str(fraction),
                            "--subsample-seed",
                            str(seed),
                            "-b",
                            "-o",
                            str(run),
                            str(bam),
                        ],
                        check=True,
                    )
                bed = sg / "results/ice_chr20.bed"
                raw = subprocess.check_output(
                    [
                        sam,
                        "depth",
                        "-a",
                        "-b",
                        str(bed),
                        "-s",
                        "-q",
                        "20",
                        "-Q",
                        "20",
                        "-G",
                        "SUPPLEMENTARY",
                        str(run),
                    ],
                    text=True,
                )
                depths = [int(s.split()[2]) for s in raw.splitlines()]
                n = sum(
                    int(b) - int(a)
                    for _, a, b in (s.split() for s in bed.read_text().splitlines())
                )
                assert len(depths) <= n
                rows.append(
                    dict(
                        sample=bam.name.split("-v1")[0],
                        fraction=fraction,
                        seed=seed,
                        bases=n,
                        mean_depth=sum(depths) / n,
                        fraction20=sum(v >= 20 for v in depths) / n,
                    )
                )
save("sequencing/alignment_subsampling.tsv", rows)
