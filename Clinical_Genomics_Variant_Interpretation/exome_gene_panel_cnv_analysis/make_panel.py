"""Intersect hereditary cancer gene spans with the source exome capture design."""

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

from run import intervals, overlap

GENES = {
    "ATM",
    "BRCA1",
    "BRCA2",
    "MLH1",
    "MSH2",
    "MSH6",
    "PALB2",
    "PMS2",
    "PTEN",
    "TP53",
}


def make_panel(genes, capture, out):
    selected = {}
    for line in genes.read_text().splitlines():
        chrom, start, end, gene = line.split()[:4]
        if gene in GENES:
            if gene in selected:
                raise ValueError(f"Duplicate gene span for {gene}")
            selected[gene] = f"{chrom}\t{start}\t{end}\n"
    if selected.keys() != GENES:
        raise ValueError(f"Missing gene spans for {GENES - selected.keys()}")
    with TemporaryDirectory() as temp:
        spans = Path(temp) / "spans.bed"
        spans.write_text("".join(selected.values()))
        captured = intervals(capture)
        overlaps = list(overlap(intervals(spans), captured))
    if not overlaps:
        raise ValueError("No selected genes overlap the capture design")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        "".join(f"chr{chrom}\t{start}\t{end}\n" for chrom, start, end in overlaps)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("genes", "capture", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    make_panel(args.genes, args.capture, args.out)
