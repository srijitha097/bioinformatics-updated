"""Cross-check normalized records against a separately intersected BED mask."""

import gzip
from bisect import bisect_right
from collections import Counter, defaultdict
from pathlib import Path

from ranking import save

data = Path("exome_gene_panel_cnv_analysis/data")


def mask(paths):
    events = defaultdict(Counter)
    for path in paths:
        spans = defaultdict(list)
        for line in path.read_text().splitlines():
            c, l, r = line.split()[:3]
            c = c.removeprefix("chr")
            start, end = int(l), int(r)
            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval: {line}")
            spans[c].append((start, end))
        for c, intervals in spans.items():
            merged = []
            for start, end in sorted(intervals):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            for start, end in merged:
                events[c][start] += 1
                events[c][end] -= 1
    boundaries = {}
    for c, counts in events.items():
        depth = 0
        edges = []
        for pos in sorted(counts):
            after = depth + counts[pos]
            if (depth >= len(paths)) != (after >= len(paths)):
                edges.append(pos)
            depth = after
        boundaries[c] = edges
    return boundaries


def calls(path, regions):
    sets = {"SNV": set(), "indel": set()}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c, pos, _, ref, alt, _, flt, _, fmt, sample = line.rstrip().split("\t")
            p = int(pos) - 1
            edges = regions.get(c, [])
            i = bisect_right(edges, p)
            if not i % 2 or p + len(ref) > edges[i] or set(ref + alt) - set("ACGT"):
                continue
            gt = (
                sample.split(":")[fmt.split(":").index("GT")]
                .replace("|", "/")
                .split("/")
            )
            if sorted(gt) not in (["0", "1"], ["1", "1"]):
                continue
            if len(ref) == len(alt) > 1:
                continue
            sets["SNV" if len(ref) == len(alt) else "indel"].add(
                (c, p, ref, alt, *sorted(gt))
            )
    return sets


def main():
    out = []
    for scope in ("exome", "expanded panel"):
        paths = [data / "verify_confident.bed", data / "verify_wes.bed"]
        if scope == "expanded panel":
            paths.append(Path("comparative_genomics/sequencing/panel_spans.bed"))
        regions = mask(paths)
        n = sum(
            sum(b - a for a, b in zip(edges[::2], edges[1::2]))
            for edges in regions.values()
        )
        t, q = [
            calls(
                next(Path("comparative_genomics/data").glob(prefix + "*.vcf.gz")),
                regions,
            )
            for prefix in ("GIAB", "converted")
        ]
        for kind in t:
            out.append(
                dict(
                    scope=scope,
                    kind=kind,
                    confident_bases=n,
                    tp=len(t[kind] & q[kind]),
                    fp=len(q[kind] - t[kind]),
                    fn=len(t[kind] - q[kind]),
                )
            )
    save("sequencing/normalized_crosscheck.tsv", out)


if __name__ == "__main__":
    main()
