"""Exact normalized allele/genotype comparisons on the GRCh37 autosomal capture."""

import os
import re
import subprocess
import sys
from pathlib import Path

from ranking import save

P = Path("exome_gene_panel_cnv_analysis")
sys.path.insert(0, str(P.resolve()))
from bisect import bisect_right

from run import intervals

d = P / "data"
bc = os.environ.get("BCFTOOLS", "bcftools")
extra = set(
    "ATM BRCA1 BRCA2 MLH1 MSH2 MSH6 PALB2 PMS2 PTEN TP53 APC MUTYH CDH1 STK11 SMAD4 BMPR1A VHL RET MEN1 RB1 NF1 NF2 TSC1 TSC2 FLCN FH SDHB SDHC SDHD WT1".split()
)
panel = Path("comparative_genomics/sequencing/panel_spans.bed")
panel.write_text(
    "".join(
        line + "\n"
        for line in (d / "Omim_Genes.bed").read_text().splitlines()
        if line.split()[3] in extra
    )
)
regions = [intervals(d / f) for f in ("verify_confident.bed", "verify_wes.bed")]
genes = intervals(panel)


def inside(bed, c, start, end):
    if c not in bed:
        return False
    starts, ends = bed[c]
    i = bisect_right(starts, start) - 1
    return i >= 0 and ends[i] >= end


def records(file):
    args = [
        bc,
        "query",
        "-R",
        str(d / "verify_wes.bed"),
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT\t%DP\t%GQ]\n",
        str(file),
    ]
    output = subprocess.check_output(args, text=True)
    calls = {}
    for line in output.splitlines():
        c, p, r, a, f, gt, dp, gq = line.split("\t")
        pos = int(p) - 1
        if (
            not re.fullmatch("[ACGT]+", a)
            or not re.fullmatch("[ACGT]+", r)
            or (len(a) == len(r) > 1)
        ):
            continue
        if not all(inside(b, c, pos, pos + len(r)) for b in regions):
            continue
        if not re.fullmatch("[01][/|][01]", gt) or "1" not in gt:
            continue
        key = (c, int(p), r, a)
        calls[key] = (
            gt.count("1"),
            f,
            float(dp) if dp != "." else -1,
            float(gq) if gq != "." else -1,
        )
    return calls


truth = records(next(d.glob("GIAB*.vcf.gz")))
query = records(next(d.glob("converted_ARUP*.vcf.gz")))
table = []
errors = []
for scope in ("exome", "expanded panel"):
    for kind in ("SNV", "indel"):

        def wanted(k):
            c, p, r, a = k
            return (len(r) == len(a) == 1) == (kind == "SNV") and (
                scope == "exome" or inside(genes, c, p - 1, p - 1 + len(r))
            )

        t = {k: v[0] for k, v in truth.items() if wanted(k) and v[1] == "PASS"}
        for dp, gq in [(0, 0), (10, 0), (20, 0), (30, 0), (0, 30), (20, 30)]:
            q = {
                k: v[0]
                for k, v in query.items()
                if wanted(k)
                and v[1] == "PASS"
                and (not dp or v[2] >= dp)
                and (not gq or v[3] >= gq)
            }
            tp = len(t.items() & q.items())
            fp = len(q) - tp
            fn = len(t) - tp
            table.append(
                dict(
                    scope=scope,
                    kind=kind,
                    DP=dp,
                    GQ=gq,
                    truth=len(t),
                    query=len(q),
                    tp=tp,
                    fp=fp,
                    fn=fn,
                    precision=tp / len(q) if q else None,
                    recall=tp / len(t) if t else None,
                    f1=2 * tp / (len(t) + len(q)) if t or q else None,
                )
            )
            if scope == "exome" and dp == gq == 0:
                for k in sorted(t.keys() | q.keys()):
                    if t.get(k) == q.get(k):
                        continue
                    reason = (
                        "genotype mismatch"
                        if k in t and k in q
                        else (
                            "query filtered"
                            if k in t and k in query
                            else "absent eligible query"
                            if k in t
                            else "absent eligible truth"
                        )
                    )
                    errors.append(
                        dict(
                            chrom=k[0],
                            pos=k[1],
                            ref=k[2],
                            alt=k[3],
                            kind=kind,
                            reason=reason,
                            truth_dosage=t.get(k),
                            query_dosage=q.get(k),
                            query_record=query.get(k),
                        )
                    )
save("sequencing/metrics.tsv", table)
save("sequencing/discordances.tsv", errors)
