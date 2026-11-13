#!/usr/bin/env python3
import argparse,csv,html,json
from pathlib import Path

def render(benchmark,coverage):
    fields=''.join(f'<tr><td>{html.escape(str(k))}</td><td>{html.escape(str(v))}</td></tr>' for k,v in benchmark.items())
    cov=''.join(f'<tr><td>{html.escape(str(r["metric"]))}</td><td>{html.escape(str(r["value"]))}</td></tr>' for r in coverage)
    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>HG002 QC report</title>
<style>body{{font:16px system-ui;max-width:950px;margin:3rem auto;padding:1rem;color:#172330}}
table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #d3dce2;padding:0.55rem;text-align:left}}th{{background:#e8eff5}}
.warning{{background:#fff0ce;padding:1rem;border-left:4px solid #cf8618}}</style></head><body>
<h1>HG002 chr20 germline regional QC</h1><p class="warning">Research and engineering demonstration. The data are a regional slice of WGS.
Exact allele overlap is NOT a clinical validation or haplotype-aware benchmark.</p>
<h2>Sequencing coverage</h2><table><tr><th>Metric</th><th>Value</th></tr>{cov}</table>
<h2>GIAB allele overlap</h2><table><tr><th>Metric</th><th>Value</th></tr>{fields}</table></body></html>"""

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--benchmark',required=True);p.add_argument('--coverage',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();b=json.loads(Path(a.benchmark).read_text())
    with open(a.coverage,newline='') as f:c=list(csv.DictReader(f,delimiter='\t'))
    Path(a.out).write_text(render(b,c))
