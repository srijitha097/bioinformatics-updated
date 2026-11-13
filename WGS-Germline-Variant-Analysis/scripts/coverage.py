#!/usr/bin/env python3
"""Summarize depth at every position in a target interval."""
import csv
import sys
from pathlib import Path

def summarize(depth_path):
    depths = []
    with open(depth_path, newline='') as inp:
        for line in inp:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError('samtools depth row expected CHROM POS DP')
            depths.append(int(parts[2]))
    if not depths:
        raise ValueError('Empty coverage file; check region/contig mismatch')
    n = len(depths)
    return {'positions': n, 'mean_depth': sum(depths)/n,
            'pct_ge_10x': 100*sum(d>=10 for d in depths)/n,
            'pct_ge_20x': 100*sum(d>=20 for d in depths)/n,
            'pct_zero': 100*sum(d==0 for d in depths)/n}

if __name__ == '__main__':
    stats=summarize(sys.argv[1]); target=Path(sys.argv[2]); target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('w',newline='') as file:
        writer=csv.writer(file,delimiter='\t'); writer.writerow(['metric','value'])
        for key,value in stats.items(): writer.writerow([key, f'{value:.3f}' if isinstance(value,float) else value])
    print(stats)
