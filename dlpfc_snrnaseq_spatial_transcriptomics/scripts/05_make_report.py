#!/usr/bin/env python3
"""Summarize actual outputs, visibly distinguishing simulation from real data."""
import argparse
from pathlib import Path
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--real',action='store_true');args=p.parse_args()
root=Path(__file__).resolve().parents[1];r=root/'results'
prop=pd.read_csv(r/'spot_celltype_proportions.csv')
meta=pd.read_csv(root/'data/processed/spots.csv')
m=pd.read_csv(r/'moran_by_section.csv')
v=pd.read_csv(r/'heldout_gene_validation.csv')
lines=['# Project 2 — analysis run summary','', '**Input:** '+('REAL LIBD public human brain data' if args.real else 'SYNTHETIC SOFTWARE TEST DATA'),'',
       f'Analyzed {len(meta):,} spots from {meta.sample_id.nunique()} sections.',
       f'Inferred {len(prop.columns)-1} broad cell-type components using NNLS.',
       f'Per-section Moran statistics: {len(m)}; FDR < 0.05: {int((m.fdr<0.05).sum())}.',
       f'Holdout gene reconstruction median Spearman: {v.heldout_gene_spearman.median():.3f}.','',
       'NNLS proportions are model estimates and are not observed cell counts. Moran I measures spatial autocorrelation, not cell–cell communication.',
       'Do not claim psychiatric disease effects from neurotypical-only samples.']
(r/'RUN_SUMMARY.md').write_text('\n'.join(lines)+'\n')
print('Created',r/'RUN_SUMMARY.md')
