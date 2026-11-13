#!/usr/bin/env python3
"""Generate an accurate computational summary without inventing biological discoveries."""
import argparse
from pathlib import Path
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--real',action='store_true');a=p.parse_args()
root=Path(__file__).resolve().parents[1]
r=root/'results'; manifest=r/'pseudobulk_manifest.csv'; summary=[]
if not manifest.exists():raise SystemExit('No pseudobulk manifest; run pseudobulk analysis first')
m=pd.read_csv(manifest)
summary.append('# Project 1 — analysis run summary')
summary.append('')
summary.append('**Input:** '+('REAL public GSE254569 measurements (exploratory subset)' if a.real else 'SYNTHETIC TEST DATA — NOT BIOLOGICAL FINDINGS'))
summary.append('')
summary.append(f'Cell types with pseudobulk data (not all necessarily pass coverage QC): {len(m)}')
summary.append('')
summary.append('| Cell type | Case donors before coverage QC | Control donors before coverage QC | Genes with FDR<0.05 |')
summary.append('|---|---:|---:|---:|')
for _,row in m.iterrows():
 f=r/'deseq2'/f"{row['folder']}_deseq2.csv"
 if f.exists():
  tab=pd.read_csv(f);sig=int((tab.padj<0.05).sum())
 else:sig='not run'
 summary.append(f"| {row['cell_type']} | {int(row['n_case_donors'])} | {int(row['n_control_donors'])} | {sig} |")
summary.extend(['','Donor counts displayed above precede the DESeq2 per-sample nucleus coverage filter; see cell-type _summary.txt files for donors retained in each model.','', 'These counts are descriptive, not biological claims. Examine gene identities, covariates, confounders, sample selection, and statistical validity before interpretation.','',
 'See `data/processed/real_data_provenance.txt` for dataset, layer, and downsampling if real data were used.'])
(r/'RUN_SUMMARY.md').write_text('\n'.join(summary)+'\n')
print('Created',r/'RUN_SUMMARY.md')
