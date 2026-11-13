#!/usr/bin/env python3
"""Optional overrepresentation testing on a supplied GMT with the tested-gene universe.

Not automatically run because MSigDB/Reactome GMT downloads may have their own terms.
Gene identifiers in the GMT must match the exported DE result identifiers exactly.
"""
from pathlib import Path
import argparse
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--gmt',required=True);p.add_argument('--fdr',type=float,default=0.05)
a=p.parse_args()
sets={}
with open(a.gmt,encoding='utf-8') as handle:
 for line in handle:
  s=line.rstrip('\n').split('\t')
  if len(s)>=4:sets[s[0]]=set(s[2:])
out=R/'results/enrichment';out.mkdir(exist_ok=True,parents=True)
files=list(sorted((R/'results/deseq2').glob('*_deseq2.csv')))
if not files:raise SystemExit('No DE results found. Run Rscript scripts/02_deseq2.R first.')
for fn in files:
 df=pd.read_csv(fn)
 eligible=df[df['pvalue'].notna()].copy()
 universe=set(eligible.gene.astype(str))
 sig=set(eligible.loc[(eligible.padj<a.fdr)&(eligible.log2FoldChange.abs()>=0.5),'gene'].astype(str))
 rows=[]
 for term,all_genes in sets.items():
  genes=all_genes & universe
  if not 10<=len(genes)<=500:continue
  overlap=sig & genes
  pv=hypergeom.sf(len(overlap)-1,len(universe),len(genes),len(sig)) if overlap else 1.0
  rows.append({'term':term,'pvalue':pv,'overlap_n':len(overlap),'term_n':len(genes),'genes':';'.join(sorted(overlap))})
 if not rows:
  print('No eligible gene sets for',fn.name,'; check gene identifier namespace')
  continue
 tab=pd.DataFrame(rows);tab['padj']=multipletests(tab.pvalue,method='fdr_bh')[1]
 tab.sort_values('padj').to_csv(out/(fn.stem+'_enrichment.csv'),index=False)
 print(fn.name,'sets tested',len(tab),'significant',(tab.padj<a.fdr).sum())
