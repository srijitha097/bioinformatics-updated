#!/usr/bin/env python3
"""Read PLINK 1.9 outputs, validate columns, and create original QC graphics."""
import argparse
import html
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def read_table(path):
    return pd.read_csv(path, sep=r'\s+', na_values=['NA','nan','-nan'])

def association(path):
    table=read_table(path)
    if 'TEST' in table.columns:
        table=table[table['TEST']=='ADD'].copy()
    required={'CHR','BP','SNP','P'}
    if not required.issubset(table):raise ValueError(f'Missing fields {required-set(table)} in {path}')
    table=table.copy()
    table['CHR']=pd.to_numeric(table['CHR'],errors='coerce')
    table['BP']=pd.to_numeric(table['BP'],errors='coerce')
    table['P']=pd.to_numeric(table['P'],errors='coerce')
    return table.dropna(subset=['CHR','BP','P']).query('P > 0 and P <= 1').sort_values(['CHR','BP'])

def plot_manhattan(df, outfile):
    fig,ax=plt.subplots(figsize=(12,4.3)); offset=0; ticks=[]; labels=[]
    for chrom,grp in df.groupby('CHR',sort=True):
        positions=grp['BP'].to_numpy()+offset
        ax.scatter(positions,-np.log10(grp['P']),s=5,alpha=0.6)
        ticks.append((positions.min()+positions.max())/2)
        labels.append(str(int(chrom)))
        offset+=grp['BP'].max()+1000000
    ax.set_xticks(ticks,labels,fontsize=7)
    ax.set_xlabel('Chromosome');ax.set_ylabel('-log10(p)')
    ax.set_title('HapMap educational GWAS — simulated phenotype')
    ax.axhline(-np.log10(5e-8),linestyle='--',linewidth=1,color='gray')
    fig.tight_layout();fig.savefig(outfile,dpi=160);plt.close(fig)

def plot_qq(df,outfile):
    p=np.sort(df['P'].to_numpy());n=len(p)
    if not n:raise ValueError('No valid association p-values')
    expected=-np.log10(np.arange(1,n+1)/(n+1))
    observed=-np.log10(p)
    fig,ax=plt.subplots(figsize=(5,5));ax.scatter(expected,observed,s=5,alpha=0.65)
    top=max(expected.max(),observed.max(),1)
    ax.plot([0,top],[0,top],'k--',linewidth=1)
    ax.set(xlabel='Expected -log10(p)',ylabel='Observed -log10(p)',
           title='QQ — simulated phenotype')
    fig.tight_layout();fig.savefig(outfile,dpi=160);plt.close(fig)

def plot_pca(path,outfile):
    df=pd.read_csv(path,sep=r'\s+',header=None)
    if df.shape[1]<5:raise ValueError('PLINK eigenvec requires FID IID PC1 PC2 PC3')
    fig,ax=plt.subplots(figsize=(6,5));ax.scatter(df.iloc[:,2],df.iloc[:,3],s=22,alpha=0.75)
    ax.set(xlabel='PC1',ylabel='PC2',title='PCA of LD-pruned HapMap SNPs')
    fig.tight_layout();fig.savefig(outfile,dpi=160);plt.close(fig)

def build_report(directory):
    d=Path(directory)
    before_samples=sum(1 for line in open(d/'raw.fam') if line.strip())
    after_samples=sum(1 for line in open(d/'filtered.fam') if line.strip())
    before_variants=sum(1 for line in open(d/'raw.bim') if line.strip())
    after_variants=sum(1 for line in open(d/'filtered.bim') if line.strip())
    if not after_samples or not after_variants:raise ValueError('No variants/samples survived filters')
    g1=association(d/'gwas_unadjusted.assoc')
    g2=association(d/'gwas_adjusted.assoc.logistic')
    if g1.empty:raise ValueError('Unadjusted GWAS has no valid p-values')
    plot_manhattan(g1,d/'manhattan.png');plot_qq(g1,d/'qq.png')
    if not g2.empty:
        plot_manhattan(g2,d/'manhattan_adjusted.png')
        plot_qq(g2,d/'qq_adjusted.png')
    plot_pca(d/'pca.eigenvec',d/'pca.png')
    imiss=read_table(d/'initial.imiss');lmiss=read_table(d/'initial.lmiss')
    mean_sample_missing=float(imiss['F_MISS'].mean())
    mean_variant_missing=float(lmiss['F_MISS'].mean())
    data={'before_samples':before_samples,'after_samples':after_samples,
          'before_variants':before_variants,'after_variants':after_variants,
          'mean_sample_missingness':mean_sample_missing,
          'mean_snp_missingness':mean_variant_missing,
          'unadjusted_tested_snps':int(len(g1)),'adjusted_tested_snps':int(len(g2)),
          'unadjusted_smallest_p':float(g1['P'].min()),
          'adjusted_smallest_p':float(g2['P'].min()) if not g2.empty else None,
          'phenotype':'SIMULATED tutorial binary disease status; not a real phenotype'}
    (d/'qc_metrics.json').write_text(json.dumps(data,indent=2)+'\n')
    rows=''.join(f'<tr><td>{html.escape(k)}</td><td>{html.escape(str(v))}</td></tr>' for k,v in data.items())
    report=f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>HapMap QC/GWAS</title>
<style>body{{font:16px system-ui;max-width:980px;margin:2em auto;padding:1em}}img{{max-width:100%}}td,th{{border:1px solid #ddd;padding:9px}}table{{border-collapse:collapse;width:100%}}.note{{background:#fff3d4;padding:1em}}</style>
</head><body><h1>Genotyping QC, population PCA and GWAS</h1>
<p class="note"><strong>Demo only.</strong> HapMap genotypes are real; phenotype was SIMULATED. Only 89 participants; not powered for disease inference.</p>
<h2>QC & retained data</h2><table>{rows}</table><h2>Manhattan — unadjusted</h2><img src="manhattan.png">
<h2>QQ — unadjusted</h2><img src="qq.png">
<h2>Manhattan — PCA-adjusted</h2><img src="manhattan_adjusted.png">
<h2>QQ — PCA-adjusted</h2><img src="qq_adjusted.png"><h2>PCA</h2><img src="pca.png">
<p>Review full adjusted results at <code>gwas_adjusted.assoc.logistic</code>. PLINK HWE and frequency filters are for teaching, not universally appropriate QC thresholds.</p></body></html>"""
    (d/'report.html').write_text(report)
    return data

if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('--results',default='results')
    a=pa.parse_args();print(json.dumps(build_report(a.results),indent=2))
