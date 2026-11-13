#!/usr/bin/env python3
"""Aggregate raw UMI counts by independent donor AND annotated cell type."""
from pathlib import Path
import re
import numpy as np, pandas as pd
from scipy import sparse
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
data=ROOT/'data'/'processed'; out=ROOT/'results'; out.mkdir(exist_ok=True)
cells=pd.read_csv(data/'cells.csv'); genes=pd.read_csv(data/'genes.csv')['gene'].astype(str).tolist()
X=sparse.load_npz(data/'counts.npz').tocsr()
assert X.shape==(len(cells),len(genes)), 'Dimension mismatch for cells / genes / counts'
needed={'donor_id','cell_type','diagnosis'}
assert needed.issubset(cells.columns), f'Missing metadata columns: {needed-set(cells.columns)}'
assert cells['donor_id'].notna().all() and cells['cell_type'].notna().all(), 'Missing donor or cell type'
assert set(cells['diagnosis'].unique())=={'case','control'}, 'Require case and control only'
assert np.all(X.data>=0) and np.allclose(X.data,np.rint(X.data)), 'Pseudobulk requires integer raw UMI counts'
meta_fields=[c for c in ['diagnosis','age','sex','pmi','rin','lib_batch'] if c in cells]
check=cells.groupby('donor_id')[meta_fields].nunique(dropna=False)
if (check>1).any().any(): raise ValueError('Conflicting donor metadata within donor_id')
composition=cells.groupby(['donor_id','diagnosis','cell_type'],observed=True).size().reset_index(name='n_nuclei')
composition.to_csv(out/'donor_celltype_counts.csv',index=False)
piv=composition.pivot_table(index='donor_id',columns='cell_type',values='n_nuclei',fill_value=0)
(piv.div(piv.sum(axis=1),axis=0)).plot(kind='bar',stacked=True,figsize=(11,4),width=0.86)
plt.ylabel('Fraction of nuclei');plt.title('Nuclei composition by donor (descriptive, not DE)')
plt.tight_layout(); (out/'figures').mkdir(exist_ok=True);plt.savefig(out/'figures'/'01_celltype_composition.png',dpi=160);plt.close()
log=[]; pb=out/'pseudobulk';pb.mkdir(exist_ok=True)
for ct in sorted(cells['cell_type'].unique()):
    ix=np.flatnonzero(cells['cell_type'].to_numpy()==ct)
    subset=cells.iloc[ix].copy(); groups=sorted(subset['donor_id'].unique())
    group_index=pd.Categorical(subset['donor_id'],categories=groups).codes
    A=sparse.csr_matrix((np.ones(len(ix),dtype=np.int32),(group_index,np.arange(len(ix)))),shape=(len(groups),len(ix)))
    agg=(A@X[ix,:]).tocsr()
    md=subset.drop_duplicates('donor_id').set_index('donor_id').loc[groups][meta_fields]
    md['n_nuclei']=subset.groupby('donor_id', observed=True).size().reindex(groups).to_numpy()
    md.index.name='sample_id'; md=md.reset_index()
    if min((md['diagnosis']=='case').sum(),(md['diagnosis']=='control').sum())<3:
        print(f'Skipping {ct}: fewer than 3 independent donors/group');continue
    safe=re.sub(r'[^A-Za-z0-9_.-]+','_',str(ct))[:60]
    target=pb/safe;target.mkdir(exist_ok=True)
    pd.DataFrame(agg.toarray().T,index=genes,columns=groups).to_csv(target/'counts.csv',index_label='gene')
    md.to_csv(target/'metadata.csv',index=False)
    log.append({'cell_type':ct,'folder':safe,'n_nuclei':len(ix),'n_donors':len(groups), 'n_case_donors':int((md['diagnosis']=='case').sum()),'n_control_donors':int((md['diagnosis']=='control').sum())})
pd.DataFrame(log).to_csv(out/'pseudobulk_manifest.csv',index=False)
if not log: raise SystemExit('No cell types had sufficient donors in both conditions; check annotations and filtering.')
print('Pseudobulk complete. NEXT: Rscript scripts/02_deseq2.R  (requires DESeq2).')
