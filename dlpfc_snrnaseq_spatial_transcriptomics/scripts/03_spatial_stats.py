#!/usr/bin/env python3
"""Within-section Moran's I with coordinate-permutation p-values + FDR."""
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.neighbors import NearestNeighbors
from statsmodels.stats.multitest import multipletests
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];d=ROOT/'data/processed';out=ROOT/'results'
meta=pd.read_csv(d/'spots.csv');prop=pd.read_csv(out/'spot_celltype_proportions.csv')
tab=meta.merge(prop,on='spot_id',how='inner',validate='one_to_one')
if len(tab)!=len(meta):raise ValueError('Spot metadata/proportions mismatch')
cts=[c for c in prop.columns if c!='spot_id'];rng=np.random.default_rng(42);rows=[]
def moran(x,edges):
    z=x-x.mean();den=(z*z).sum()
    if den==0:return float('nan')
    src,dst=edges;return len(x)*float((z[src]*z[dst]).sum())/(len(src)*den)
for sample,g in tab.groupby('sample_id'):
    if len(g)<12:continue
    xy=g[['x','y']].to_numpy(); k=min(7,len(g))
    neigh=NearestNeighbors(n_neighbors=k).fit(xy).kneighbors(return_distance=False)
    src=np.repeat(np.arange(len(g)),neigh.shape[1]); dst=neigh.ravel(); edges=(src,dst)
    for ct in cts:
      values=g[ct].to_numpy(float);value=moran(values,edges)
      if not np.isfinite(value):continue
      null=np.array([moran(rng.permutation(values),edges) for _ in range(199)])
      pv=(1+np.sum(np.abs(null)>=abs(value)))/(len(null)+1)
      rows.append({'sample_id':sample,'cell_type':ct,'n_spots':len(g),'moran_I':value,'permutation_p':pv})
res=pd.DataFrame(rows)
if res.empty:raise ValueError('No valid section-level Moran statistics')
res['fdr']=multipletests(res.permutation_p,method='fdr_bh')[1];res.to_csv(out/'moran_by_section.csv',index=False)
if 'layer' in tab and tab.layer.notna().any():
   av=tab.groupby(['sample_id','layer'])[cts].mean().reset_index();av.to_csv(out/'layer_celltype_estimates.csv',index=False)
   overall=tab.groupby('layer')[cts].mean().T
   overall.plot(kind='bar',figsize=(11,4));plt.ylabel('Average predicted proportion');plt.title('Exploratory layer patterns (section confounding possible)');plt.tight_layout()
   plt.savefig(out/'figures/03_layer_celltype_composition.png',dpi=160);plt.close()
print('Saved Moran I and exploratory layer results:',out)
