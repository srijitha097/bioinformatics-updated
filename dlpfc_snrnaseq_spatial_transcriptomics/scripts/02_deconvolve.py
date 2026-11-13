#!/usr/bin/env python3
"""NNLS cell-type proportion estimates using independent snRNA-seq reference signatures."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy import io,sparse
from scipy.optimize import nnls
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];d=ROOT/'data/processed';out=ROOT/'results';out.mkdir(exist_ok=True)
ref=pd.read_csv(d/'reference_signatures.csv').set_index('gene');genes=pd.read_csv(d/'genes.csv')['gene'].astype(str)
meta=pd.read_csv(d/'spots.csv'); Y=io.mmread(d/'spots.mtx');Y=Y.tocsr() if sparse.issparse(Y) else sparse.csr_matrix(Y)
if Y.shape!=(len(genes),len(meta)):raise ValueError(f'Matrix orientation mismatch: {Y.shape} vs {len(genes)} genes and {len(meta)} spots')
if np.any(Y.data<0):raise ValueError('Visium counts must be non-negative')
if meta['spot_id'].duplicated().any():raise ValueError('Spot IDs must be unique')
common=genes[genes.isin(ref.index)].tolist()
if len(common)<30:raise ValueError('Too few shared genes; check gene identifier harmonization')
ix=pd.Index(genes).get_indexer(common)
Y2=Y[ix,:].toarray().astype(float)
Y2=1e6*Y2/np.maximum(Y2.sum(axis=0),1)
A=ref.loc[common].to_numpy(dtype=float)
A/=np.maximum(A.sum(axis=0,keepdims=True),1)
A*=1e6
# Remove features absent from the reference and spatial observations.
good=(A.max(axis=1)>0)&(Y2.sum(axis=1)>0)
A=A[good];Y2=Y2[good,:]
if len(A)<30:raise ValueError('Insufficient informative genes')
# NNLS reference weight estimation, followed by simplex normalization.
W=np.vstack([nnls(A,Y2[:,i])[0] for i in range(Y2.shape[1])]);W/=np.maximum(W.sum(axis=1,keepdims=True),1e-12)
prop=pd.DataFrame(W,columns=ref.columns);prop.insert(0,'spot_id',meta.spot_id.values)
prop.to_csv(out/'spot_celltype_proportions.csv',index=False)
(out/'figures').mkdir(exist_ok=True)
avg=prop[ref.columns].mean().sort_values(ascending=False)
avg.plot(kind='bar',figsize=(9,4));plt.ylabel('Estimated proportion');plt.title('Mean inferred cell-type composition (NNLS)');plt.tight_layout()
plt.savefig(out/'figures/01_mean_composition.png',dpi=160);plt.close()
for ct in avg.index[:4]:
 fig,axs=plt.subplots(1,len(meta.sample_id.unique()),figsize=(4.3*len(meta.sample_id.unique()),4),squeeze=False)
 for j,sample in enumerate(sorted(meta.sample_id.unique())):
    keep=(meta.sample_id==sample).to_numpy(); im=axs[0,j].scatter(meta.loc[keep,'x'],meta.loc[keep,'y'],c=prop.loc[keep,ct],s=10,cmap='viridis',vmin=0,vmax=1)
    axs[0,j].set_title(str(sample));axs[0,j].set_aspect('equal');axs[0,j].invert_yaxis()
 fig.colorbar(im,ax=axs.ravel().tolist(),label=f'{ct} estimated fraction',fraction=.015)
 fig.suptitle(f'Spatial distribution: {ct}');fig.savefig(out/f'figures/02_spatial_{ct.replace("/","_").replace(" ","_")}.png',dpi=160,bbox_inches='tight');plt.close(fig)
print('Deconvolution complete:',W.shape[0],'spots x',W.shape[1],'cell types. NOTE: NNLS is a baseline model, not single-cell ground truth.')
