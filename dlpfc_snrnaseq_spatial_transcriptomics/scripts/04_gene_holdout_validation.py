#!/usr/bin/env python3
"""Independent held-out gene expression reconstruction; avoids validation on training genes."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy import io,sparse
from scipy.optimize import nnls
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];d=ROOT/'data/processed';out=ROOT/'results';rng=np.random.default_rng(31)
ref=pd.read_csv(d/'reference_signatures.csv').set_index('gene');genes=pd.read_csv(d/'genes.csv')['gene'].astype(str).tolist()
Y=io.mmread(d/'spots.mtx');Y=Y.tocsr() if sparse.issparse(Y) else sparse.csr_matrix(Y)
common=[g for g in genes if g in ref.index]
if len(common)<60:raise ValueError('Need >=60 shared genes for holdout validation')
A=ref.loc[common].to_numpy(float);A=A/np.maximum(A.sum(axis=0),1)*1e6
Y=Y[[genes.index(g) for g in common],:].toarray().astype(float)
Y=Y/np.maximum(Y.sum(axis=0,keepdims=True),1)*1e6
ind=rng.permutation(len(common));held=ind[:max(20,len(common)//5)];train=ind[len(held):]
# Scaling factors held fixed across both matrices; train on 80% genes and evaluate on held-out 20%.
W=np.vstack([nnls(A[train,:],Y[train,i])[0] for i in range(Y.shape[1])]).T
prediction=A[held,:]@W
# Spearman is evaluated per spot, across held-out genes; this is predictive correspondence,
# not a clinical or biological ground truth metric.
rows=[]
for i in range(Y.shape[1]):
 observed=Y[held,i]; pred=prediction[:,i]
 rho=spearmanr(observed,pred).statistic if np.std(observed)>0 and np.std(pred)>0 else np.nan
 rows.append({'spot_number':i,'heldout_gene_spearman':rho})
tab=pd.DataFrame(rows);tab.to_csv(out/'heldout_gene_validation.csv',index=False)
print('Heldout gene reconstruction: median Spearman=',float(tab.heldout_gene_spearman.median()))
