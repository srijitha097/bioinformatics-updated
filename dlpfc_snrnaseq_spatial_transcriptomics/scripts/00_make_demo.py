#!/usr/bin/env python3
"""Make a spatially structured *synthetic* reference and Visium-like spot matrix."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse,io
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'data/processed';out.mkdir(parents=True,exist_ok=True)
rng=np.random.default_rng(101)
cts=['Excitatory','Inhibitory','Astrocyte','Oligodendrocyte','Microglia']
markers=['SLC17A7','SATB2','GAD1','GAD2','GFAP','AQP4','MBP','PLP1','P2RY12','CX3CR1','SYN1','SNAP25','GRIN2A','NRXN1']
genes=markers+[f'GENE_{k:04d}' for k in range(386)]
ref=rng.gamma(1.2,0.3,(len(genes),len(cts)))
for c,idx in enumerate([[0,1,10,11,12,13],[2,3,10,11],[4,5],[6,7],[8,9]]):ref[idx,c]+=25
ref=ref/ref.sum(axis=0,keepdims=True)*1e6
pd.DataFrame(ref,index=genes,columns=cts).to_csv(out/'reference_signatures.csv',index_label='gene')
cols=[];counts=[];truth=[]
for sample in range(3):
 for yy in range(18):
  for xx in range(18):
   layer= 'upper' if yy<6 else ('middle' if yy<12 else 'deep')
   p=np.array([.51,.23,.09,.12,.05]) if layer=='upper' else (np.array([.32,.39,.11,.13,.05]) if layer=='middle' else np.array([.30,.13,.10,.40,.07]))
   p=rng.dirichlet(50*p)
   depth=rng.integers(1200,2600)
   spot=rng.poisson(depth*(ref/1e6@p))
   counts.append(spot);truth.append(p)
   cols.append({'spot_id':f'S{sample+1}_{xx:02d}_{yy:02d}','sample_id':f'S{sample+1}','x':xx*11+(sample*2),'y':yy*11,'layer':layer})
pd.DataFrame(cols).to_csv(out/'spots.csv',index=False)
pd.DataFrame({'gene':genes}).to_csv(out/'genes.csv',index=False)
io.mmwrite(out/'spots.mtx',sparse.csr_matrix(np.array(counts,dtype=np.int32).T))
pd.DataFrame(truth,columns=cts).to_csv(out/'DEMO_TRUE_PROPORTIONS.csv',index=False)
(out/'SYNTHETIC_DATA_DO_NOT_REPORT_AS_DISCOVERY.txt').write_text('Synthetic expression and spatial locations. Biological conclusions are invalid.\n')
print('SYNTHETIC demo:',len(cols),'spots,',len(genes),'genes in',out)
