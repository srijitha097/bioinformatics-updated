#!/usr/bin/env python3
"""Create reproducible *synthetic* cell-by-gene UMI counts for pipeline testing."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse
ROOT = Path(__file__).resolve().parents[1]
out = ROOT / 'data' / 'processed'; out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(23)
cell_types = ['Excitatory', 'Inhibitory', 'Astrocyte', 'Oligodendrocyte']
genes = ['SLC17A7','GAD1','GAD2','GFAP','AQP4','MBP','PLP1','SYN1','SNAP25','GRIN2A','NRXN1','DLG4'] + [f'GENE_{i:04d}' for i in range(288)]
rows, count_chunks = [], []
for d in range(24):
    diag = 'case' if d >= 12 else 'control'
    age = int(rng.integers(35, 75)); sex = 'F' if d % 2 else 'M'; pmi = round(float(rng.uniform(8, 28)), 1)
    for ct in cell_types:
        n = int(rng.integers(28, 43))
        mu = np.full(len(genes), 0.09)
        if ct == 'Excitatory': mu[[0,7,8,9,10,11]] = 1.7
        if ct == 'Inhibitory': mu[[1,2,8]] = 1.8
        if ct == 'Astrocyte': mu[[3,4]] = 2.5
        if ct == 'Oligodendrocyte': mu[[5,6]] = 2.2
        if diag == 'case' and ct == 'Excitatory': mu[12:26] *= 2.7
        if diag == 'case' and ct == 'Astrocyte': mu[40:50] *= 2.1
        mu *= rng.lognormal(0, 0.1)
        mat = rng.poisson(mu * rng.lognormal(0,0.25,(n,1))).astype(np.int32)
        count_chunks.append(sparse.csr_matrix(mat))
        for j in range(n): rows.append({'cell_id':f'D{d:02d}_{ct}_{j:03d}', 'donor_id':f'D{d:02d}', 'diagnosis':diag, 'cell_type':ct, 'age':age, 'sex':sex, 'pmi':pmi})
X = sparse.vstack(count_chunks,format='csr')
pd.DataFrame(rows).to_csv(out/'cells.csv',index=False)
pd.DataFrame({'gene':genes}).to_csv(out/'genes.csv',index=False)
sparse.save_npz(out/'counts.npz',X)
(out/'SYNTHETIC_DATA_DO_NOT_REPORT_AS_DISCOVERY.txt').write_text('Synthetic simulated expression, not real human measurements. Use only to test software.\n')
print(f'SYNTHETIC DEMO: {X.shape[0]} nuclei x {X.shape[1]} genes -> {out}')
