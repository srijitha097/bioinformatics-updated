#!/usr/bin/env python3
"""Audit GSE254569's real H5AD schema before selecting donor, diagnosis, and count layer."""
import argparse
from pathlib import Path
import anndata
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser()
p.add_argument('--h5ad', default=str(ROOT/'data/raw/GSE254569_adata_RNA.h5ad'))
p.add_argument('--sample', type=int, default=5)
a=p.parse_args()
adata=anndata.read_h5ad(a.h5ad, backed='r')
try:
    print('\n========== REAL DATA AUDIT ==========')
    print('Cells × genes:', adata.shape)
    print('obs fields:', list(adata.obs.columns))
    print('var fields:', list(adata.var.columns))
    print('layers:', list(adata.layers.keys()), 'raw available:', adata.raw is not None)
    print('\nFirst metadata rows (may contain sensitive donor identifiers):')
    print(adata.obs.head(a.sample).to_string(max_cols=45))
    for c in adata.obs.columns:
        if any(s in c.lower() for s in ('diagnos','disease','disorder','condition','group','classif','celltype','cell_type','sample','donor')):
            series=adata.obs[c].astype(str)
            print('\n',repr(c), 'unique:',series.nunique(), 'examples:',series.value_counts().head(12).to_dict())
    print('\nIMPORTANT: verify donor IDs, case/control categories, and that chosen matrix contains genuine integer UMIs.')
finally:
    adata.file.close()
