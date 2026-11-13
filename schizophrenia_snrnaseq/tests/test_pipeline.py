from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse
ROOT=Path(__file__).resolve().parents[1]
def test_synthetic_inputs_and_integer_counts():
    d=ROOT/'data/processed'; x=sparse.load_npz(d/'counts.npz'); cells=pd.read_csv(d/'cells.csv');genes=pd.read_csv(d/'genes.csv')
    assert x.shape==(len(cells),len(genes))
    assert x.data.min()>=0 and np.allclose(x.data,np.rint(x.data))
    assert cells.donor_id.nunique()==24
    assert set(cells.diagnosis)=={'case','control'}
def test_pseudobulk_equals_raw_sum():
    d=ROOT/'data/processed';x=sparse.load_npz(d/'counts.npz');cells=pd.read_csv(d/'cells.csv')
    manifest=pd.read_csv(ROOT/'results/pseudobulk_manifest.csv')
    for _,row in manifest.iterrows():
        counts=pd.read_csv(ROOT/'results/pseudobulk'/row['folder']/'counts.csv',index_col=0)
        sub=cells.cell_type==row['cell_type']
        assert int(counts.to_numpy().sum())==int(x[sub.to_numpy()].sum())
