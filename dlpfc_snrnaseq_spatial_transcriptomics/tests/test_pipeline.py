from pathlib import Path
import numpy as np,pandas as pd
from scipy import io
ROOT=Path(__file__).resolve().parents[1]
def test_demo_spatial_input_and_proportions():
    d=ROOT/'data/processed';r=ROOT/'results';metadata=pd.read_csv(d/'spots.csv')
    genes=pd.read_csv(d/'genes.csv');X=io.mmread(d/'spots.mtx')
    assert X.shape==(len(genes),len(metadata))
    prop=pd.read_csv(r/'spot_celltype_proportions.csv')
    assert len(prop)==len(metadata) and len(metadata.sample_id.unique())==3
    assert (prop.drop(columns=['spot_id']).to_numpy()>=-1e-10).all()
    assert np.allclose(prop.drop(columns=['spot_id']).sum(axis=1),1,atol=1e-6)
def test_spatial_fdr_and_validation():
    r=ROOT/'results';v=pd.read_csv(r/'heldout_gene_validation.csv');m=pd.read_csv(r/'moran_by_section.csv')
    assert len(v)>100 and v.heldout_gene_spearman.notna().sum()>100
    assert m['fdr'].between(0,1).all()
