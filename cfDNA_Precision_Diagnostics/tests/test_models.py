import numpy as np
from molecular_quantification.run import collapse
from prenatal_genomic_screening.run import estimate_ff
from single_gene_screening.run import prior, posterior, probabilities

def test_umi_error_correction():
    rows=[(0,'AAAAAAAAAAAA',1)]*5+[(0,'AAAAAAAAAAAT',0)]+[(1,'AAAAAAAAAAAA',0)]*3
    assert sorted(a for _,a in collapse(rows))==[0,1]
def test_tied_consensus_rejected():
    assert collapse([(0,'AAAA',0),(0,'AAAA',1)])==[]
def test_fetal_fraction():
    f,lo,hi=estimate_ff(np.full(100,51),np.full(100,1000));assert abs(f-.1)<.002 and lo<f<hi
def test_inheritance():
    assert np.allclose(prior(1,1),[.25,.5,.25]);assert np.allclose(prior(0,0),[1,0,0])
def test_posterior():
    p=posterior(5500,10000,1,.1,1,rho=.00001);assert np.isclose(p.sum(),1) and p.argmax()==2
