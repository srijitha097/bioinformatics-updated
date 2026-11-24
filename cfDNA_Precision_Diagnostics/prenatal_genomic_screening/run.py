"""SNP-based fetal fraction and chromosome dosage simulation."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, brentq
from scipy.special import xlogy,xlog1py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def estimate_ff(alt,depth,error=.001):
    # Informative SNPs: maternal AA, fetal AB, known in this simulation.
    def nll(f):
        p=error+(1-2*error)*f/2
        return -float(np.sum(xlogy(alt,p)+xlog1py(depth-alt,-p)))
    fit=minimize_scalar(nll,bounds=(0,.5),method='bounded');f=fit.x;target=fit.fun+1.920729
    lo=brentq(lambda x:nll(x)-target,0,f) if nll(0)>target else 0
    hi=brentq(lambda x:nll(x)-target,f,.5) if nll(.5)>target else .5
    return f,lo,hi

def run(out,seed=42):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(seed)
    chrom=np.repeat([1,2,3,4,5,13,18,21,22],100);gc=rng.uniform(.3,.7,len(chrom));bias=np.exp(-9*(gc-.45)**2);mapq=rng.uniform(.85,1,len(chrom));base=bias*mapq
    def profile(ff,tri,depth):
        w=base*np.where(chrom==tri,1+ff/2,1);c=rng.multinomial(depth,w/w.sum());corrected=c/base
        return np.array([corrected[chrom==k].sum()/corrected[np.isin(chrom,[1,2,3,4,5])].sum() for k in [13,18,21]])
    records=[]
    for depth in [100000,500000,1000000]:
        ref=np.array([profile(0,0,depth) for _ in range(100)]);mu=ref.mean(0);sd=ref.std(0,ddof=1)
        for ff in [.01,.02,.04,.08,.12,.2]:
            for tri in [0,13,18,21]:
                for rep in range(15):
                    snp_depth=rng.poisson(depth/10000,100);alt=rng.binomial(snp_depth,.001+.998*ff/2);est,lo,hi=estimate_ff(alt,snp_depth);z=(profile(ff,tri,depth)-mu)/sd
                    call='inconclusive' if lo<.04 else str([13,18,21][int(np.argmax(z))]) if max(z)>=3 else '0'
                    records.append(dict(depth=depth,true_ff=ff,estimated_ff=est,ff_lower=lo,ff_upper=hi,truth=str(tri),call=call,z13=z[0],z18=z[1],z21=z[2]))
    df=pd.DataFrame(records);df.to_csv(out/'benchmark.csv',index=False);df['correct']=df.truth==df.call;df['inconclusive']=df.call=='inconclusive';df['called']=~df.inconclusive
    s=df.groupby(['depth','true_ff']).agg(accuracy_all=('correct','mean'),no_call_rate=('inconclusive','mean'),n=('truth','size'),n_called=('called','sum'),n_correct=('correct','sum')).reset_index();s['accuracy_called']=s.n_correct/s.n_called.replace(0,np.nan);s.to_csv(out/'performance.csv',index=False)
    fig,ax=plt.subplots();ax.scatter(df.true_ff,df.estimated_ff,s=4,alpha=.2);ax.plot([0,.22],[0,.22],'k--');ax.set(xlabel='True fetal fraction',ylabel='Estimated fetal fraction');fig.tight_layout();fig.savefig(out/'fetal_fraction.png');plt.close(fig)
    fig,ax=plt.subplots();
    for d,g in s.groupby('depth'):ax.plot(g.true_ff,g.accuracy_all,marker='o',label=str(d))
    ax.set(xlabel='True fetal fraction',ylabel='Correct calls / all samples',ylim=(0,1.05));ax.legend(title='Reads');fig.tight_layout();fig.savefig(out/'screening_performance.png');plt.close(fig);print(s.to_string(index=False))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',default='results/prenatal');p.add_argument('--seed',type=int,default=42);a=p.parse_args();run(a.out,a.seed)
