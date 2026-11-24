"""Research-only fetal genotype likelihoods from synthetic molecular counts."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import betabinom
from scipy.special import logsumexp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def prior(maternal,paternal=None,q=.02):
    m=maternal/2;p=q if paternal is None else paternal/2
    return np.array([(1-m)*(1-p),m*(1-p)+(1-m)*p,m*p])
def probabilities(maternal,ff,error=.001):
    return error+(1-2*error)*((1-ff)*maternal/2+ff*np.arange(3)/2)
def posterior(alt,n,maternal,ff,paternal=None,q=.02,rho=.0005,ff_sd=0):
    pr=prior(maternal,paternal,q);grid=np.clip(ff+ff_sd*np.linspace(-3,3,13),.0001,.5);weights=-.5*np.linspace(-3,3,13)**2;weights-=logsumexp(weights)
    p=np.array([probabilities(maternal,f) for f in grid]);k=1/rho-1
    ll=logsumexp(betabinom.logpmf(alt,n,p*k,(1-p)*k)+weights[:,None],axis=0)
    with np.errstate(divide='ignore'):lp=ll+np.log(pr)
    return np.exp(lp-logsumexp(lp))
def run(out,seed=42):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(seed);rows=[];rho=.0005;k=1/rho-1
    # Illustrative single-variant recessive loci, not real patient genotypes or clinical variant annotations.
    for gene in ['CFTR','HBB']:
        for ff in [.02,.04,.08,.12,.2]:
            for n in [500,2000,10000]:
                for truth in range(3):
                    for rep in range(25):
                        p=probabilities(1,ff)[truth];alt=rng.binomial(n,rng.beta(p*k,(1-p)*k));est=float(np.clip(rng.normal(ff,.005),.001,.5))
                        for mode,paternal in [('carrier_known',1),('paternal_unknown',None)]:
                            post=posterior(alt,n,1,est,paternal,ff_sd=.005);call=int(np.argmax(post)) if max(post)>=.95 else -1
                            rows.append(dict(gene=gene,true_ff=ff,estimated_ff=est,molecules=n,truth=truth,alt=alt,mode=mode,call=call,p0=post[0],p1=post[1],p2=post[2]))
    df=pd.DataFrame(rows);df.to_csv(out/'benchmark.csv',index=False);df['correct']=df.call==df.truth;df['no_call']=df.call==-1;df['tp']=(df.truth==2)&(df.call==2);df['fp']=(df.truth!=2)&(df.call==2);df['affected']=df.truth==2;df['unaffected']=df.truth!=2
    s=df.groupby(['mode','true_ff','molecules']).agg(accuracy_all=('correct','mean'),no_call_rate=('no_call','mean'),tp=('tp','sum'),fp=('fp','sum'),affected=('affected','sum'),unaffected=('unaffected','sum')).reset_index();s['affected_detection_rate']=s.tp/s.affected;s['false_positive_rate']=s.fp/s.unaffected;s.to_csv(out/'performance.csv',index=False)
    fig,ax=plt.subplots();
    for (mode,n),g in s.groupby(['mode','molecules']):ax.plot(g.true_ff,g.affected_detection_rate,marker='o',label=f'{mode}, {n}')
    ax.set(xlabel='True fetal fraction',ylabel='Affected-genotype detection rate',ylim=(0,1.05));ax.legend(fontsize=7);fig.tight_layout();fig.savefig(out/'detection_performance.png');plt.close(fig);print(s.to_string(index=False))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',default='results/single_gene');p.add_argument('--seed',type=int,default=42);a=p.parse_args();run(a.out,a.seed)
