"""Simulate targeted UMI read families and compare read and molecule counts."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASES='ACGT'
def mutate(s, rate, rng):
    return ''.join(rng.choice([b for b in BASES if b!=x]) if rng.random()<rate else x for x in s)
def collapse(rows):
    groups=defaultdict(Counter)
    for locus, umi, allele in rows: groups[locus][umi]+=1
    mapping={}
    for locus, counts in groups.items():
        for u in sorted(counts,key=lambda x:(-counts[x],x)):
            candidates=[v for v in counts if counts[v]>counts[u] and counts[v]>=2*counts[u]-1 and sum(a!=b for a,b in zip(u,v))==1]
            parent=max(candidates,key=lambda v:(counts[v],v)) if candidates else u
            mapping[locus,u]=mapping.get((locus,parent),parent)
    families=defaultdict(Counter)
    for locus,u,a in rows: families[locus,mapping[locus,u]][a]+=1
    consensus=[(k,c.most_common(1)[0][0]) for k,c in families.items() if c.total()>=2 and c.most_common(1)[0][1]/c.total()>=0.75]
    return consensus

def run(out,seed=42,n=1200):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(seed); summaries=[]; all_reads=[]
    for vaf in [0,0.005,0.01,0.03,0.05,0.1]:
        for rep in range(5):
            rows=[];truth=[];observed=set()
            for i in range(n):
                locus=i%12;u=''.join(rng.choice(list(BASES),12));a=int(rng.random()<vaf);truth.append(a)
                copies=rng.poisson(5 if a==0 else 8)
                if copies: observed.add(i)
                for _ in range(copies):
                    allele=a if rng.random()>=0.008 else 1-a
                    rows.append((locus,mutate(u,0.002,rng),allele))
            cons=collapse(rows);raw=np.mean([r[2] for r in rows]); corrected=np.mean([a for _,a in cons]); actual=np.mean(truth)
            summaries.append(dict(target_vaf=vaf,replicate=rep,true_vaf=actual,raw_vaf=raw,consensus_vaf=corrected,input_molecules=n,observed_molecules=len(observed),reads=len(rows),consensus_families=len(cons),family_count_ratio=len(cons)/n,raw_abs_error=abs(raw-actual),consensus_abs_error=abs(corrected-actual)))
            if rep==0: all_reads.extend((vaf,*r) for r in rows)
    df=pd.DataFrame(summaries);df.to_csv(out/'benchmark.csv',index=False);pd.DataFrame(all_reads,columns=['target_vaf','locus','umi','allele']).to_csv(out/'example_reads.csv.gz',index=False)
    fig,ax=plt.subplots();ax.scatter(df.true_vaf,df.raw_vaf,label='Raw reads');ax.scatter(df.true_vaf,df.consensus_vaf,label='Consensus');ax.plot([0,.13],[0,.13],'k--');ax.set(xlabel='True molecular VAF',ylabel='Estimated VAF');ax.legend();fig.tight_layout();fig.savefig(out/'vaf_concordance.png');plt.close(fig)
    print(df[['raw_abs_error','consensus_abs_error','family_count_ratio']].mean().to_string())
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',default='results/molecular');p.add_argument('--seed',type=int,default=42);p.add_argument('--molecules',type=int,default=1200);a=p.parse_args();run(a.out,a.seed,a.molecules)
