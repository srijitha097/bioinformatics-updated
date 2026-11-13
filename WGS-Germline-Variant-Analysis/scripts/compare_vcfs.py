#!/usr/bin/env python3
"""Educational exact allele overlap on normalized VCFs inside GIAB confident BED.
Not a substitute for hap.py / vcfeval, which account for complex variants.
"""
import argparse
import bisect
import gzip
import json
from collections import defaultdict
from pathlib import Path

def open_text(path):
    return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path,encoding='utf-8')

def load_bed(path, chrom, start, end):
    intervals=[]
    with open_text(path) as handle:
        for row in handle:
            if not row.strip() or row.startswith('#'): continue
            fields=row.split()
            if fields[0] != chrom: continue
            # BED zero-based, half-open. Convert to 1-based closed interval.
            left,right=int(fields[1])+1,int(fields[2])
            left,right=max(start,left),min(end,right)
            if left<=right:intervals.append((left,right))
    intervals.sort()
    merged=[]
    for left,right in intervals:
        if merged and left<=merged[-1][1]+1:merged[-1]=(merged[-1][0],max(right,merged[-1][1]))
        else:merged.append((left,right))
    return merged

def alleles(path, chrom, start, end, intervals):
    points=[v[0] for v in intervals]
    keys=set()
    with open_text(path) as handle:
        for row in handle:
            if row.startswith('#'):continue
            fields=row.rstrip('\n').split('\t')
            if len(fields)<8:raise ValueError('Malformed VCF row')
            c,pos,_,ref,alts,_,fil,_=fields[:8]
            if c!=chrom:continue
            pos=int(pos)
            if pos<start or pos>end:continue
            ix=bisect.bisect_right(points,pos)-1
            if ix<0:continue
            # Include full reference span only inside high-confidence interval.
            if pos+len(ref)-1>intervals[ix][1]:continue
            if fil not in ('PASS','.'):
                continue
            for alt in alts.split(','):
                if alt!='.' and not alt.startswith('<') and '[' not in alt and ']' not in alt:
                    keys.add((c,pos,ref,alt))
    return keys

def compare(calls,truth,bed,region):
    chrom,coords=region.split(':')
    start,end=map(int,coords.replace(',','').split('-'))
    intervals=load_bed(bed,chrom,start,end)
    if not intervals: raise ValueError(f'No confident BED intervals for {region}')
    c=alleles(calls,chrom,start,end,intervals)
    t=alleles(truth,chrom,start,end,intervals)
    tp=len(c&t);fp=len(c-t);fn=len(t-c)
    precision=tp/(tp+fp) if tp+fp else None
    recall=tp/(tp+fn) if tp+fn else None
    f1=2*precision*recall/(precision+recall) if precision is not None and recall is not None and precision+recall else None
    return {'region':region,'truth_bed_covered_bases':sum(b-a+1 for a,b in intervals),
            'call_alleles':len(c),'truth_alleles':len(t),'TP':tp,'FP':fp,'FN':fn,
            'precision':precision,'recall':recall,'f1':f1,
            'method':'Exact normalized REF/ALT allele matching within GIAB confident BED, PASS/. only',
            'limitations':'Not genotype-aware or haplotype-aware. Indel representation may still differ. Region is 1Mb of WGS source, not entire WGS.'}

if __name__=='__main__':
    pa=argparse.ArgumentParser()
    for opt in ('calls','truth','bed','region','output'):pa.add_argument('--'+opt,required=True)
    x=pa.parse_args()
    result=compare(x.calls,x.truth,x.bed,x.region)
    Path(x.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
