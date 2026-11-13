#!/usr/bin/env python3
"""Import public human brain snRNA-seq with explicit, auditable metadata validation.

No RNA expression is silently converted to counts. If diagnosis cannot be linked to
nucleus donor IDs, generate donor_mapping_template.csv and STOP for human review.
"""
import argparse
import re
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]

def guess(columns, candidates):
    index = {str(c).casefold(): c for c in columns}
    return next((index[t.casefold()] for t in candidates if t.casefold() in index), None)

def key(s):
    return re.sub(r'[^a-z0-9]', '', str(s).casefold())

from import_helpers import diagnosis_group

def cli():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--h5ad',default=str(ROOT/'data/raw/GSE254569_adata_RNA.h5ad'))
    p.add_argument('--geo-metadata', default=str(ROOT/'data/raw/geo_donor_metadata.csv'))
    p.add_argument('--donor-map',help='Optional manually reviewed CSV with donor_id, diagnosis (case/control/other), and optionally age/sex/pmi')
    p.add_argument('--donor-col');p.add_argument('--celltype-col');p.add_argument('--diagnosis-col')
    p.add_argument('--geo-key-col', default='sample_title',help='Column in GEO SOFT metadata matching nucleus donor labels')
    p.add_argument('--geo-diagnosis-col', help='Column in GEO SOFT metadata containing diagnosis; if omitted, auto-scan common names')
    p.add_argument('--count-layer',default='auto',help='auto|X|raw|layer name; must be genuine integer UMI counts')
    p.add_argument('--max-nuclei-per-donor',type=int,default=150,help='Exploratory subsample only; use larger settings for robust final analysis')
    p.add_argument('--max-genes',type=int,default=8000,help='Top genes by pilot prevalence. Use 0 for all genes.')
    p.add_argument('--max-celltypes',type=int,default=6)
    p.add_argument('--sampling-mode', choices=['donor','donor_celltype'], default='donor_celltype', help='Use donor_celltype for a balanced, cell-type-aware analysis')
    p.add_argument('--seed',type=int,default=42)
    args=p.parse_args()
    if args.max_nuclei_per_donor<1 or args.max_genes<0 or args.max_celltypes<1: p.error('Invalid positive sampling limits')
    return args

def export_template(obs, destination):
    d=obs[['donor_id']].drop_duplicates().sort_values('donor_id').copy()
    d['diagnosis']='';d['age']='';d['sex']='';d['pmi']=''
    destination.parent.mkdir(parents=True,exist_ok=True)
    d.to_csv(destination,index=False)
    return destination

def main():
    args=cli();rng=np.random.default_rng(args.seed)
    ad=anndata.read_h5ad(args.h5ad,backed='r')
    try:
        obs=ad.obs.copy()
        donor=args.donor_col or guess(obs.columns,['donor_id','donor','sample_id','sample','sample_name','individual','subject_id','patient_id'])
        ct=args.celltype_col or guess(obs.columns,['cell_type_broad','cellType_broad','celltype_broad','cell_type','cellType','celltype','major_celltype','annotation','cluster_annotation'])
        if donor is None or ct is None:
            raise ValueError('Missing donor or cell-type column. Run: python scripts/00_inspect_h5ad.py; then pass --donor-col and --celltype-col.')
        obs['donor_id']=obs[donor].astype(str).str.strip()
        obs['cell_type']=obs[ct].astype(str).str.strip()
        if (obs.donor_id.str.lower().isin(['nan','','none'])).any():raise ValueError('Missing donor values; cannot proceed.')
        dia_col=args.diagnosis_col or guess(obs.columns,['diagnosis','classification','diagnostic_group','disease','disease_status','condition','Dx'])
        labels=None; label_source=None
        if args.donor_map:
            metadata=pd.read_csv(args.donor_map,dtype={'donor_id':str,'diagnosis':str})
            if not {'donor_id','diagnosis'}.issubset(metadata):raise ValueError('Donor map must include donor_id and diagnosis')
            if metadata.donor_id.duplicated().any():raise ValueError('Duplicate donor IDs in donor map')
            # Explicit map requires exact case/control/other labels to avoid misleading regex guesses.
            valid={'case','control','other'}
            if not metadata.diagnosis.dropna().isin(valid).all():raise ValueError('Donor-map diagnosis must be exactly case/control/other')
            merged=obs['donor_id'].map(metadata.set_index('donor_id')['diagnosis'])
            if merged.isna().any():raise ValueError('Donor map missing one or more observed donor IDs')
            obs['diagnosis']=merged
            for name in ('age','sex','pmi'):
                if name in metadata:
                    obs[name]=obs.donor_id.map(metadata.set_index('donor_id')[name])
            label_source='reviewed donor map'
        elif dia_col:
            obs['source_diagnosis']=obs[dia_col].astype(str)
            obs['diagnosis']=obs.source_diagnosis.map(diagnosis_group)
            label_source=f'AnnData obs[{dia_col}]'
        else:
            path=Path(args.geo_metadata)
            template=export_template(obs,ROOT/'data/raw/donor_mapping_template.csv')
            if not path.exists():raise ValueError(f'Diagnosis absent in AnnData and GEO metadata missing. Generated template: {template}')
            md=pd.read_csv(path,dtype=str)
            if args.geo_key_col not in md.columns:
                raise ValueError(f'GEO key {args.geo_key_col!r} absent. Columns={list(md.columns)}. Fill {template} with reviewed donor labels and use --donor-map.')
            lab=args.geo_diagnosis_col or guess(md.columns,['classification','diagnosis','disease_status','group','disease','condition'])
            if lab is None:raise ValueError(f'Cannot identify GEO diagnosis column. Columns={list(md.columns)}. Use --geo-diagnosis-col or reviewed --donor-map {template}')
            pairs=pd.DataFrame({'key':md[args.geo_key_col].map(key),'diagnosis':md[lab].map(diagnosis_group)})
            if pairs.groupby('key').diagnosis.nunique().gt(1).any():raise ValueError('Ambiguous GEO mapping: same key has conflicting labels; use --donor-map')
            mapping=pairs.drop_duplicates('key').set_index('key')['diagnosis']
            obs['diagnosis']=obs.donor_id.map(lambda x:mapping.get(key(x),None))
            missing=obs.loc[obs.diagnosis.isna(),'donor_id'].drop_duplicates().head(8).tolist()
            if missing:raise ValueError(f'GEO donor match failed. Sample unmatched donor IDs: {missing}. Generated {template}; create reviewed map and use --donor-map.')
            label_source=f'GEO SOFT {args.geo_key_col}→{lab}'
        per_donor=obs.groupby('donor_id').diagnosis.nunique(dropna=False)
        if per_donor.gt(1).any():raise ValueError('Conflicting diagnosis labels for one donor; cannot run donor-level DE')
        print('Diagnosis by independent donor:\n',obs[['donor_id','diagnosis']].drop_duplicates().diagnosis.value_counts().to_string())
        keep=obs.diagnosis.isin(['case','control']) & obs.cell_type.notna() & (~obs.cell_type.str.casefold().isin(['','nan','unassigned','unknown','none']))
        if not keep.any():raise ValueError('No valid schizophrenia cases and healthy controls. Inspect source diagnosis labels.')
        top=obs.loc[keep,'cell_type'].value_counts().head(args.max_celltypes).index
        keep &= obs.cell_type.isin(top)
        selected=[]
        grouping = ['donor_id', 'cell_type'] if args.sampling_mode == 'donor_celltype' else ['donor_id']
        for _,g in obs.loc[keep].groupby(grouping,sort=True):
            positions=obs.index.get_indexer(g.index)
            if (positions<0).any(): raise ValueError('Unexpected duplicate AnnData nucleus barcodes')
            picked=rng.choice(positions,size=min(args.max_nuclei_per_donor,len(positions)),replace=False)
            selected.extend(picked.tolist())
        selected=np.array(sorted(set(selected)),dtype=int)
        subset=obs.iloc[selected].copy()
        donor_groups=subset[['donor_id','diagnosis']].drop_duplicates()
        group_n=donor_groups.diagnosis.value_counts()
        if any(group_n.get(x,0)<3 for x in ['case','control']):
            raise ValueError('Fewer than three independent donors in a comparison group. Check mapping and filtering.')
        if args.count_layer=='auto':
            if 'counts' in ad.layers:source=ad.layers['counts'];layer_name='layers[counts]'
            elif ad.raw is not None:source=ad.raw.X;layer_name='raw.X'
            else:source=ad.X;layer_name='X'
        elif args.count_layer=='X':source=ad.X;layer_name='X'
        elif args.count_layer=='raw':
            if ad.raw is None:raise ValueError('No adata.raw; inspect count sources')
            source=ad.raw.X;layer_name='raw.X'
        else:
            if args.count_layer not in ad.layers:raise ValueError(f'Unknown layer {args.count_layer}; available={list(ad.layers)}')
            source=ad.layers[args.count_layer];layer_name=f'layers[{args.count_layer}]'
        var=ad.raw.var if layer_name=='raw.X' else ad.var
        # Never accept duplicated gene names; stable IDs are safer than symbols with collisions.
        symbol_col=guess(var.columns,['gene_name','gene_symbols','gene_symbol','feature_name','symbol'])
        syms=var[symbol_col].astype(str).to_numpy() if symbol_col else var.index.astype(str).to_numpy()
        ids=var.index.astype(str).to_numpy()
        names=syms if len(set(syms))==len(syms) and all(x and x!='nan' for x in syms) else ids
        if len(set(names)) != len(names):raise ValueError('Even stable gene IDs are not unique. Stop rather than collapse distinct genes.')
        pilot=selected[np.unique(np.linspace(0,len(selected)-1,min(len(selected),600),dtype=int))]
        test=source[pilot,:]
        prevalence=np.asarray((test>0).sum(axis=0)).reshape(-1)
        if args.max_genes and args.max_genes<len(names):
            genes_ix=np.sort(np.argsort(-prevalence,kind='stable')[:args.max_genes])
        else:genes_ix=np.arange(len(names))
        arrays=[]
        for start in range(0,len(selected),128):
            i=selected[start:start+128]
            x=source[i,:]
            if sparse.issparse(x):x=x[:,genes_ix].tocsr()
            else:x=sparse.csr_matrix(np.asarray(x)[:,genes_ix])
            if x.data.size and (np.any(~np.isfinite(x.data)) or np.any(x.data<0) or not np.allclose(x.data,np.rint(x.data),atol=1e-5)):
                raise ValueError(f'{layer_name} has noninteger, negative, or nonfinite entries; DESeq2 requires raw UMI counts. Stop and specify valid --count-layer.')
            arrays.append(x.astype(np.int64));print('Loaded',min(start+128,len(selected)),'/',len(selected),'nuclei')
        counts=sparse.vstack(arrays,format='csr')
        cols=['donor_id','cell_type','diagnosis']
        for standard,cands in [('age',['age']),('sex',['sex','gender']),('pmi',['pmi','postmortem_interval']),('rin',['rin','rna_integrity_number']),('lib_batch',['lib_batch','library_batch','batch'])]:
            if standard in subset:cols.append(standard)
            else:
                source_col=guess(subset.columns,cands)
                if source_col:
                    subset[standard]=subset[source_col];cols.append(standard)
        output=ROOT/'data/processed';output.mkdir(parents=True,exist_ok=True)
        subset[cols].reset_index().rename(columns={'index':'cell_id'}).to_csv(output/'cells.csv',index=False)
        pd.DataFrame({'gene':names[genes_ix]}).to_csv(output/'genes.csv',index=False)
        sparse.save_npz(output/'counts.npz',counts)
        (output/'real_data_provenance.txt').write_text('\n'.join([
            'DATA: GSE254569 - public, previously generated human brain measurements',
            f'AnnData: {args.h5ad}',f'Donor field: {donor}',f'Cell type field: {ct}',
            f'Label source: {label_source}',f'Count layer: {layer_name}',
            f'Subset nuclei: {len(selected)}',f'Genes retained: {len(genes_ix)}',
            f'Max nuclei/sampling group: {args.max_nuclei_per_donor}',f'Sampling mode: {args.sampling_mode}',f'Random seed: {args.seed}',
            'WARNING: downsampling is a research subset; not a full dataset analysis.',
            'WARNING: integer values alone do not prove UMI provenance; confirm the dataset processing record.'
        ])+'\n')
        (output/'SYNTHETIC_DATA_DO_NOT_REPORT_AS_DISCOVERY.txt').unlink(missing_ok=True)
        print('Imported real exploratory subset',counts.shape,'source=',layer_name)
    finally:
        ad.file.close()

if __name__=='__main__':main()
