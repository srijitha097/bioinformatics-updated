#!/usr/bin/env python3
"""Get sample/donor metadata from GEO SOFT, independently of the large AnnData."""
import gzip,re
from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[1]; path=root/'data/raw/GSE254569_family.soft.gz'
if not path.exists(): raise SystemExit(f'Missing {path}. Run scripts/00_download_real.sh or download GEO SOFT.')
rows=[]; current=None
with gzip.open(path,'rt',errors='replace') as file:
 for line in file:
    line=line.strip()
    if line.startswith('^SAMPLE = '):
        if current: rows.append(current)
        current={'gsm':line.split('=',1)[1].strip()}
    elif current is not None and line.startswith('!Sample_title = '): current['sample_title']=line.split('=',1)[1].strip()
    elif current is not None and line.startswith('!Sample_characteristics_ch1 = '):
        content=line.split('=',1)[1].strip()
        if ':' in content:
            k,v=content.split(':',1);current[re.sub('[^a-zA-Z0-9]+','_',k.strip().lower())]=v.strip()
if current:rows.append(current)
df=pd.DataFrame(rows); out=root/'data/raw/geo_donor_metadata.csv';df.to_csv(out,index=False)
print('Wrote',out,'\nColumns:',df.columns.tolist(),'\nRows:',len(df))
