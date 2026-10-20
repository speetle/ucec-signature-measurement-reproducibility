"""Fetch versioned source files directly from LinkedOmics, then reconstruct the 138-case inputs.
No third-party raw or processed LinkedOmics matrices need be redistributed by this code archive.
"""
from pathlib import Path
import argparse,hashlib,json,urllib.request,time,re
import numpy as np,pandas as pd
from excel_reader import read_xlsx
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--source-dir',type=Path);p.add_argument('--output-dir',type=Path);p.add_argument('--allow-existing-source',action='store_true');a=p.parse_args();source=a.source_dir or R/'input/external_source';out=a.output_dir or R/'input/external';source.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
record=json.loads((R/'external_sources.json').read_text())
for item in record:
 target=source/item['file']
 if not target.exists():
  tmp=target.with_suffix(target.suffix+'.part')
  for retry in range(3):
   try:
    with urllib.request.urlopen(urllib.request.Request(item['url'],headers={'User-Agent':'UCEC-reproduction/1.0'}),timeout=60) as response,tmp.open('wb') as f:
     while True:
      block=response.read(1024*1024)
      if not block:break
      f.write(block)
    tmp.replace(target);break
   except Exception:
    if retry==2:raise
    time.sleep(2)
 if hashlib.sha256(target.read_bytes()).hexdigest()!=item['sha256']:raise RuntimeError('Source version changed: '+item['file']+'; do not substitute without audit')
r=pd.read_csv(source/'RNA.cct',sep='\t',index_col=0);p=pd.read_csv(source/'protein.cct',sep='\t',index_col=0);m=read_xlsx(source/'metadata.xlsx');m['Plex']=pd.to_numeric(m['Plex'],errors='coerce').astype(float)  # preserve the original categorical label encoding, e.g. 16.0
h=pd.read_csv(R/'input/hgnc.txt',sep='\t',low_memory=False);h=h[h.status.eq('Approved')];approved=set(h.symbol);aliases={}
for _,row in h.iterrows():
 for col in ['alias_symbol','prev_symbol']:
  if pd.notna(row[col]):
   for x in str(row[col]).split('|'):
    if x not in approved:aliases.setdefault(x,set()).add(row.symbol)
alias={g:next(iter(v)) for g,v in aliases.items() if len(v)==1};coding=set(h.loc[h.locus_group.eq('protein-coding gene'),'symbol'])
def canonical_matrix(x):
 x.index=[alias.get(str(g),str(g)) for g in x.index];x=x.groupby(level=0).median();return x.loc[x.index.isin(coding)]
r=canonical_matrix(r);p=canonical_matrix(p);assert not m.Idx.duplicated().any();m=m.set_index('Idx');eligible=m[m.Group.eq('Tumor')&m.Case_excluded.eq('No')&m.discovery_study.eq('No')];ids=sorted(set(eligible.index)&set(r.columns)&set(p.columns));assert len(ids)==138
old=pd.read_csv(R/'input/discovery/CPTAC_model.csv',index_col=0).index;canon=lambda x:re.sub(r'([._-](N|T|A))$','',str(x));assert not set(map(canon,ids))&set(map(canon,old))
r[ids].to_csv(out/'RNA.csv.gz',compression={'method':'gzip','mtime':0});p[ids].to_csv(out/'protein.csv.gz',compression={'method':'gzip','mtime':0});eligible.loc[ids].to_csv(out/'metadata.csv')
print('Rebuilt external inputs from verified public source files',len(ids),str(out))
