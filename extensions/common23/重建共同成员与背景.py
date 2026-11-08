# coding: utf-8
from pathlib import Path
import gzip,json,csv
import numpy as np,pandas as pd,openpyxl
import argparse
ap=argparse.ArgumentParser();ap.add_argument('project',type=Path);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--release',type=Path);ap.add_argument('--cervical',type=Path);args=ap.parse_args();root=args.project/'UCEC 生物信息学预投稿包 20261029';out=args.output;out.mkdir(parents=True,exist_ok=True);src=(args.release if args.release else args.project/'UCEC_v15_逐条实质修订_20261020/public_release')/'input';panel=json.loads((src/'discovery/plan.json').read_text())['formal_source_response_panel'];orig=pd.read_csv(src/'discovery/protein.csv.gz',index_col=0);fixed=[g for g in panel['up']+panel['down'] if g in orig.index and orig.loc[g].notna().mean()>=.8];assert len(fixed)==29
r1=pd.read_csv(src/'discovery/RNA.csv.gz',index_col=0);p1=orig;r2=pd.read_csv(src/'external/RNA.csv.gz',index_col=0);p2=pd.read_csv(src/'external/protein.csv.gz',index_col=0)
w=openpyxl.load_workbook((args.cervical if args.cervical else root/'独立队列探索 20261106/宫颈癌补充表6.xlsx'),read_only=True,data_only=True)
def load(n):
 it=w[n].iter_rows(values_only=True);head=next(it); rows=[];idx=[]
 for row in it:idx.append(row[0]);rows.append(row[1:])
 return pd.DataFrame(rows,index=idx,columns=head[1:]).apply(pd.to_numeric,errors='coerce')
r3=load('a.RNA-seq');p3=load('b.TMT_Protein'); samples=[x for x in r3.columns if x in p3.columns and str(x).startswith('T')];r3=r3[samples];p3=p3[samples];print('sizes',r1.shape,r2.shape,r3.shape,flush=True)
allsets=set(fixed)&set(p1.index)&set(p2.index)&set(p3.index)&set(r1.index)&set(r2.index)&set(r3.index);common=[g for g in fixed if g in allsets and all(p.loc[g].notna().mean()>=.8 for p in (p1,p2,p3))];up=[g for g in common if g in panel['up']];down=[g for g in common if g in panel['down']];print('common',len(common),len(up),len(down),common,flush=True)
if len(up)<2 or len(down)<2: raise ValueError('insufficient members in signed branches')
bg=set(r1.index)&set(r2.index)&set(r3.index);bg=[g for g in r1.index if g in bg and r1.loc[g].notna().all() and r2.loc[g].notna().all() and r3.loc[g].notna().all()];assert set(panel['up']+panel['down']).issubset(bg);print('background',len(bg),flush=True)

pd.DataFrame({'gene':panel['up']+panel['down'],'direction':['up']*len(panel['up'])+['down']*len(panel['down']),'common':[x in common for x in panel['up']+panel['down']]}).to_csv(out/'复核成员清单.csv',index=False)
for label,r,pro in [('发现队列',r1,p1),('独立患者队列',r2,p2)]:
 samples=[x for x in r.columns if x in pro.columns];r.loc[bg,samples].to_csv(out/(label+' 复核RNA.csv.gz'),compression='gzip')
print('上游输入重建完成',flush=True)

for label,r,pro in [('发现队列',r1,p1),('独立患者队列',r2,p2),('宫颈癌队列',r3,p3)]:
 names=[c for c in r.columns if c in pro.columns]; R=r.loc[bg,names].rank(axis=0,pct=True,method='average'); pp=pro.loc[common,names]; Z=pp.sub(pp.mean(axis=1),axis=0).div(pp.std(axis=1,ddof=1),axis=0); full=R.loc[panel['up']].mean()-R.loc[panel['down']].mean();sub=R.loc[up].mean()-R.loc[down].mean();protein=Z.loc[up].mean()-Z.loc[down].mean();pd.DataFrame({'sample':names,'full48':full.reindex(names).values,'common':sub.reindex(names).values,'protein':protein.reindex(names).values,'protein_observed':pp.notna().sum(axis=0).values}).to_csv(out/(label+' 统一子集评分.csv'),index=False)
