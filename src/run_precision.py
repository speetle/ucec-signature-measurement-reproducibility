from pathlib import Path
import json
import numpy as np,pandas as pd,patsy
from scipy.stats import rankdata
r=Path(__file__).resolve().parents[1];old=r/'input/discovery';out=r/'results/precision';out.mkdir(parents=True,exist_ok=True);plan=json.loads((old/'plan.json').read_text());panel=plan['formal_source_response_panel'];new=json.loads((r/'input/adjudicated_panel.json').read_text());panels={'historical48':panel,'adjudicated51':{'up':new['primary_up_symbols'],'down':new['primary_down_symbols']}}
a0=pd.read_csv(old/'RNA.csv.gz',index_col=0);b0=pd.read_csv(old/'protein.csv.gz',index_col=0);a1=pd.read_csv(r/'input/external/RNA.csv.gz',index_col=0);b1=pd.read_csv(r/'input/external/protein.csv.gz',index_col=0);background=sorted(set(a0.index)&set(a1.index));fixed=[g for g in panel['up']+panel['down'] if g in b0.index and b0.loc[g].notna().mean()>=.8];rows=[]
def cp(yy,x):
 if x is not None:yy=yy-x@np.linalg.lstsq(x,yy,rcond=None)[0]
 z=rankdata(yy,axis=0);z-=z.mean(axis=0);z/=np.sqrt((z*z).sum(axis=0));return z[:,:3].T@z[:,3]
for cohort,a,b in [('discovery95',a0,b0),('external138',a1,b1)]:
 df=pd.read_csv(old/'CPTAC_model.csv',index_col=0) if cohort=='discovery95' else pd.read_csv(r/'input/external/metadata.csv',index_col=0)
 a=a[df.index];b=b[df.index]
 if cohort=='external138':
  df['age10']=pd.to_numeric(df.age,errors='coerce')/10;df['grade']=df.Histologic_grade;df['histology']=df.Histologic_Type;df['subtype']=df.Genomic_subtype;df['batch']=df.Batch
 formulas={'unadjusted':'1','clinical_subtype_batch':'age10+C(grade)+C(histology)+C(subtype)+C(batch)'}
 if cohort=='discovery95':formulas['fully_adjusted_original']='age10+C(grade)+C(histology)+C(subtype)+C(batch)+immune_logratio+stroma_logratio+proliferation+inflammation+interferon+hypoxia+oxidative_stress+dissociation'
 else:formulas['xCell_sensitivity']='age10+C(grade)+C(histology)+C(subtype)+C(batch)+xCell_immune_score+xCell_stroma_score'
 for pname,p in panels.items():
  for bg in ['native','shared_background']:
   aa=a if bg=='native' else a.loc[background];rank=aa.rank(pct=True);rz=a.loc[fixed].sub(a.loc[fixed].mean(axis=1),axis=0).div(a.loc[fixed].std(axis=1),axis=0);pz=b.loc[fixed].sub(b.loc[fixed].mean(axis=1),axis=0).div(b.loc[fixed].std(axis=1),axis=0);up=[g for g in fixed if g in p['up']];dn=[g for g in fixed if g in p['down']];score=lambda mat,u,d:mat.loc[u].mean()-mat.loc[d].mean();y=pd.DataFrame({'complete':score(rank,p['up'],p['down']),'common':score(rank,up,dn),'z':score(rz,up,dn),'protein':score(pz,up,dn)})
   for model,formula in formulas.items():
    if bg=='shared_background' and model not in ['unadjusted','clinical_subtype_batch']:continue
    design=patsy.dmatrix(formula,df,return_type='dataframe');ids=design.index;X=None if formula=='1' else design.values;Y=y.loc[ids].values;point=cp(Y,X);rng=np.random.default_rng(2026101801);B=20000 if bg=='native' else 4000;draw=[];bad=0
    for k in range(B):
     take=rng.integers(0,len(ids),len(ids));xx=None if X is None else X[take]
     if xx is not None and np.linalg.matrix_rank(xx)<xx.shape[1]:bad+=1
     c=cp(Y[take],xx);draw.append([c[1]-c[0],c[2]-c[1]])
    lo,hi=np.quantile(draw,[.025,.975],axis=0);row={'cohort':cohort,'panel':pname,'background':bg,'background_genes':len(aa),'model':model,'n':len(ids),'complete_rank_rho':point[0],'common_rank_rho':point[1],'common_z_rho':point[2],'membership_delta':point[1]-point[0],'membership_low':lo[0],'membership_high':hi[0],'scale_delta':point[2]-point[1],'scale_low':lo[1],'scale_high':hi[1],'requested':B,'successful':len(draw),'rank_deficient_draws_retained':bad,'seed':2026101801};rows.append(row);print(row,flush=True);pd.DataFrame(rows).to_csv(out/'precision_background_results.csv',index=False)
