from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.stats import rankdata,beta
ROOT=Path(__file__).resolve().parents[1]
SNAP=ROOT/"input";I=SNAP/"discovery";E=SNAP/"external";O=ROOT/"results/controls";O.mkdir(parents=True,exist_ok=True)
V=ROOT/"results"
plan=json.loads((I/'plan.json').read_text());panel=plan['formal_source_response_panel'];u=panel['up'];d=panel['down'];legacy=plan['legacy_signed_panel']['up']+plan['legacy_signed_panel']['down']
D=pd.read_csv(I/'CPTAC_model.csv',index_col=0);EM=pd.read_csv(E/'metadata.csv',index_col=0)
R=pd.read_csv(I/'RNA.csv.gz',index_col=0)[D.index];P=pd.read_csv(I/'protein.csv.gz',index_col=0)[D.index];ER=pd.read_csv(E/'RNA.csv.gz',index_col=0)[EM.index];EP=pd.read_csv(E/'protein.csv.gz',index_col=0)[EM.index]
obs=[g for g in u+d if g in P.index and P.loc[g].notna().mean()>=.8];uo=[g for g in u if g in obs];do=[g for g in d if g in obs]
cohorts={'discovery95':(R,P,D),'external138':(ER,EP,EM)}
gmt=(SNAP/'hallmark.gmt') if SNAP.exists() else W/'UCEC_证据补强_v8_20261018/benchmark/h.all.v2025.1.Hs.symbols.gmt';sets={a[0]:a[2:] for a in [l.split('\t') for l in gmt.read_text().splitlines()]}
def rho(a,b):
 a=rankdata(np.round(np.asarray(a,dtype=float),12));b=rankdata(np.round(np.asarray(b,dtype=float),12));a=a-a.mean();b=b-b.mean();return float(a@b/np.sqrt((a@a)*(b@b)))
def signed(T,a,b):return T.loc[a].mean().values-T.loc[b].mean().values
def cp(k,n):return (float(beta.ppf(.025,k,n-k+1)) if k else 0.,float(beta.ppf(.975,k+1,n-k)) if k<n else 1.)
def bh(p):
 p=np.asarray(p);q=np.full(len(p),np.nan);ok=np.flatnonzero(np.isfinite(p));i=ok[np.argsort(p[ok])];q[i]=np.minimum(1,np.minimum.accumulate((p[i]*len(i)/np.arange(1,len(i)+1))[::-1])[::-1]);return q
