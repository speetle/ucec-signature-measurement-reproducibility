from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,time,sys
from scipy.stats import rankdata,spearmanr
R=Path(__file__).resolve().parents[2];P=R;O=R/'results'/'v17_review';O.mkdir(parents=True,exist_ok=True);start=time.time();used=[]
def read(p,**kw):used.append(p);return pd.read_csv(p,**kw)
def cor(x,y):
 a=rankdata(x,axis=-1);b=rankdata(y,axis=-1);a-=a.mean(axis=-1,keepdims=True);b-=b.mean(axis=-1,keepdims=True)
 return np.sum(a*b,axis=-1)/np.sqrt(np.sum(a*a,axis=-1)*np.sum(b*b,axis=-1))
rows=[]
for co in ['discovery95','external138']:
 d=read(P/f'results/null/{co}_scores.csv').sort_values('patient');n=len(d);x=np.round(d.full48.to_numpy(),12);y=np.round(d.common29.to_numpy(),12);p=d.protein29.to_numpy();assert len(np.unique(x))==n and len(np.unique(y))==n;ix=np.triu_indices(n,1);dx=np.sign(x[ix[0]]-x[ix[1]]);dy=np.sign(y[ix[0]]-y[ix[1]]);ok=(dx!=0)&(dy!=0);discord=int(np.sum(dx[ok]!=dy[ok]));k=int(np.ceil(n*.25));topx=np.argsort(-x,kind='stable')[:k];topy=np.argsort(-y,kind='stable')[:k];inter=len(set(topx)&set(topy));seed=2026102201 if co=='discovery95' else 2026102202;rng=np.random.default_rng(seed);draw=rng.integers(0,n,(2000,n));rho=cor(x,y);rr=cor(x[draw],y[draw]);lo,hi=np.quantile(rr,[.025,.975]);dr=cor(y,p)-cor(x,p)
 rows.append(dict(cohort=co,n=n,full_common_rho=float(rho),rho_low=float(lo),rho_high=float(hi),untied_pairs=int(ok.sum()),discordant_pairs=discord,discordance_percent=100*discord/int(ok.sum()),top_n=k,retained=inter,retention_percent=100*inter/k,protein_delta=float(dr),seed=seed,bootstrap=2000))
 d.assign(full_quartile=np.isin(np.arange(n),topx),common_quartile=np.isin(np.arange(n),topy)).to_csv(O/f'{co}_ordering_patients.csv',index=False)
pd.DataFrame(rows).to_csv(O/'main_ordering_tradeoff.csv',index=False)
# Same 25 sets, same 422 participants, shared draws, original frozen score scales.
s=read(P/'results/rppa/RPPA_Hallmark_patient_scores.csv.gz');s=s[s.method.eq('within_sample_gene_rank')];names=sorted(s['set'].unique());ids=sorted(s.patient.unique());assert len(names)==25 and len(ids)==422
rng=np.random.default_rng(20261018);draw=rng.integers(0,len(ids),(2000,len(ids)));delta=[];original=[]
for name in names:
 d=s[s['set'].eq(name)].set_index('patient').loc[ids];assert np.isfinite(d[['full_RNA_score','common_RNA_score','common_protein_gene_z_score']].to_numpy()).all()
 x=d.full_RNA_score.to_numpy();y=d.common_RNA_score.to_numpy();p=d.common_protein_gene_z_score.to_numpy();delta.append(cor(y[draw],p[draw])-cor(x[draw],p[draw]));original.append(cor(y,p)-cor(x,p))
a=np.array(delta).T;med=np.median(a,axis=1);assert np.isfinite(med).all();lo,hi=np.quantile(med,[.025,.975]);point=float(np.median(original));ref=read(P/'results/rppa/RPPA_Hallmark_dependent_correlations.csv').query('method=="within_sample_gene_rank"').set_index('set').loc[names]
assert np.allclose(original,ref.paired_delta.to_numpy(),atol=1e-12)
pd.DataFrame(a,columns=names).assign(median_delta=med).to_csv(O/'rppa_shared_patient_bootstrap.csv.gz',index=False)
z=dict(n=422,sets=25,point=point,low=float(lo),high=float(hi),draws=2000,finite=2000,seed=20261018,scale='conditional on original cohort gene-z scales and frozen complete-observed score vectors',unit='patient; same index resample for all overlapping sets',positive=int((np.array(original)>0).sum()),negative=int((np.array(original)<0).sum()))
(O/'rppa_aggregate.json').write_text(json.dumps(z,indent=2));print(json.dumps({'main_ordering':rows,'rppa':z},indent=2),flush=True)
used=list(dict.fromkeys(used));rec=dict(date='2026-10-22',scope='new v17 ordering and fixed-scale RPPA aggregate analyses only; not a complete upstream rerun',python=sys.version,seconds=time.time()-start,inputs={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in used},script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),outputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file()})
(O/'new_analysis_run.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2))
