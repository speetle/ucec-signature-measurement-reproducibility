from pathlib import Path
import sys,json,hashlib,platform,time
import numpy as np,pandas as pd
from scipy.stats import rankdata
from scipy.spatial.distance import cdist
root=Path(__file__).resolve().parents[1];out=root/'results/null';out.mkdir(parents=True,exist_ok=True);v8=root/'input';old=v8/'discovery';rng=np.random.default_rng(20261018);start=time.time()
panel=json.loads((old/'plan.json').read_text())['formal_source_response_panel'];u=panel['up'];d=panel['down'];R=pd.read_csv(old/'RNA.csv.gz',index_col=0);P=pd.read_csv(old/'protein.csv.gz',index_col=0);D=pd.read_csv(old/'CPTAC_model.csv',index_col=0);R=R[D.index];P=P[D.index];obs=[g for g in u+d if g in P.index and P.loc[g].notna().mean()>=.8];uo=[g for g in u if g in obs];do=[g for g in d if g in obs];um=[g for g in u if g not in obs];dm=[g for g in d if g not in obs];assert [len(x) for x in [u,d,uo,do,um,dm]]==[32,16,19,10,13,6]
def rho(a,b):
 a=rankdata(a);b=rankdata(b);a=a-a.mean();b=b-b.mean();return float(a@b/np.sqrt((a@a)*(b@b)))
def residual(a,x):return a-x@np.linalg.lstsq(x,a,rcond=None)[0]
def score(r,p,pu,pd_,ou,od):
 rr=r.rank(pct=True);z=p.sub(p.mean(axis=1),axis=0).div(p.std(axis=1,ddof=1),axis=0)
 return np.column_stack([rr.loc[pu].mean()-rr.loc[pd_].mean(),len(ou)/len(pu)*rr.loc[ou].mean()-len(od)/len(pd_)*rr.loc[od].mean(),rr.loc[ou].mean()-rr.loc[od].mean(),z.loc[ou].mean()-z.loc[od].mean()])
E=pd.read_csv(v8/'external/metadata.csv',index_col=0);ER=pd.read_csv(v8/'external/RNA.csv.gz',index_col=0)[E.index];EP=pd.read_csv(v8/'external/protein.csv.gz',index_col=0)[E.index];E['age10']=pd.to_numeric(E.age,errors='coerce')/10;E['batch']=E.Batch;E['subtype']=E.Genomic_subtype;E['histology']=E.Histologic_Type;E['grade']=E.Histologic_grade
cohorts={'discovery95':(R,P,D),'external138':(ER,EP,E)};rows=[];draws=[];perms=[];strata=[]
for name,(r,p,df) in cohorts.items():
 y=score(r,p,u,d,uo,do);pd.DataFrame(y,index=df.index,columns=['full48','observed_fixed_branch_weight','common29','protein29']).to_csv(out/f'{name}_scores.csv',index_label='patient')
 for model,formula in [('unadjusted','1'),('clinical_subtype_batch','age10+C(grade)+C(histology)+C(subtype)+C(batch)')]:
  X=pd.DataFrame({'Intercept':np.ones(len(df))},index=df.index) if formula=='1' else pd.concat([pd.Series(1.,index=df.index,name='Intercept'),df[['age10']],pd.get_dummies(df[['grade','histology','subtype','batch']].astype(str),drop_first=True,dtype=float)],axis=1).loc[df[['age10','grade','histology','subtype','batch']].notna().all(axis=1)];X=X.loc[:,X.ne(0).any(axis=0)];ix=df.index.get_indexer(X.index);ys=y[ix];x=X.values;good=np.isfinite(ys).all(axis=1);ys=ys[good];x=x[good];n=len(ys)
  def calc(a,b):
   z=a if model=='unadjusted' else residual(a,b);c=[rho(z[:,j],z[:,3]) for j in range(3)];return np.array(c+[c[1]-c[0],c[2]-c[1],c[2]-c[0]])
  val=calc(ys,x);bs=[];rankdef=0
  for it in range(2000):
   j=rng.integers(n,size=n);bs.append(calc(ys[j],x[j]));rankdef+=np.linalg.matrix_rank(x[j])<x.shape[1]
  bs=np.array(bs);lo,hi=np.quantile(bs,[.025,.975],axis=0)
  for k,metric in enumerate(['rho_full48','rho_fixed_weight_observed','rho_common29','deletion_step','reweight_step','total_delta']):rows.append(dict(cohort=name,model=model,n=n,metric=metric,value=val[k],low=lo[k],high=hi[k],draws=2000,rank_deficient_draws=rankdef,design_columns=x.shape[1],design_rank=int(np.linalg.matrix_rank(x))))
  for it,b in enumerate(bs):draws.append([name,model,it,*b])
 # Protein scores are standardized before permutations; move whole patient vectors/scores together. Original batches only, singleton fixed.
 batches=df['batch'].astype(str).values;groups=[np.where(batches==b)[0] for b in pd.unique(batches)]
 for g in groups:strata.append(dict(cohort=name,batch=batches[g[0]],n=len(g),exchangeable=len(g)>1,policy='permute within original batch' if len(g)>1 else 'fixed singleton; no cross-batch merge'))
 target=y[:,3];rv=[rho(y[:,j],target) for j in [0,2]]
 for it in range(1000):
  idx=np.arange(len(df))
  for g in groups:idx[g]=rng.permutation(g)
  a=rho(y[:,0],target[idx]);b=rho(y[:,2],target[idx]);perms.append(dict(cohort=name,iteration=it,rho_full48=a,rho_common29=b,total_delta=b-a,observed_delta=rv[1]-rv[0]))
 print('bootstrap+permutation',name,flush=True)
pd.DataFrame(rows).to_csv(out/'weight_contrasts_bootstrap.csv',index=False);pd.DataFrame(draws,columns=['cohort','model','iteration','rho_full48','rho_observed','rho_common29','deletion_step','reweight_step','total_delta']).to_csv(out/'bootstrap_draws.csv.gz',index=False);pd.DataFrame(perms).to_csv(out/'within_batch_permutations.csv',index=False);pd.DataFrame(strata).to_csv(out/'permutation_strata.csv',index=False)
# Matching parameters are estimated solely in discovery. Target-specific rank and background never refit on external patients.
genes=R.index[R.notna().all(axis=1)&(R.std(axis=1)>0)];det=P.notna().mean(axis=1).reindex(genes).fillna(0);features=pd.DataFrame({'mean':R.loc[genes].mean(axis=1),'log_sd':np.log(R.loc[genes].std(axis=1)),'RNA_detection':(R.loc[genes]>0).mean(axis=1),'protein_detection':det})
rr=R.loc[genes].rank(axis=1).values;rr-=rr.mean(axis=1,keepdims=True);rr/=np.sqrt((rr*rr).sum(axis=1,keepdims=True))
for c in ['proliferation','inflammation','interferon','hypoxia','oxidative_stress','dissociation']:
 a=rankdata(D[c].values);a-=a.mean();a/=np.linalg.norm(a);features[c]=rr@a
sd=features.std().replace(0,1);F=(features-features.mean())/sd;features.to_csv(out/'discovery_matching_features.csv');eligible=F.index.difference(u+d);pools={True:eligible[det.loc[eligible]>=.8],False:eligible[det.loc[eligible]<.8]};templates=u+d;neigh={}
for g in templates:
 pp=pools[g in obs];dist=np.square(F.loc[pp].values-F.loc[g].values).mean(axis=1);order=np.argsort(dist)[:100];neigh[g]=(pp[order].tolist(),dist[order])
print('candidate sizes', {str(k):len(v) for k,v in pools.items()},flush=True)
rankmat={c:r.rank(pct=True) for c,(r,p,df) in cohorts.items()};zmat={c:p.sub(p.mean(axis=1),axis=0).div(p.std(axis=1,ddof=1),axis=0) for c,(r,p,df) in cohorts.items()};null=[];members=[];balance=[]
for it in range(1000):
 mapping={};used=set()
 for g in rng.permutation(templates):
  cand,dist=neigh[g];w=np.exp(-(dist-dist.min())/max(np.median(dist),.05));w=np.array([ww if gg not in used else 0 for gg,ww in zip(cand,w)]);gg=rng.choice(cand,p=w/w.sum());mapping[g]=gg;used.add(gg);members.append(dict(panel=it,template=g,gene=gg,direction='up' if g in u else 'down',observed=g in obs))
 pu=[mapping[g] for g in u];pdn=[mapping[g] for g in d];ou=[mapping[g] for g in uo];od=[mapping[g] for g in do]
 dif=F.loc[[mapping[g] for g in templates]].values-F.loc[templates].values
 for k,c in enumerate(F.columns):balance.append(dict(panel=it,feature=c,mean_standardized_difference=dif[:,k].mean(),mean_absolute_standardized_difference=np.abs(dif[:,k]).mean()))
 for name,(r,p,df) in cohorts.items():
  missingRNA=set(pu+pdn)-set(r.index);missingP=set(ou+od)-set(p.index)
  if missingRNA or missingP:null.append(dict(cohort=name,panel=it,computable=False,missing_RNA=';'.join(sorted(missingRNA)),missing_protein=';'.join(sorted(missingP))));continue
  q=rankmat[name];z=zmat[name];a=(q.loc[pu].mean()-q.loc[pdn].mean()).values;b=(q.loc[ou].mean()-q.loc[od].mean()).values;t=(z.loc[ou].mean()-z.loc[od].mean()).values;good=np.isfinite(a+b+t);c1=rho(a[good],t[good]);c2=rho(b[good],t[good]);null.append(dict(cohort=name,panel=it,computable=True,n=int(good.sum()),rho_full48=c1,rho_common29=c2,total_delta=c2-c1))
pd.DataFrame(null).to_csv(out/'matched_panel_null.csv',index=False);pd.DataFrame(members).to_csv(out/'frozen_random_panels.csv.gz',index=False);pd.DataFrame(balance).to_csv(out/'matching_balance.csv',index=False)
summ=[]
for c in cohorts:
 obsdelta=next(z['value'] for z in rows if z['cohort']==c and z['model']=='unadjusted' and z['metric']=='total_delta')
 for typ,t in [('matched_panel',pd.DataFrame(null)),('within_batch_permutation',pd.DataFrame(perms))]:
  v=t.loc[t.cohort==c,'total_delta'].dropna().values;summ.append(dict(cohort=c,null_type=typ,observed_delta=obsdelta,computable=len(v),median=float(np.median(v)),low=float(np.quantile(v,.025)),high=float(np.quantile(v,.975)),upper_tail_fraction_plus1=float((1+(v>=obsdelta).sum())/(1+len(v)))))
pd.DataFrame(summ).to_csv(out/'null_summary.csv',index=False)
inputs=[old/'plan.json',old/'RNA.csv.gz',old/'protein.csv.gz',old/'CPTAC_model.csv',v8/'external/RNA.csv.gz',v8/'external/protein.csv.gz',v8/'external/metadata.csv']
(out/'provenance.json').write_text(json.dumps({'seed':20261018,'runtime_seconds':time.time()-start,'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'matching':'discovery-only standardized RNA mean/log SD/detection, protein detection and six ordinary-program correlations; 100 nearest candidates sampled distance-weighted without replacement within panel; enforce detection groups; full covariance not matched','external':'same frozen panel identities; no re-selection; absent required members => noncomputable','bootstrap':'paired patient resampling, nuisance coefficients refit; gene z scales fixed; rank-deficient draws minimum-norm retained','permutation':'whole protein patient vector equivalent score permutation within original batches; singleton fixed; no cross-batch merge','interpretation':'exploratory structural/comparison null, not causal decomposition or independent biological replicates'},indent=2))
print(pd.DataFrame(summ).to_string(index=False),flush=True)
