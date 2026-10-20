from controls_common import *
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits
threadpool_limits(2)
rng=np.random.default_rng(2026101902);f=pd.read_csv(V/'null/discovery_matching_features.csv',index_col=0);F=(f-f.mean())/f.std().replace(0,1);templates=u+d;eligible=F.index.difference(templates).intersection(ER.index)
# Joint presence only is used in external; no external expression features or correlations enter selection.
ep=EP.index[EP.notna().any(axis=1)];avail={True:eligible[(f.loc[eligible,'protein_detection']>=.8)&eligible.isin(ep)],False:eligible[f.loc[eligible,'protein_detection']<.8]}
genes=F.index.intersection(R.index);rr=R.loc[genes].rank(axis=1).values;rr-=rr.mean(axis=1,keepdims=True);rr/=np.sqrt((rr*rr).sum(axis=1,keepdims=True));embed=PCA(n_components=12,random_state=20261019).fit_transform(rr);mods=KMeans(n_clusters=12,random_state=20261019,n_init=10).fit_predict(embed);mod=pd.Series(mods,index=genes);pd.DataFrame({'module':mod}).to_csv(O/'discovery_coexpression_modules.csv');gn={g:i for i,g in enumerate(genes)};tar=rr[[gn[g] for g in templates]];targetcov=tar@tar.T;tri=np.triu_indices(48,1)
# Candidate pools preserve sign slot and discovery protein-detection class; module count is also exact.
neighbor={}
for scheme in ['feature','module']:
 neighbor[scheme]={}
 for g in templates:
  pp=avail[g in obs]
  if scheme=='module':pp=pp[mod.loc[pp].values==mod[g]]
  dd=np.square(F.loc[pp].values-F.loc[g].values).mean(axis=1);ix=np.argsort(dd)[:100];assert len(ix)>=10,(scheme,g,len(ix));neighbor[scheme][g]=(pp[ix].tolist(),dd[ix])
records=[];maps=[]
for scheme,n in [('feature',1000),('module',5000)]:
 for it in range(n):
  mapping={};used=set()
  for g in rng.permutation(templates):
   cand,dist=neighbor[scheme][g];w=np.exp(-(dist-dist.min())/max(np.median(dist),.05));w=np.array([ww if gg not in used else 0 for gg,ww in zip(cand,w)]);assert w.sum()>0;gg=rng.choice(cand,p=w/w.sum());mapping[g]=gg;used.add(gg)
  mg=[mapping[g] for g in templates];x=rr[[gn[g] for g in mg]];cv=x@x.T;err=float(np.sqrt(np.mean((cv[tri]-targetcov[tri])**2)));records.append(dict(generator=scheme,panel=it,covariance_RMSE=err,mean_absolute_rho=np.abs(cv[tri]).mean(),target_mean_absolute_rho=np.abs(targetcov[tri]).mean()));maps.append((scheme,it,mapping))
 print('generated',scheme,n,flush=True)
meta=pd.DataFrame(records);meta['ensemble']='not_selected';meta.loc[meta.generator=='feature','ensemble']='feature_only'
ix=meta[meta.generator=='module'].index[:1000];meta.loc[ix,'ensemble']='module_matched'
# An overlapping conditional ensemble uses the closest 1000/5000 covariance matrices, independent of all protein effects.
sel=set(meta[meta.generator=='module'].nsmallest(1000,'covariance_RMSE').panel)
meta['selected_covariance_ensemble']=(meta.generator=='module')&meta.panel.isin(sel);meta.to_csv(O/'coexpression_panel_selection.csv',index=False)
rankmat={co:r.rank(pct=True) for co,(r,p,m) in cohorts.items()};zmat={co:p.sub(p.mean(axis=1),axis=0).div(p.std(axis=1,ddof=1),axis=0) for co,(r,p,m) in cohorts.items()};out=[];mem=[];balance=[]
for scheme,it,mapping in maps:
 ensembles=['feature_only'] if scheme=='feature' else (['module_matched'] if it<1000 else [])+(['covariance_constrained'] if it in sel else [])
 if not ensembles:continue
 pu=[mapping[g] for g in u];pdn=[mapping[g] for g in d];ou=[mapping[g] for g in uo];od=[mapping[g] for g in do]
 vals={}
 for co in cohorts:
  q=rankmat[co];z=zmat[co];a=signed(q,pu,pdn);b=signed(q,ou,od);t=signed(z,ou,od);assert np.isfinite(a+b+t).all();vals[co]=(rho(a,t),rho(b,t))
 for ens in ensembles:
  for g in templates:mem.append(dict(ensemble=ens,panel=it,template=g,gene=mapping[g],direction='up' if g in u else 'down',observed=g in obs))
  dif=F.loc[[mapping[g] for g in templates]].values-F.loc[templates].values
  for k,c in enumerate(F.columns):balance.append(dict(ensemble=ens,panel=it,feature=c,mean_absolute_standardized_difference=np.abs(dif[:,k]).mean()))
  for co,(a,b) in vals.items():out.append(dict(ensemble=ens,panel=it,cohort=co,rho_full=a,rho_common=b,delta=b-a,covariance_RMSE=meta[(meta.generator==scheme)&(meta.panel==it)].covariance_RMSE.iloc[0]))
out=pd.DataFrame(out);out.to_csv(O/'coexpression_null_results.csv',index=False);pd.DataFrame(mem).to_csv(O/'coexpression_frozen_panels.csv.gz',index=False);pd.DataFrame(balance).to_csv(O/'coexpression_matching_balance.csv',index=False);summary=[]
for (ens,co),g in out.groupby(['ensemble','cohort']):
 s=pd.read_csv(V/f'null/{co}_scores.csv');od=rho(s.common29,s.protein29)-rho(s.full48,s.protein29);k=int((g.delta>=od).sum());lo,hi=cp(k,len(g)) if ens!='covariance_constrained' else (np.nan,np.nan);summary.append(dict(ensemble=ens,cohort=co,n=len(g),observed_delta=od,median_delta=g.delta.median(),tail_plus1=(k+1)/(len(g)+1),binomial_low=lo,binomial_high=hi,median_covariance_RMSE=g.covariance_RMSE.median()))
pd.DataFrame(summary).to_csv(O/'coexpression_null_summary.csv',index=False)
(O/'coexpression_policy.json').write_text(json.dumps({'seed':2026101902,'stage':'post-results exploratory structural reference, not calibrated confirmatory null','modules':'discovery-only Spearman profile PCA12 then KMeans12; exact module labels preserved for each signed/detection slot','covariance':'select 1000 lowest template-ordered pairwise rank correlation RMSE among 5000 module-matched panels; finite approximation, not exact covariance matching','ensembles_overlap':True,'external_used':'row existence and any protein observation only; no outcome, effect, expression feature or external covariance matching','full_covariance_exact':False,'reference_tail_CI':'Clopper-Pearson Monte Carlo interval conditional on generator; selection makes covariance ensemble dependent; not population/biological confidence interval'},indent=2))
print(pd.DataFrame(summary).to_string(index=False),flush=True)
