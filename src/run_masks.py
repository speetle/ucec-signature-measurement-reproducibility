from pathlib import Path
import numpy as np,pandas as pd,json,hashlib,sys
from scipy.stats import rankdata,kendalltau
W=Path(__file__).resolve().parents[1];O=W/'results/masks';O.mkdir(parents=True,exist_ok=True);I=W/'input/discovery'
A=pd.read_csv(I/'RNA.csv.gz',index_col=0)[pd.read_csv(I/'CPTAC_model.csv',index_col=0).index];B=pd.read_csv(W/'input/external/RNA.csv.gz',index_col=0)
source=json.loads((I/'plan.json').read_text())['formal_source_response_panel'];P=pd.read_csv(I/'protein.csv.gz',index_col=0);f29={d:[g for g in source[d] if g in P.index and P.loc[g].notna().mean()>=.8] for d in ['up','down']}
gmt=W/'input/hallmark.gmt';sets={a[0]:a[2:] for a in [l.split('\t') for l in gmt.read_text().splitlines()]}
panels=json.loads((W/'frozen/measurement_gene_lists.json').read_text());(O/'measurement_gene_lists.json').write_text(json.dumps(panels,indent=2))
inputs=[I/'RNA.csv.gz',I/'CPTAC_model.csv',I/'plan.json',I/'protein.csv.gz',W/'input/external/RNA.csv.gz',gmt,W/'frozen/measurement_gene_lists.json']
def metrics(x,y):
 x=np.round(np.asarray(x),12);y=np.round(np.asarray(y),12)  # eliminate float pseudoranks far below score resolution
 if not np.isfinite(x).all() or not np.isfinite(y).all() or np.std(x)==0 or np.std(y)==0:return dict(rho=np.nan,kendall_tau=np.nan,discordant_fraction_among_untied_pairs=np.nan,top_quartile_overlap_fraction=np.nan)
 rx=rankdata(x);ry=rankdata(y);rho=np.corrcoef(rx,ry)[0,1];tau=kendalltau(x,y).statistic
 ix,jx=np.triu_indices(len(x),1);prod=(x[ix]-x[jx])*(y[ix]-y[jx]);nz=prod!=0;dis=np.mean(prod[nz]<0) if nz.any() else np.nan
 k=int(np.ceil(.25*len(x)));over=len(set(np.argsort(x,kind='stable')[-k:])&set(np.argsort(y,kind='stable')[-k:]))/k
 return dict(rho=rho,kendall_tau=tau,discordant_fraction_among_untied_pairs=dis,top_quartile_overlap_fraction=over)
def signed(T,up,down):
 if len(up)<2 or len(down)<2:return np.full(T.shape[1],np.nan)
 return T.loc[up].mean().values-T.loc[down].mean().values
rows=[];scores=[]
for cohort,X in [('discovery95',A),('external138',B)]:
 full=X.rank(pct=True)
 for signature,branches in [('historical48',source),('fixed29',f29)]:
  ref=signed(full,branches['up'],branches['down'])
  for pn,genes in panels.items():
   bg=sorted(set(X.index)&set(genes));partial={d:[g for g in branches[d] if g in bg] for d in ['up','down']}
   # Midrank empirical CDF against measured background: (#< x + .5#=x)/N. Exact for observed members; permits explicitly unobservable counterfactual genes.
   vals=X.loc[bg].values; cdf=np.empty((len(X),X.shape[1]))
   for j in range(X.shape[1]):
    v=np.sort(vals[:,j]);x=X.iloc[:,j].values;cdf[:,j]=(np.searchsorted(v,x,'left')+np.searchsorted(v,x,'right'))/(2*len(bg))
   TF=pd.DataFrame(cdf,index=X.index,columns=X.columns)
   for member in ['full','observed']:
    br=branches if member=='full' else partial
    for background,T in [('full_transcriptome',full),('measured_panel',TF)]:
     z=signed(T,br['up'],br['down']);deployable=member=='observed' or all(g in bg for d in branches.values() for g in d)
     reason='' if np.isfinite(z).all() else 'fewer_than_2_members_in_either_direction'
     rows.append(dict(cohort=cohort,signature=signature,panel=pn,member_state=member,background=background,n=len(X.columns),background_genes=len(bg),up_observed=len(partial['up']),up_defined=len(branches['up']),down_observed=len(partial['down']),down_defined=len(branches['down']),deployable_members=deployable and np.isfinite(z).all(),counterfactual_fullRNA_required=(background=='full_transcriptome') or not deployable,status='computable' if not reason else 'NA',reason=reason,**metrics(ref,z)))
     for patient,value in zip(X.columns,z):scores.append(dict(cohort=cohort,signature=signature,panel=pn,member_state=member,background=background,patient=patient,score=value))
pd.DataFrame(rows).to_csv(O/'real_panel_factorial_metrics.csv',index=False);pd.DataFrame(scores).to_csv(O/'real_panel_patient_scores.csv.gz',index=False)
print('Real panels done',flush=True)
# Frozen common eligible genes avoid changing gene membership in the external set. This is exploratory and both existing datasets have already been viewed.
common=A.index.intersection(B.index);mu=A.loc[common].mean(axis=1);sd=A.loc[common].std(axis=1,ddof=1);eligible=common[sd.gt(0)&np.isfinite(sd)]
A0=A.loc[eligible];B0=B.loc[eligible];mu=mu.loc[eligible];sd=sd.loc[eligible]
# Source units are both published log2 normalized RNA but UQ equivalence is not established. Frozen-z is therefore an explicitly nonharmonized stress test, not calibrated transport.
methods={}
for co,X in [('discovery95',A),('external138',B)]:
 methods[co]={'percentile_rank':X.rank(pct=True).loc[eligible], 'discovery_frozen_z':X.loc[eligible].sub(mu,axis=0).div(sd,axis=0), 'top1500_rank_area':((1501-X.rank(ascending=False))/1500).clip(lower=0).loc[eligible]}
trainpars=pd.DataFrame({'discovery_mean':mu,'discovery_sd':sd});trainpars.to_csv(O/'frozen_z_parameters.csv.gz')
rng=np.random.default_rng(20261018);maskrows=[];out=[]
allsets={**{k:{'up':[g for g in v if g in eligible],'down':[]} for k,v in sets.items()},'source48':{d:[g for g in source[d] if g in eligible] for d in ['up','down']}}
for sn,branches in allsets.items():
 genes=branches['up']+branches['down'];is_signed=bool(branches['down']);mm=mu.loc[genes].values;vv=sd.loc[genes].values**2
 pools={'uniform':np.arange(len(genes)),'low_expression':np.arange(len(genes)),'high_variance':np.arange(len(genes))}
 if is_signed:pools.update({'up_only':np.arange(len(branches['up'])),'down_only':np.arange(len(branches['up']),len(genes))})
 for frac in [.1,.3,.5]:
  nr=int(round(len(genes)*frac))
  for mech,pool in pools.items():
   if nr>len(pool):
    maskrows.append(dict(signature=sn,mechanism=mech,fraction=frac,replicate=-1,status='NA',reason='requested_removal_exceeds_branch',removed=''));continue
   prob=None
   if mech in ['low_expression','high_variance']:
    t=mm if mech=='low_expression' else np.log(vv+1e-12);t=(t-t.mean())/(t.std()+1e-12);prob=np.exp(np.clip((-t if mech=='low_expression' else t),-3,3));prob/=prob.sum()
   keep=np.ones((100,len(genes)));removed=[]
   for rep in range(100):
    ix=rng.choice(pool,nr,replace=False,p=prob);keep[rep,ix]=0;removed.append([genes[i] for i in sorted(ix)]);maskrows.append(dict(signature=sn,mechanism=mech,fraction=frac,replicate=rep,status='computed',reason='',removed=';'.join(removed[-1])))
   w=keep.copy()
   if is_signed:
    nu=len(branches['up']);ku=keep[:,:nu].sum(axis=1);kd=keep[:,nu:].sum(axis=1);w[:,:nu]/=ku[:,None];w[:,nu:]/=-kd[:,None];bad=(ku<2)|(kd<2);w[bad]=np.nan
   else:w/=w.sum(axis=1)[:,None]
   for co in methods:
    for method,T in methods[co].items():
     tg=T.loc[genes].values;full=tg[:len(branches['up'])].mean(axis=0)-(tg[len(branches['up']):].mean(axis=0) if is_signed else 0);part=w@tg
     for rep in range(100):out.append(dict(cohort=co,signature=sn,method=method,mechanism=mech,fraction=frac,replicate=rep,members=len(genes),removed=nr,status='computed' if np.isfinite(part[rep]).all() else 'NA',**metrics(full,part[rep])))
 print(sn,len(out),flush=True)
D=pd.DataFrame(out);D.to_csv(O/'frozen_transfer_scenarios.csv.gz',index=False);pd.DataFrame(maskrows).to_csv(O/'frozen_masks.csv.gz',index=False)
keys=['cohort','method','mechanism','fraction'];D.groupby(keys).agg(rows=('rho','size'),finite=('rho','count'),median_rho=('rho','median'),q10_rho=('rho',lambda x:x.quantile(.1)),q90_rho=('rho',lambda x:x.quantile(.9)),median_discordance=('discordant_fraction_among_untied_pairs','median'),median_top_overlap=('top_quartile_overlap_fraction','median')).reset_index().to_csv(O/'frozen_transfer_summary.csv',index=False)
res=[]
for keys0,g in D.groupby(['method','mechanism','fraction']):
 Q=g.pivot(index=['signature','replicate'],columns='cohort',values='rho').dropna()
 for threshold in [.8,.9,.95]:
  y=Q.external138.values<threshold;pred=Q.discovery95.values<threshold;tp=sum(y&pred);fp=sum(~y&pred);fn=sum(y&~pred);tn=sum(~y&~pred)
  res.append(dict(zip(['method','mechanism','fraction'],keys0),threshold=threshold,n=len(Q),tp=tp,fp=fp,fn=fn,tn=tn,sensitivity=tp/(tp+fn) if tp+fn else np.nan,specificity=tn/(tn+fp) if tn+fp else np.nan,requires_full_reference=True))
pd.DataFrame(res).to_csv(O/'illustrative_threshold_sensitivity.csv',index=False)
(O/'input_sha256.json').write_text(json.dumps({str(p.relative_to(W)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2))
(O/'interpretation.json').write_text(json.dumps({'scenario_rows':len(D),'real_factorial_rows':len(rows),'frozen_definition':'All discovery means, SDs, selection weights, mask gene identities are reused unchanged for external138. No external fitting. All95 discovery and all138 evaluation; exploratory viewed datasets, no prospective test.','source_scale_caveat':'Both are published processed RNA values, but exact UQ/normalization calibration equivalence has not been established; frozen_z is nonharmonized stress sensitivity, not calibrated absolute transport. Rank methods primary.','methods':'Explicit formulas, not named software implementations. Each compares its own full reference. Unsigned Hallmarks positive mean; signed source up mean minus down mean.','NA_policy':'At least2 genes per directional branch; branch removal>size NA; missing masks never altered to force output.','dependence':'Masks, signatures and pathways dependent; no scenario P values. Thresholds illustrative and require complete reference.','RPPA':'Gene background projected on RNA, not actual RNA measurement by RPPA.','full_by_panel':'Unmeasured genes use full RNA only for explicit empirical-CDF counterfactual, never deployable.','software':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__}},indent=2))
print('DONE',len(D),flush=True)

D['signature_family']=np.where(D.signature=='source48','source48','Hallmark50')
D.groupby(['signature_family','cohort','method','mechanism','fraction']).agg(rows=('rho','size'),finite=('rho','count'),median_rho=('rho','median'),q10_rho=('rho',lambda x:x.quantile(.1)),q90_rho=('rho',lambda x:x.quantile(.9)),median_discordance=('discordant_fraction_among_untied_pairs','median'),median_top_overlap=('top_quartile_overlap_fraction','median')).reset_index().to_csv(O/'frozen_transfer_summary_by_family.csv',index=False)
D.groupby(['signature','cohort','method','mechanism','fraction']).agg(median_rho=('rho','median'),median_discordance=('discordant_fraction_among_untied_pairs','median'),median_top_overlap=('top_quartile_overlap_fraction','median')).reset_index().to_csv(O/'signature_specific_summary.csv',index=False)
