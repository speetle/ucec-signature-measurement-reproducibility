from pathlib import Path
import gzip, csv, io, json, hashlib, sys
import numpy as np,pandas as pd
from scipy.stats import spearmanr,kendalltau
W=Path(__file__).resolve().parents[1];R=W/'input/arrays';V=W/'input';out=W/'results/arrays';out.mkdir(parents=True,exist_ok=True)
plan=json.loads((V/'discovery/plan.json').read_text())['formal_source_response_panel']
ort=json.loads((V/'adjudicated_panel.json').read_text())
panels={'historical48':plan,'adjudicated51':{'up':ort['primary_up_symbols'],'down':ort['primary_down_symbols']}}
P=pd.read_csv(V/'discovery/protein.csv.gz',index_col=0)
fixed29=[g for g in plan['up']+plan['down'] if g in P.index and P.loc[g].notna().mean()>=.8]
panels['fixed_discovery29']={d:[g for g in fixed29 if g in plan[d]] for d in ['up','down']}
rows=[];coverage=[];rng=np.random.default_rng(20261018)
for name,gpl,ann in [('GSE17025','GPL570','GPL570.cached_complete.annot.gz'),('GSE63678','GPL571','GPL571.annot.gz')]:
 G=pd.read_csv(R/f'{name}_gene_expression.csv.gz',index_col=0);D=pd.read_csv(R/f'{name}_sample_metadata.csv',index_col=0)
 G.to_csv(out/f'{name}_gene_expression.csv.gz')
 if name=='GSE17025':tum=D.index[D['histology'].isin(['endometrioid','serous','papillary serous'])]
 else:tum=D.index[(D['tissue']=='endometrium') & (D['disease state']=='carcinoma')]
 for group,ids in [('tumor_primary',tum),('all_endometrium_sensitivity',D.index if name=='GSE17025' else D.index[D.tissue=='endometrium'])]:
  X=G[ids]; rank=X.rank(axis=0,pct=True)
  z=X.sub(X.mean(axis=1),axis=0).div(X.std(axis=1,ddof=1).replace(0,np.nan),axis=0)
  S={}
  for pn,pa in panels.items():
   up=[g for g in pa['up'] if g in G.index];down=[g for g in pa['down'] if g in G.index]
   for direction in ['up','down']:
    for g in pa[direction]:coverage.append({'cohort':name,'group':group,'panel':pn,'gene':g,'direction':direction,'present':g in G.index})
   for method,T in [('within_sample_gene_rank',rank),('cohort_gene_z',z)]:
    S[(pn,method)]=T.loc[up].mean()-T.loc[down].mean()
  scores=pd.DataFrame({f'{pn}_{me}':ss for (pn,me),ss in S.items()})
  scores.to_csv(out/f'{name}_{group}_scores.csv',index_label='GSM')
  for method in ['within_sample_gene_rank','cohort_gene_z']:
   for ref,target in [('historical48','fixed_discovery29'),('historical48','adjudicated51')]:
    x,y=S[(ref,method)].values,S[(target,method)].values;n=len(x)
    rho=spearmanr(x,y).statistic;tau=kendalltau(x,y).statistic
    boots=[]
    for _ in range(2000):
     ix=rng.integers(0,n,n);boots.append(spearmanr(x[ix],y[ix]).statistic)
    lo,hi=np.nanquantile(boots,[.025,.975]);k=max(1,int(np.ceil(n*.25)));u=set(np.argsort(x)[-k:]);v=set(np.argsort(y)[-k:])
    rows.append({'cohort':name,'platform':gpl,'group':group,'n':n,'method':method,'reference_panel':ref,'target_panel':target,'reference_present':sum(g in X.index for d in panels[ref].values() for g in d),'target_present':sum(g in X.index for d in panels[target].values() for g in d),'spearman_rho':rho,'rho_bootstrap_low':lo,'rho_bootstrap_high':hi,'bootstrap_draws':2000,'bootstrap_finite_draws':int(np.isfinite(boots).sum()),'kendall_tau':tau,'discordant_pair_fraction_no_ties':(1-tau)/2,'top_quartile_n':k,'top_quartile_overlap':len(u&v),'top_quartile_jaccard':len(u&v)/len(u|v),'role':'RNA-only cross-technology cohort sensitivity; not RNA-protein validation'})
  print(name,group,len(ids),flush=True)
pd.DataFrame(rows).to_csv(out/'array_panel_sensitivity.csv',index=False)
pd.DataFrame(coverage).to_csv(out/'array_gene_coverage.csv',index=False)
(out/'analysis_contract.json').write_text(json.dumps({'date':'2026-10-18','seed':20261018,'bootstrap':2000,'input_level':'GSE17025 MAS5 normalized linear signal -> log2(max(signal,1)); GSE63678 RMA already log2 kept unchanged; no CEL reconstruction','mapping':'single-Entrez probes only, joined to frozen HGNC Entrez IDs; ambiguous multi-Entrez probes excluded; median across probes per gene','scoring':'up mean minus down mean; within-sample rank over all unambiguously HGNC-mapped genes; cohort gene z uses sample SD ddof=1','primary_population':'tumors only; pooled endometrium sensitivity reported separately','inference':'descriptive patient ordering distortion; panels overlap and correlation is expected; not a biological validation, platform-equivalence test, or RNA-protein replication','fixed29':fixed29,'environment':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__}},indent=2))
