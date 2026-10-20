"""Audit actual protein denominators without redefining the primary endpoint.
Complete-coverage comparisons are post-review, selected-subpopulation sensitivities.
"""
from controls_common import *
from scipy.stats import rankdata
Q=ROOT/'results/review';Q.mkdir(parents=True,exist_ok=True)
patients=[];summaries=[];complete=[]
for co,(r,p,m) in cohorts.items():
 z=p.loc[obs].sub(p.loc[obs].mean(axis=1),axis=0).div(p.loc[obs].std(axis=1),axis=0)
 batch='batch' if co=='discovery95' else 'Plex'
 counts=pd.DataFrame({'patient':m.index,'cohort':co,'batch':m[batch].astype(str).values,'up_observed':z.loc[uo].notna().sum().values,'down_observed':z.loc[do].notna().sum().values})
 counts['total_observed']=counts.up_observed+counts.down_observed;counts['complete29']=counts.total_observed.eq(29);patients.append(counts)
 summaries.append(dict(cohort=co,n=len(m),up_min=counts.up_observed.min(),up_median=counts.up_observed.median(),up_max=counts.up_observed.max(),down_min=counts.down_observed.min(),down_median=counts.down_observed.median(),down_max=counts.down_observed.max(),total_min=counts.total_observed.min(),total_median=counts.total_observed.median(),total_max=counts.total_observed.max(),complete29=int(counts.complete29.sum())))
 take=np.flatnonzero(counts.complete29);scores=pd.read_csv(ROOT/'results/primary'/f'{co}_scores.csv');Y=scores[['full48','common29','protein29']].to_numpy()[take]
 def corr_pair(A):
  a=rankdata(A,axis=-2);a-=a.mean(axis=-2,keepdims=True);a/=np.sqrt(np.sum(a*a,axis=-2,keepdims=True));cc=np.sum(a[...,:2]*a[...,2,None],axis=-2);return cc
 def rs(A):
  cc=corr_pair(A);return cc[...,1]-cc[...,0]
 row=dict(cohort=co,n=len(take),resamples=0,seed=2026102101,delta=np.nan,low=np.nan,high=np.nan,scope='post-review complete29 selected patients; original cohort protein scales fixed; not missing-at-random identification')
 if len(take)>=30:
  rng=np.random.default_rng(2026102101);jj=rng.integers(len(take),size=(20000,len(take)));draw=np.concatenate([rs(Y[t]) for t in np.array_split(jj,20)]);cc=corr_pair(Y);row.update(rho_full=float(cc[0]),rho_common=float(cc[1]),delta=float(rs(Y)),low=float(np.quantile(draw,.025)),high=float(np.quantile(draw,.975)),resamples=20000)
 else:row['scope']+='; fewer than 30, inference not computed'
 complete.append(row)
allp=pd.concat(patients,ignore_index=True);allp.to_csv(Q/'patient_protein_branch_coverage.csv',index=False);pd.DataFrame(summaries).to_csv(Q/'patient_protein_coverage_summary.csv',index=False);pd.DataFrame(complete).to_csv(Q/'complete29_patient_sensitivity.csv',index=False)
allp.groupby(['cohort','batch'],sort=True).agg(n=('patient','size'),up_min=('up_observed','min'),up_median=('up_observed','median'),up_max=('up_observed','max'),down_min=('down_observed','min'),down_median=('down_observed','median'),down_max=('down_observed','max'),total_min=('total_observed','min'),total_median=('total_observed','median'),total_max=('total_observed','max'),complete29=('complete29','sum')).reset_index().to_csv(Q/'protein_coverage_by_batch.csv',index=False)
assert len(allp)==233 and (allp.up_observed>=1).all() and (allp.down_observed>=1).all()
print(pd.DataFrame(summaries).to_string(index=False));print(pd.DataFrame(complete).to_string(index=False))
