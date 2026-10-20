from controls_common import *
rows=[]
for co,(r,p,m) in cohorts.items():
 q=r.rank(pct=True);z=p.sub(p.mean(axis=1),axis=0).div(p.std(axis=1,ddof=1),axis=0);a=signed(q,u,d);b=signed(q,uo,do);g='LOXL3';missing=z.loc[g].isna().values;known=z.loc[g].dropna().values;n=missing.sum()
 for bound,lo,hi in [('observed_q05_q95',*np.quantile(known,[.05,.95])),('observed_mean_plusminus2SD',known.mean()-2*known.std(ddof=1),known.mean()+2*known.std(ddof=1))]:
  for scenario in ['all_low','all_high','aligned_common','inverse_common','aligned_full','inverse_full']:
   vals=np.full(n,lo if scenario=='all_low' else hi)
   if scenario not in ['all_low','all_high']:
    score=b if 'common' in scenario else a;rk=rankdata(score[missing]);vals=lo+(rk-1)/(n-1)*(hi-lo)
    if scenario.startswith('inverse'):vals=hi+lo-vals
   zz=z.copy();zz.loc[g,zz.columns[missing]]=vals;tt=signed(zz,uo,do);rows.append(dict(cohort=co,gene=g,missing_assigned=int(n),bound=bound,scenario=scenario,assumed_low=lo,assumed_high=hi,delta=rho(b,tt)-rho(a,tt),scope='deterministic assumed LOXL3 missing-z stress scenario; not imputed experimental observations or identification bounds'))
pd.DataFrame(rows).to_csv(O/'LOXL3_assumption_stress.csv',index=False);print(pd.DataFrame(rows).groupby('cohort').delta.agg(['min','max']))
