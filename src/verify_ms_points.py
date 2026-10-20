from common import *
from scipy.stats import spearmanr
import hashlib,sys
x=pd.read_csv(O/'hallmark_actual_correspondence.csv');errors=[]
for co,(r,p,m) in cohorts.items():
 rank=r.rank(pct=True);f=pd.read_csv(O/'hallmark_frozen_members.csv');
 for _,row in x[x.cohort.eq(co)].iterrows():
  fm=f[f['set'].eq(row['set'])];full=fm.gene.tolist();com=fm[fm.protein_common].gene.tolist();z=p.loc[com].apply(lambda v:(v-v.mean())/v.std(ddof=1),axis=1);t=z.mean();a=rank.loc[full].mean();b=rank.loc[com].mean();valid=(z.notna().sum()>=int(np.ceil(.8*len(com))))&t.notna();v=spearmanr(a[valid].round(12),t[valid].round(12)).statistic;w=spearmanr(b[valid].round(12),t[valid].round(12)).statistic;errors.append(abs(w-v-row.delta))
source=[]
for co,(r,p,m) in cohorts.items():
 q=r.rank(pct=True);z=p.loc[obs].sub(p.loc[obs].mean(axis=1),axis=0).div(p.loc[obs].std(axis=1),axis=0);a=signed(q,u,d);b=signed(q,uo,do);t=signed(z,uo,do);val=spearmanr(np.round(b,12),np.round(t,12)).statistic-spearmanr(np.round(a,12),np.round(t,12)).statistic;source.append({'cohort':co,'delta':val})
# Identical pathway scope, preserving common patient bootstrap draws.
shared=sorted(set(x[x.cohort.eq('discovery95')]['set'])&set(x[x.cohort.eq('external138')]['set']));rows=[]
for co in cohorts:
 z=x[x.cohort.eq(co)&x['set'].isin(shared)];b=pd.read_csv(O/f'{co}_hallmark_bootstraps.csv.gz')[shared];draw=b.median(axis=1);rows.append(dict(cohort=co,sets=len(shared),positive=int((z.delta>0).sum()),negative=int((z.delta<0).sum()),median_delta=z.delta.median(),low=draw.quantile(.025),high=draw.quantile(.975),interpretation='Fixed shared eligible scope; conditional protein-scale intervals, overlapping sets not independent'))
pd.DataFrame(rows).to_csv(O/'hallmark_shared_scope_summary.csv',index=False)
# Attribute every non-transportable member, no external reselection.
f=pd.read_csv(O/'hallmark_frozen_members.csv');attrs=[]
for name,g in f.groupby('set'):
 attrs.append(dict(set=name,missing_external_RNA=';'.join([v for v in g.gene if v not in ER.index]),missing_external_protein=';'.join([v for v in g[g.protein_common].gene if v not in EP.index])))
pd.DataFrame(attrs).to_csv(O/'hallmark_external_missing_members.csv',index=False)
rec={'independent_implementation':'scipy spearmanr + pandas per-row z; same input matrices, not independent experiment','Hallmark_rows_checked':len(errors),'max_absolute_delta_error':max(errors),'source_primary':source,'shared_scope_sets':len(shared),'all_passed':bool(max(errors)<1e-12)}
(O/'independent_point_validation.json').write_text(json.dumps(rec,indent=2));assert rec['all_passed'];print(json.dumps(rec,indent=2));print(pd.DataFrame(rows).to_string(index=False))
