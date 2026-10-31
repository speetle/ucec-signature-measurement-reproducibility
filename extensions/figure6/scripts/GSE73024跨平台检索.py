from pathlib import Path
import pandas as pd,numpy as np
from scipy.stats import rankdata
P=Path(__file__).resolve().parents[1];(P/'results').mkdir(exist_ok=True);z=np.load(P/'assets'/'Tahoe板2_AN3官方响应.npz');g=pd.read_csv(P/'assets'/'GSE73024_AN3CA基因表达与响应.csv',index_col=0);map=pd.read_csv(P/'assets'/'Tahoe源48及共有29基因映射.csv');valid=pd.read_csv(P/'assets'/'Tahoe板2官方响应完整性.csv');assert valid.complete.all() and len(valid)==92
assert list(z['drug'])==list(valid.drug); genes=z['genes'];index={x:i for i,x in enumerate(genes)};common=g.index.intersection(index);ss=set(map.gene);cc=set(map.loc[map.member_common29.astype(str)=='True','gene']);sets={'source48':ss,'common29':cc,'full':set(common)};minimum={'source48':36,'common29':22,'full':2000};positive={'Infigratinib','Erdafitinib','Pemigatinib','Futibatinib'};assert positive<=set(z['drug'])
rows=[]
def pearson_rank(x,y):
 x=rankdata(x);y=rankdata(y);x-=x.mean();y-=y.mean();den=np.sqrt(x@x*(y@y));return float(x@y/den) if den else np.nan
cases=['mean_drug_minus_control']+[x for x in g.columns if ' minus ' in x];assert len(cases)==5
for case in cases:
 for panel,gset in sets.items():
  gg=sorted(set(common)&gset);tidx=np.array([index[x] for x in gg]);query=g.loc[gg,case].to_numpy();
  for drugid,name in enumerate(z['drug']):
   reference=z['values'][drugid,tidx];mask=np.isfinite(query)&np.isfinite(reference);n=int(mask.sum());score=pearson_rank(query[mask],reference[mask]) if n>=minimum[panel] else np.nan
   rows.append({'case':case,'panel':panel,'drug':name,'common_features':n,'score':score,'FGFR_positive':name in positive,'full_cells':int(valid.iloc[drugid].cells_trt)})
r=pd.DataFrame(rows);output=[];summary=[]
for scenario in ['all92','ge100']:
 for (case,panel),part in r.groupby(['case','panel'],sort=False):
  part=part.copy();part['eligible']=True if scenario=='all92' else part.full_cells>=100;part['rank']=part.score.where(part.eligible).rank(ascending=False,method='average');part['scenario']=scenario;output.append(part)
  pos=part[part.FGFR_positive & part.eligible]; summary.append({'scenario':scenario,'case':case,'panel':panel,'eligible_candidates':int(part.eligible.sum()),'scorable_candidates':int(part['rank'].notna().sum()),'positive_available':len(pos),'positive_best_rank':float(pos['rank'].min()) if pos['rank'].notna().any() else np.nan,'FGFR_in_top5':int((pos['rank']<=5).sum()),'FGFR_in_top10':int((pos['rank']<=10).sum()),'positive_ranks':'; '.join(f"{a}: {b:.1f}" if np.isfinite(b) else f"{a}: NA" for a,b in zip(pos.drug,pos['rank'])),'note':'independent study, same cell line and drug class, different compounds/doses/platforms; descriptive'})
pd.concat(output).to_csv(P/'results'/'GSE73024跨平台全部候选排名.csv',index=False,encoding='utf-8-sig');pd.DataFrame(summary).to_csv(P/'results'/'GSE73024跨平台检索摘要.csv',index=False,encoding='utf-8-sig');print(pd.DataFrame(summary).query("case=='mean_drug_minus_control'").to_string(index=False))
