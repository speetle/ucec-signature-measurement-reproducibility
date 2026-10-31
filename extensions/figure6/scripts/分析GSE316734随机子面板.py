from pathlib import Path
import runpy, numpy as np, pandas as pd, json
from scipy.stats import rankdata, spearmanr
P=Path(__file__).resolve().parents[1]
e=runpy.run_path(str(P/'scripts'/'分析GSE316734.py'))
m=e['mapping']; rank=e['rank']; scores=e['out']; common=m[m.member_common29.astype(str)=='True']
up=m.loc[m.direction=='up','gene'].to_list(); down=m.loc[m.direction=='down','gene'].to_list()
nu=(common.direction=='up').sum(); nd=(common.direction=='down').sum()
assert nu+nd==29 and len(up)+len(down)==48 and set(up+down)<=set(rank.index)
rng=np.random.default_rng(20261029); full=scores.source48.to_numpy(); pairs=np.triu_indices(6,1)
rows=[]
for i in range(5000):
 u=rng.choice(up,size=nu,replace=False); d=rng.choice(down,size=nd,replace=False)
 s=rank.loc[u].mean().to_numpy()-rank.loc[d].mean().to_numpy()
 rho=float(spearmanr(full,s).statistic); inv=int(np.sum(np.sign(full[pairs[0]]-full[pairs[1]])*np.sign(s[pairs[0]]-s[pairs[1]])<0))
 rows.append((i+1,rho,inv))
a=pd.DataFrame(rows,columns=['draw','spearman','discordant_pairs']);a.to_csv(P/'results'/'GSE316734随机子面板5000次.csv',index=False,encoding='utf-8-sig')
summary={'draws':len(a),'seed':20261029,'up_genes':int(nu),'down_genes':int(nd),'original_rho':float(spearmanr(full,scores.common29).statistic),'proportion_perfect_order':float((a.discordant_pairs==0).mean()),'rho_median':float(a.spearman.median()),'rho_2_5_percentile':float(a.spearman.quantile(.025)),'rho_97_5_percentile':float(a.spearman.quantile(.975)),'interpretation':'post hoc descriptive background; nested panels and six samples; not a confirmatory P value'}
(P/'results'/'GSE316734随机子面板摘要.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False));print(json.dumps(summary,indent=2,ensure_ascii=False))
