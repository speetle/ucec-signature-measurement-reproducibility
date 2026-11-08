# coding:utf-8
from pathlib import Path
import numpy as np,pandas as pd,json
root=Path(__file__).resolve().parent;out=root/'output';out.mkdir(exist_ok=True);work=root/'input';base=root/'input';src=base;m=pd.read_csv(work/'复核成员清单.csv');tmpl=m.gene.tolist();common=m[m.common].gene.tolist();f=pd.read_csv(base/'discovery_matching_features.csv',index_col=0);F=(f-f.mean())/f.std().replace(0,1);mods=pd.read_csv(base/'discovery_coexpression_modules.csv',index_col=0)['module'];R1=pd.read_csv(work/'发现队列 复核RNA.csv.gz',index_col=0);R2=pd.read_csv(work/'独立患者队列 复核RNA.csv.gz',index_col=0);P1=pd.read_csv(src/'discovery/protein.csv.gz',index_col=0);P2=pd.read_csv(src/'external/protein.csv.gz',index_col=0);P1=P1.loc[:,R1.columns];P2=P2.loc[:,R2.columns];eligible=R1.index.intersection(R2.index).intersection(F.index).intersection(mods.index).difference(tmpl);protein_ok=eligible.intersection(P1.index).intersection(P2.index);protein_ok=protein_ok[(P1.loc[protein_ok].notna().mean(axis=1)>=.8)&(P2.loc[protein_ok].notna().mean(axis=1)>=.8)&(P1.loc[protein_ok].std(axis=1)>0)&(P2.loc[protein_ok].std(axis=1)>0)];obs=set(common);neighbor={}
for mode in ['feature','module']:
 neighbor[mode]={}
 for g in tmpl:
  pool=protein_ok if g in obs else eligible[(f.loc[eligible,'protein_detection']>=.8)==(f.loc[g,'protein_detection']>=.8)]
  if mode=='module':pool=pool[mods.loc[pool].values==mods[g]]
  dist=np.square(F.loc[pool].to_numpy()-F.loc[g].to_numpy()).mean(axis=1);ix=np.argsort(dist)[:min(100,len(pool))]
  if len(ix)<10:raise RuntimeError(f'{mode} {g} only {len(ix)} candidates')
  neighbor[mode][g]=(pool[ix].tolist(),dist[ix])
print('candidate pools valid; common-protein pool',len(protein_ok),flush=True)
rank_gene=R1.rank(axis=1).to_numpy(float);rank_gene-=rank_gene.mean(axis=1,keepdims=True);den=np.linalg.norm(rank_gene,axis=1);den[den==0]=1;rank_gene/=den[:,None];gi={g:i for i,g in enumerate(R1.index)};real=rank_gene[[gi[g] for g in common]];C=real@real.T;tri=np.triu_indices(len(common),1);rng=np.random.default_rng(2026110803);panels=[]
for mode,n in [('feature',1000),('module',3000)]:
 for it in range(n):
  used=set();mapping={}
  for g in rng.permutation(tmpl):
   cand,dist=neighbor[mode][g];weights=np.exp(-(dist-dist.min())/max(float(np.median(dist)),.05));weights=np.array([w if x not in used else 0 for x,w in zip(cand,weights)]);mapping[g]=rng.choice(cand,p=weights/weights.sum());used.add(mapping[g])
  names=[mapping[g] for g in common];Q=rank_gene[[gi[g] for g in names]];err=float(np.sqrt(np.mean(((Q@Q.T)[tri]-C[tri])**2)));dif=F.loc[names].to_numpy()-F.loc[common].to_numpy(); panels.append({'mode':mode,'num':it,'mapping':mapping,'covariance_rmse':err,'expression_masd':float(np.abs(dif[:,list(F.columns).index('mean')]).mean()),'detection_masd':float(np.abs(dif[:,list(F.columns).index('protein_detection')]).mean())})
 print('generated',mode,n,flush=True)
module=[x for x in panels if x['mode']=='module'];chosen=set(x['num'] for x in sorted(module,key=lambda x:x['covariance_rmse'])[:1000]); selected=[(x,'feature_only') for x in panels if x['mode']=='feature']+[(x,'module_matched') for x in module[:1000]]+[(x,'covariance_constrained') for x in module if x['num'] in chosen]
refs=[];slots=[];for_index={g:i for i,g in enumerate(tmpl)}; rankm={'discovery':R1.rank(axis=0,pct=True),'external':R2.rank(axis=0,pct=True)}; zmat={'discovery':P1.sub(P1.mean(axis=1),axis=0).div(P1.std(axis=1,ddof=1),axis=0),'external':P2.sub(P2.mean(axis=1),axis=0).div(P2.std(axis=1,ddof=1),axis=0)}
for x,ens in selected:
 for g in tmpl:slots.append({'ensemble':ens,'panel':x['num'],'template':g,'gene':x['mapping'][g],'common':g in obs})
 for co,R,P in [('discovery',R1,P1),('external',R2,P2)]:
  RR=rankm[co];Z=zmat[co];mp=x['mapping'];gu=[mp[g] for g in m.loc[m.direction.eq('up'),'gene']];gd=[mp[g] for g in m.loc[m.direction.eq('down'),'gene']];su=[mp[g] for g in m.loc[m.direction.eq('up') & m.common,'gene']];sd=[mp[g] for g in m.loc[m.direction.eq('down') & m.common,'gene']];a=RR.loc[gu].mean()-RR.loc[gd].mean();s=RR.loc[su].mean()-RR.loc[sd].mean();pro=Z.loc[su].mean()-Z.loc[sd].mean();delta=s.round(12).rank().corr(pro.rank())-a.round(12).rank().corr(pro.rank());refs.append({'cohort':co,'ensemble':ens,'panel':x['num'],'delta':delta,'covariance_rmse':x['covariance_rmse'],'expression_masd':x['expression_masd'],'detection_masd':x['detection_masd']})
pd.DataFrame(refs).to_csv(out/'共同背景新参照结果.csv',index=False);pd.DataFrame(slots).to_csv(out/'共同背景新参照成员.csv.gz',index=False,compression='gzip');print('done',len(refs),flush=True)
