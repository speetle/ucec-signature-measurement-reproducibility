"""Recompute displayed donor/mixture analyses from frozen gene-by-lineage counts.
Raw count-matrix QC and clustering are upstream, separately archived and not rerun here.
"""
from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
from scipy.stats import rankdata,wilcoxon
R=Path(__file__).resolve().parents[1];I=R/'input/cells';O=R/'results/cells';O.mkdir(parents=True,exist_ok=True)
plan=json.loads((I/'plan.json').read_text());meta=json.loads((I/'samples.json').read_text());hgnc=pd.read_csv(R/'input/hgnc.txt',sep='\t',low_memory=False);approved=set(hgnc.symbol);protein=set(hgnc.loc[hgnc.locus_group.eq('protein-coding gene'),'symbol']);aliases={}
for _,a in hgnc.iterrows():
 for col in ['prev_symbol','alias_symbol']:
  if pd.notna(a[col]):
   for g in str(a[col]).split('|'):aliases.setdefault(g,set()).add(a.symbol)
rename={g:next(iter(v)) for g,v in aliases.items() if len(v)==1 and g not in approved}
def read(s,correction):
 x=pd.read_csv(I/s/(correction+'_pseudobulk.csv.gz'),index_col=0);x.index=[rename.get(g,g) for g in x.index];return x.groupby(level=0,sort=False).sum()
corrected={a['sample']:read(a['sample'],'corrected') for a in meta};common=sorted(protein.intersection(set.intersection(*(set(x.index) for x in corrected.values()))))
(O/'rank_background.json').write_text(json.dumps(common));rows=[];profiles={};cellcounts=[];qcs=[];thresholds=[];evidence=[]
for a in meta:
 s=a['sample'];counts=pd.read_csv(I/s/'患者谱系细胞数.csv');counts['dataset']=a['dataset'];cellcounts.append(counts);qcs.append(pd.read_csv(I/s/'QC汇总.csv'));thresholds.append(pd.read_csv(I/s/'QC阈值.csv'));ev=pd.read_csv(I/s/'注释证据.csv');ev['sample']=s;evidence.append(ev)
 for correction in ['raw','corrected']:
  x=(read(s,'raw') if correction=='raw' else corrected[s]).loc[common];cpm=x.div(x.sum(),axis=1)*1e6;rank=cpm.rank(axis=0,pct=True)
  if correction=='raw':profiles[s]=cpm
  for lin in x.columns:
   z={**a,'lineage':lin,'n_cells':int(counts.loc[counts.lineage.eq(lin),'cells'].iloc[0]),'correction':correction}
   for panel,key in [('legacy','legacy_signed_panel'),('response','formal_source_response_panel')]:
    p=plan[key];assert set(p['up']+p['down']).issubset(rank.index);z[panel]=rank.loc[p['up'],lin].mean()-rank.loc[p['down'],lin].mean()
   rows.append(z)
d=pd.DataFrame(rows);d.to_csv(O/'lineage_pseudobulk.csv',index=False)
# Primary family: 2 prespecified scores x 5 epithelial-to-other lineage contrasts in the 5 known donors.
# Enumerating all sign assignments of ranked absolute nonzero donor differences supplies exact two-sided signed-rank p.
def exact_signed_rank(z):
 z=np.asarray(z,dtype=float);z=z[np.isfinite(z)]
 if (z==0).any():raise ValueError('Zero donor difference: exact no-zero policy requires a declared alternative')
 ranks=rankdata(np.abs(z));w=float(ranks[z>0].sum());total=float(ranks.sum());vals=np.array([np.dot(v,ranks) for v in itertools.product([0,1],repeat=len(z))]);tail=np.minimum(vals,total-vals);p=float(np.mean(tail<=min(w,total-w)+1e-12))
 return p,w,len(vals),bool(len(np.unique(np.abs(z)))<len(z))
tests=[]
for panel in ['legacy','response']:
 f=d.query('dataset=="GSE173682" and correction=="corrected" and n_cells>=30').pivot(index='sample',columns='lineage',values=panel)
 for lin in ['Myeloid','Stromal','Lymphoid','Endothelial','Perivascular']:
  if lin not in f:continue
  z=(f.Epithelial-f[lin]).dropna();p,w,k,tied=exact_signed_rank(z.values);p2=float(wilcoxon(z.values,method='exact',alternative='two-sided',zero_method='wilcox',correction=False).pvalue);assert abs(p-p2)<1e-12
  tests.append(dict(panel=panel,contrast='Epithelial minus '+lin,n=len(z),donors=';'.join(z.index),median_difference=z.median(),p=p,positive_rank_sum=w,enumerated_sign_assignments=k,zero_differences=0,tied_absolute_differences=tied,algorithm='two-sided exact signed-rank, all 2^n sign assignments',family='2 scores x 5 lineages = 10 whole-cell comparisons'))
t=pd.DataFrame(tests);assert len(t)==10
order=np.argsort(t.p.to_numpy());p=t.p.to_numpy();q=np.empty(10);q[order]=np.minimum(1,np.minimum.accumulate((p[order]*10/np.arange(1,11))[::-1])[::-1]);t['BH10']=q
p18=np.r_[p,np.ones(8)];order=np.argsort(p18);qq=np.empty(18);qq[order]=np.minimum(1,np.minimum.accumulate((p18[order]*18/np.arange(1,19))[::-1])[::-1]);t['historical_BH18_sensitivity']=qq[:10];t.to_csv(O/'exact_donor_tests.csv',index=False)
mix=[]
for a in meta:
 if a['dataset']!='GSE173682':continue
 s=a['sample'];x=profiles[s]
 for fraction in np.linspace(0,1,21):
  y=x.Epithelial*fraction+x.Stromal*(1-fraction);rank=y.rank(pct=True)
  for panel,key in [('legacy','legacy_signed_panel'),('response','formal_source_response_panel')]:
   pp=plan[key];mix.append(dict(sample=s,epithelial_RNA_fraction=fraction,panel=panel,score=rank.loc[pp['up']].mean()-rank.loc[pp['down']].mean()))
pd.DataFrame(mix).to_csv(O/'admixture_curves.csv',index=False);pd.concat(cellcounts).to_csv(O/'lineage_cell_counts.csv',index=False);pd.concat(qcs).to_csv(O/'upstream_QC_record.csv',index=False);pd.concat(thresholds).to_csv(O/'upstream_QC_thresholds.csv',index=False);pd.concat(evidence).to_csv(O/'upstream_annotation_evidence.csv',index=False)
summary={'input_level':'raw and corrected processed gene-by-lineage counts; not raw per-cell sequencing or new QC','coding_rank_background':len(common),'known_whole_cell_donors':5,'nuclear_donor_ids_resolved':False,'primary_test_family':10,'conservative_historical_family':18,'primary_significance_threshold':.05,'zero_policy':'reject execution if any zero difference; no zero observed','ties_policy':'average ranks; enumerate all sign patterns; no tied absolute differences observed in these tests','exact_test':'sum of sign-pattern probability where min(W+,W-) <= observed, two-sided','software_scipy_check':True,'raw_cell_reprocessing_rerun':False}
(O/'cell_reconstruction_scope.json').write_text(json.dumps(summary,indent=2));print(t.to_string(index=False));print(summary)
