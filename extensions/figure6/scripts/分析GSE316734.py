from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,itertools
from scipy.stats import spearmanr
P=Path(__file__).resolve().parents[1];p=P/'input'/'GSE316734_Raw_count_file_rsemgenes.txt.gz';hpath=P.parents[1]/'input'/'hgnc.txt';mapping=pd.read_csv(P/'assets'/'Tahoe源48及共有29基因映射.csv');(P/'results').mkdir(exist_ok=True);assert len(mapping)==48
cols=['1_AN3CA_1','2_AN3CA_2','3_AN3CA_3','4_AN3CA_DHX9i_1','5_AN3CA_DHX9i_2','6_AN3CA_DHX9i_3'];d=pd.read_csv(p,sep='\t',usecols=['gene_id']+cols);assert d.gene_id.is_unique and d[cols].notna().all().all() and (d[cols]>=0).all().all();d['ensembl']=d.gene_id.str.split('_',n=1).str[0].str.split('.').str[0]
h=pd.read_csv(hpath,sep='\t',dtype=str).dropna(subset=['ensembl_gene_id','symbol']);ng=h.groupby('ensembl_gene_id').symbol.nunique();ns=h.groupby('symbol').ensembl_gene_id.nunique();h=h[h.ensembl_gene_id.isin(ng[ng==1].index)&h.symbol.isin(ns[ns==1].index)].drop_duplicates('ensembl_gene_id');e=h.set_index('ensembl_gene_id').symbol;d['gene']=d.ensembl.map(e);d=d[d.gene.notna()];duplicated=d.gene.duplicated(keep=False);d=d[~duplicated];d=d[d[cols].sum(axis=1)>0].set_index('gene');assert d.index.is_unique
rank=d[cols].rank(axis=0,pct=True);results={};missing={}
for name,g in [('source48',mapping),('common29',mapping[mapping.member_common29.astype(str)=='True'])]:
 up=[x for x in g[g.direction=='up'].gene if x in rank.index];down=[x for x in g[g.direction=='down'].gene if x in rank.index];assert len(up)>=2 and len(down)>=2
 missing[name]=sorted(set(g.gene)-set(up)-set(down));results[name]=rank.loc[up].mean()-rank.loc[down].mean()
out=pd.DataFrame(results);out['condition']=['DMSO']*3+['ATX968']*3;out.index.name='sample';out.to_csv(P/'results'/'GSE316734_AN3两个面板样本分数.csv',encoding='utf-8-sig')
x=out.source48.to_numpy();y=out.common29.to_numpy();discord=ties=0
for i,j in itertools.combinations(range(6),2):
 a=np.sign(x[i]-x[j]);b=np.sign(y[i]-y[j]);ties+=int(a==0 or b==0);discord+=int(a*b<0)
perms=list(itertools.combinations(range(6),3));metric={}
for name in ['source48','common29']:
 v=out[name].to_numpy();delta=float(v[3:].mean()-v[:3].mean());null=[]
 for treatment in perms:
  ind=np.array([i in treatment for i in range(6)]);null.append(v[ind].mean()-v[~ind].mean())
 metric[name]={'treated_minus_control':delta,'exact_two_sided_permutation_p':sum(abs(t)>=abs(delta)-1e-12 for t in null)/len(null),'permutations':len(null)}
meta={'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'samples':cols,'mapped_nonzero_genes':len(d),'ambiguous_gene_rows_removed':int(duplicated.sum()),'missing_panel_genes':missing,'rank_spearman':float(spearmanr(x,y).statistic),'discordant_pairs':discord,'tied_pairs':ties,'untied_pairs':15-ties,'treatment_effects':metric,'limitation':'one non-overlapping drug and six independent culture samples per GEO design; exact permutation assumes label exchangeability, no causal or clinical inference'}
(P/'results'/'GSE316734_AN3覆盖与扰动摘要.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False));print(json.dumps(meta,indent=2,ensure_ascii=False))
