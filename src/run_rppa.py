from pathlib import Path
import zipfile,json,re,hashlib,sys
import pandas as pd,numpy as np
from scipy.stats import rankdata,spearmanr
w=Path(__file__).resolve().parents[1];r=w/'input/rppa';v=w/'input';o=w/'results/rppa';o.mkdir(parents=True,exist_ok=True)
rnafile=r/'ucec_pan_rsem.tsv.gz'
a=pd.read_csv(rnafile,sep='\t');a=a.drop(columns=['Entrez_Gene_Id']).set_index('Hugo_Symbol');a=a.groupby(level=0).median()
def read_s3(path):
 import xml.etree.ElementTree as ET
 ns={'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
 with zipfile.ZipFile(path) as zz:
  shared=[]
  if 'xl/sharedStrings.xml' in zz.namelist():
   shared=[''.join(t.text or '' for t in si.iter('{'+ns['x']+'}t')) for si in ET.fromstring(zz.read('xl/sharedStrings.xml'))]
  wb=ET.fromstring(zz.read('xl/workbook.xml'));sheets=wb.find('x:sheets',ns);sheet=next(s for s in sheets if s.attrib['name']=='Table S3');rid=sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
  rels=ET.fromstring(zz.read('xl/_rels/workbook.xml.rels'));target=next(x.attrib['Target'] for x in rels if x.attrib['Id']==rid);target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
  rows=[]
  for row in ET.fromstring(zz.read(target)).findall('.//x:sheetData/x:row',ns):
   values={}
   for c in row:
    label=re.match(r'[A-Z]+',c.attrib['r']).group();j=0
    for cc in label:j=j*26+ord(cc)-64
    v=c.find('x:v',ns);typ=c.attrib.get('t');val=v.text if v is not None else ''
    if typ=='s':val=shared[int(val)]
    if typ=='inlineStr':val=''.join(t.text or '' for t in c.findall('.//x:t',ns))
    values[j-1]=val
   rows.append(values)
  header=rows[1];width=max(header)+1;names=[header.get(k,'Unnamed_'+str(k)) for k in range(width)]
  return pd.DataFrame([[row.get(k,None) for k in range(width)] for row in rows[2:]],columns=names)
ann=read_s3(r/'RPPA500_supplementary_tables.xlsx')
ann['tokens']=ann.Gene_symbol.fillna('').apply(lambda s:re.split(r'[\s,;/|]+',str(s).strip()))
ann=ann[(ann.Type=='total') & ann.tokens.apply(lambda x:len(x)==1 and x[0]!='')].copy()
ann['gene']=ann.tokens.str[0];ann['matrix_marker']=ann.Protein.astype(str).str.replace(r'[^A-Za-z0-9]','',regex=True).str.upper()
ann.loc[ann.matrix_marker.str.match(r'^\d'),'matrix_marker']='X'+ann.loc[ann.matrix_marker.str.match(r'^\d'),'matrix_marker']
z=zipfile.ZipFile(r/'TCPA_TCGA_RPPA500.zip');b=pd.read_csv(z.open(z.namelist()[0]),sep='\t',index_col=0)
ann['matrix_present']=ann.matrix_marker.isin(b.columns);ann.to_csv(o/'RPPA500_total_mapping_audit.csv',index=False)
ann=ann[ann.matrix_present]
pieces=[]
for gene,g in ann.groupby('gene'):
 pieces.append(pd.Series(b[g.matrix_marker.tolist()].median(axis=1),name=gene))
p=pd.concat(pieces,axis=1).T
rna_ids=sorted(c for c in a.columns if c.startswith('TCGA-') and c[13:15]=='01');prot_ids=sorted(c for c in p.columns if c.startswith('TCGA-') and c[13:15]=='01')
def pick(ids):
 d={};excluded=[]
 for s in ids:
  if s[:12] in d:excluded.append({'patient':s[:12],'excluded_sample':s,'retained_sample':d[s[:12]],'reason':'deterministic lexical-first primary-tumor sample per patient'})
  else:d[s[:12]]=s
 return d,excluded
rs,rex=pick(rna_ids);ps,pex=pick(prot_ids);patients=sorted(set(rs)&set(ps))
audit=pd.DataFrame([{'patient':q,'RNA_sample':rs[q],'RPPA_sample':ps[q],'RNA_sample_type':rs[q][13:15],'RPPA_sample_type':ps[q][13:15],'match_level':'same participant and primary sample type; RNA vial absent in processed label','selected':True} for q in patients]);audit.to_csv(o/'TCGA_RNA_RPPA_patient_matching.csv',index=False)
pd.DataFrame(rex+pex).to_csv(o/'TCGA_duplicate_sample_exclusions.csv',index=False)
A=a[[rs[q] for q in patients]].copy();A.columns=patients;P=p[[ps[q] for q in patients]].copy();P.columns=patients
# gene z is fit to the matched finite cohort and fixed for score-vector bootstrap.
zz=lambda t:t.sub(t.mean(axis=1),axis=0).div(t.std(axis=1,ddof=1).replace(0,np.nan),axis=0)
rank=A.rank(axis=0,pct=True);Z=zz(A);Q=zz(P)
gmt=v/'hallmark.gmt';sets={v[0]:v[2:] for v in [l.rstrip().split('\t') for l in gmt.read_text().splitlines()]}
rows=[];cover=[];scoretables=[];rng=np.random.default_rng(20261018);draw=rng.integers(0,len(patients),(2000,len(patients)))
def rc(x,y):
 x=rankdata(x,axis=1);y=rankdata(y,axis=1);x-=x.mean(axis=1)[:,None];y-=y.mean(axis=1)[:,None];return (x*y).sum(axis=1)/np.sqrt((x*x).sum(axis=1)*(y*y).sum(axis=1))
for name,genes in sets.items():
 full=[g for g in genes if g in A.index and A.loc[g].notna().all() and A.loc[g].std()>0]
 common=[g for g in full if g in Q.index and Q.loc[g].notna().all()]
 eligible=len(common)>=5
 cover.append({'set':name,'defined':len(genes),'RNA_measured_variable_complete':len(full),'RPPA_common_total_genes':len(common),'eligible_ge5':eligible,'common_genes':';'.join(common),'exclusion':'' if eligible else 'fewer_than5_RNA_RPPA_common_total_protein_genes'})
 if not eligible:continue
 endpoint=Q.loc[common].mean().values
 for method,M in [('within_sample_gene_rank',rank),('cohort_gene_z',Z)]:
  fullscore=M.loc[full].mean().values;subset=M.loc[common].mean().values
  rfull=spearmanr(fullscore,endpoint).statistic;rcommon=spearmanr(subset,endpoint).statistic
  diff=rc(subset[draw],endpoint[draw])-rc(fullscore[draw],endpoint[draw]);lo,hi=np.quantile(diff,[.025,.975])
  # Stricter additional interval: refit cohort gene-z transforms in each paired resample.
  # Ranks across genes within one specimen need no refit when patients are resampled.
  def resampled_gene_z_mean(mat,ix):
   q=mat[:,ix];sd=q.std(axis=2,ddof=1);sd[sd==0]=np.nan
   return np.nanmean((q-q.mean(axis=2)[:,:,None])/sd[:,:,None],axis=0)
  refit=[]
  for start in range(0,2000,50):
   ix=draw[start:start+50];ep=resampled_gene_z_mean(P.loc[common].values,ix)
   if method=='within_sample_gene_rank':xx=fullscore[ix];yy=subset[ix]
   else:xx=resampled_gene_z_mean(A.loc[full].values,ix);yy=resampled_gene_z_mean(A.loc[common].values,ix)
   refit.extend(rc(yy,ep)-rc(xx,ep))
  rlo,rhi=np.nanquantile(refit,[.025,.975])
  rows.append({'set':name,'method':method,'n_patients':len(patients),'RNA_full_measured':len(full),'RPPA_common_total_genes':len(common),'full_RNA_protein_rho':rfull,'common_RNA_protein_rho':rcommon,'paired_delta':rcommon-rfull,'delta_low':lo,'delta_high':hi,'delta_refit_low':rlo,'delta_refit_high':rhi,'bootstrap':2000,'scope':'descriptive dependent correlation contrasts; delta_low/high conditional fixed scales, delta_refit_low/high restandardized within each bootstrap; intervals unadjusted; not source48 validation'})
  scoretables.extend({'set':name,'method':method,'patient':q,'full_RNA_score':fullscore[i],'common_RNA_score':subset[i],'common_protein_gene_z_score':endpoint[i]} for i,q in enumerate(patients))
pd.DataFrame(cover).to_csv(o/'RPPA_Hallmark_coverage.csv',index=False);D=pd.DataFrame(rows);D.to_csv(o/'RPPA_Hallmark_dependent_correlations.csv',index=False);pd.DataFrame(scoretables).to_csv(o/'RPPA_Hallmark_patient_scores.csv.gz',index=False)
summary={'matched_patients':len(patients),'RNA_primary_samples':len(rna_ids),'RNA_unique_primary_patients':len(rs),'RPPA_all_cancer_primary_samples':len(prot_ids),'RPPA_unique_total_genes':len(P),'eligible_pathways':int(sum(v['eligible_ge5'] for v in cover)),'all_pathways':50,'bootstrap_seed':20261018,'bootstrap_interval_columns':'delta_refit_low/high refit all cohort gene-z transforms within paired resamples; delta_low/high retained as conditional fixed-scale sensitivity','RNA_range':[float(np.nanmin(A.values)),float(np.nanmax(A.values))],'RNA_transform':'none; already batch-normalized RSEM, negative values retained','limitations':['participant ID and primary-tumor type match; exact aliquot/vial match cannot be established from RNA labels','TCGA cohort not new relative to prior TCGA survival analysis, but independent source from CPTAC','RNA and total-antibody RPPA endpoints; not whole-proteome measurements','intervals are not multiple-testing-adjusted; no per-pathway significance claims','same signatures are dependent; summary medians have no cross-pathway inferential p value','this benchmark does not validate the mouse-derived48 or51 panel, which have zero RPPA gene coverage']}
summary['by_method']={m:{'n_pathways':len(g),'median_full_rho':g.full_RNA_protein_rho.median(),'median_common_rho':g.common_RNA_protein_rho.median(),'median_paired_delta':g.paired_delta.median(),'delta_positive':int((g.paired_delta>0).sum()),'delta_negative':int((g.paired_delta<0).sum()),'delta_min':g.paired_delta.min(),'delta_max':g.paired_delta.max()} for m,g in D.groupby('method')}
(o/'RPPA_Hallmark_summary.json').write_text(json.dumps(summary,indent=2))
files=[rnafile,gmt,Path(__file__)]
(o/'RPPA_Hallmark_input_sha256.json').write_text(json.dumps({str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},indent=2))
print(json.dumps(summary,indent=2))
