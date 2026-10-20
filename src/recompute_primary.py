from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr
ROOT=Path(__file__).resolve().parents[1];I=ROOT/'input';F=ROOT/'frozen';O=ROOT/'results/primary';O.mkdir(parents=True,exist_ok=True)
plan=json.loads((I/'discovery/plan.json').read_text());panel=plan['formal_source_response_panel'];u,d=panel['up'],panel['down'];pd0=pd.read_csv(I/'discovery/protein.csv.gz',index_col=0);obs=[g for g in u+d if g in pd0.index and pd0.loc[g].notna().mean()>=.8];uo=[g for g in u if g in obs];do=[g for g in d if g in obs];assert (len(u),len(d),len(uo),len(do))==(32,16,19,10)
(O/'membership.json').write_text(json.dumps({'full_up':u,'full_down':d,'common_up':uo,'common_down':do},indent=2))
rows=[];errs=[];patient_sets=[]
def corr(Y):
 z=rankdata(Y,axis=0);z-=z.mean(axis=0);z/=np.sqrt(np.sum(z*z,axis=0));return z[:,:3].T@z[:,3]
for co,dir in [('discovery95','discovery'),('external138','external')]:
 m=pd.read_csv(I/dir/('CPTAC_model.csv' if dir=='discovery' else 'metadata.csv'),index_col=0);r=pd.read_csv(I/dir/'RNA.csv.gz',index_col=0)[m.index];p=pd.read_csv(I/dir/'protein.csv.gz',index_col=0)[m.index];patient_sets.append(set(m.index));q=r.rank(pct=True);pz=p.loc[obs].sub(p.loc[obs].mean(axis=1),axis=0).div(p.loc[obs].std(axis=1),axis=0);rz=r.loc[obs].sub(r.loc[obs].mean(axis=1),axis=0).div(r.loc[obs].std(axis=1),axis=0)
 signed=lambda a,x,y:(a.loc[x].mean()-a.loc[y].mean()).values
 full=signed(q,u,d);common=signed(q,uo,do);Y=np.column_stack([full,common,signed(rz,uo,do),signed(pz,uo,do)]);point=corr(Y);rng=np.random.default_rng(2026101801);jj=rng.integers(len(m),size=(20000,len(m)));draw=[]
 for take in np.array_split(jj,20):
  z=rankdata(Y[take],axis=1);z-=z.mean(axis=1,keepdims=True);z/=np.sqrt(np.sum(z*z,axis=1,keepdims=True));cs=np.sum(z[:,:,:3]*z[:,:,3,None],axis=1);draw.append(cs[:,1]-cs[:,0])
 bs=np.concatenate(draw);lo,hi=np.quantile(bs,[.025,.975]);row={'cohort':co,'n':len(m),'rho_full':point[0],'rho_common':point[1],'rho_common_z':point[2],'delta':point[1]-point[0],'low':lo,'high':hi,'resamples':20000,'seed':2026101801};rows.append(row)
 score=pd.DataFrame(Y,index=m.index,columns=['full48','common29','common_z','protein29']);score['original_weight_observed']=len(uo)/len(u)*q.loc[uo].mean()-len(do)/len(d)*q.loc[do].mean();score.to_csv(O/f'{co}_scores.csv',index_label='patient');pd.DataFrame({'iteration':range(20000),'delta':bs}).to_csv(O/f'{co}_primary_draws.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 expected=pd.read_csv(F/'precision_background_results.csv').query('cohort==@co and panel=="historical48" and background=="native" and model=="unadjusted"').iloc[0]
 for a,b in [('rho_full','complete_rank_rho'),('rho_common','common_rank_rho'),('delta','membership_delta'),('low','membership_low'),('high','membership_high')]:errs.append(abs(row[a]-expected[b]))
 print(co,row,flush=True)
assert not patient_sets[0]&patient_sets[1];pd.DataFrame(rows).to_csv(O/'primary_results.csv',index=False)
# Reconstruct antibody-platform point estimates directly from the frozen individual-level scores.
rp=pd.read_csv(F/'RPPA_Hallmark_patient_scores.csv.gz',float_precision='round_trip');reference=pd.read_csv(F/'RPPA_Hallmark_dependent_correlations.csv');rr=[]
for (name,method),z in rp.groupby(['set','method'],sort=False):
 a=spearmanr(z.full_RNA_score,z.common_protein_gene_z_score).statistic;b=spearmanr(z.common_RNA_score,z.common_protein_gene_z_score).statistic;ref=reference.query('set==@name and method==@method').iloc[0];errs.append(abs((b-a)-ref.paired_delta));rr.append({'set':name,'method':method,'n':len(z),'rho_full':a,'rho_common':b,'delta':b-a})
pd.DataFrame(rr).to_csv(O/'RPPA_points.csv',index=False)
validation={'primary_from_processed_matrices':True,'primary_resamples_each':20000,'patient_overlap':0,'RPPA_point_rows_from_individual_scores':len(rr),'max_absolute_error':float(max(errs)),'tolerance':1e-10,'passed':bool(max(errs)<1e-10),'unperformed':'raw FASTQ alignment, raw MS quantification, independent hardware'};(O/'validation.json').write_text(json.dumps(validation,indent=2));assert validation['passed'],validation
