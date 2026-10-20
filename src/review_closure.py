"""Descriptive reviewer-requested audits; no post-hoc significance selection."""
from controls_common import *
from scipy.stats import spearmanr
import patsy
A=ROOT/'results/review';A.mkdir(parents=True,exist_ok=True)
# M4. Fixed slot match balance and error-conditional external reference summaries.
bal=pd.read_csv(O/'coexpression_matching_balance.csv');selection=pd.read_csv(O/'coexpression_panel_selection.csv');null=pd.read_csv(O/'coexpression_null_results.csv')
bal.groupby(['ensemble','feature']).mean_absolute_standardized_difference.agg(['count','median','min','max']).reset_index().to_csv(A/'matching_feature_balance.csv',index=False)
null.groupby(['ensemble','cohort']).covariance_RMSE.agg(['count','median','min','max']).reset_index().to_csv(A/'matching_covariance_balance.csv',index=False)
rows=[]
primary=pd.read_csv(ROOT/'results/primary/primary_results.csv')
for (ens,co),v in null.groupby(['ensemble','cohort']):
 v=v.copy();v['error_quartile']=pd.qcut(v.covariance_RMSE,4,labels=False)+1;obsdelta=float(primary.query('cohort==@co').delta.iloc[0])
 for quartile,g in v.groupby('error_quartile'):
  rows.append(dict(ensemble=ens,cohort=co,error_quartile=quartile,panels=len(g),min_error=g.covariance_RMSE.min(),max_error=g.covariance_RMSE.max(),median_error=g.covariance_RMSE.median(),median_delta=g.delta.median(),observed_delta=obsdelta,tail_plus1=(1+sum(g.delta>=obsdelta))/(len(g)+1),interpretation='dependent post-results descriptive strata, not calibrated P or equivalence test'))
pd.DataFrame(rows).to_csv(A/'matching_error_strata.csv',index=False)
# Candidate universe and exact signed/module/detection slot eligibility, no effect-based external filtering.
f=pd.read_csv(ROOT/'results/null/discovery_matching_features.csv',index_col=0);modules=pd.read_csv(O/'discovery_coexpression_modules.csv',index_col=0).iloc[:,0];eligible=f.index.difference(u+d).intersection(ER.index);ep=EP.index[EP.notna().any(axis=1)];availability={True:eligible[(f.loc[eligible,'protein_detection']>=.8)&eligible.isin(ep)],False:eligible[f.loc[eligible,'protein_detection']<.8]};rows=[]
for g in u+d:
 pool=availability[g in obs];mpool=pool[modules.loc[pool].values==modules[g]]
 rows.append(dict(template=g,direction='up' if g in u else 'down',protein_common=g in obs,feature_candidate_pool=len(pool),module_candidate_pool=len(mpool),used_nearest_candidates=min(100,len(mpool)),source_excluded=True,external_filter='RNA row presence; observed slots additionally any external protein observation; no external correlation/outcome'))
pd.DataFrame(rows).to_csv(A/'matching_slot_eligibility.csv',index=False)
# M5. Describe all eligible AND unavailable pathways without restoring missing external members.
elig=pd.read_csv(ROOT/'results/ms/hallmark_eligibility.csv');members=pd.read_csv(ROOT/'results/ms/hallmark_frozen_members.csv');hm=pd.read_csv(ROOT/'results/ms/hallmark_actual_correspondence.csv');rows=[]
for name,g in members.groupby('set',sort=True):
 full=g.gene.tolist();common=g.loc[g.protein_common,'gene'].tolist();ext=elig.query('cohort=="external138" and set==@name').iloc[0];absent=[x for x in common if x not in EP.index]
 rows.append(dict(set=name,external_eligible=bool(ext.eligible),reason=ext.reason,full_members=len(full),frozen_protein_members=len(common),discovery_coverage=len(common)/len(full),discovery_median_RNA_expression=float(R.loc[full].mean(axis=1).median()),discovery_median_RNA_variance=float(R.loc[full].var(axis=1).median()),discovery_median_detection_all_members=g.discovery_detection.median(),external_absent_frozen_members=len(absent),external_absent_genes=';'.join(absent),external_n=int(ext.n)))
q=pd.DataFrame(rows);q.to_csv(A/'pathway_eligibility_all50.csv',index=False);vars=['full_members','frozen_protein_members','discovery_coverage','discovery_median_RNA_expression','discovery_median_RNA_variance','discovery_median_detection_all_members','external_absent_frozen_members'];rows=[]
for eligible,g in q.groupby('external_eligible'):
 for col in vars:rows.append(dict(external_eligible=bool(eligible),sets=len(g),feature=col,median=g[col].median(),q25=g[col].quantile(.25),q75=g[col].quantile(.75),min=g[col].min(),max=g[col].max(),inference='descriptive only; overlapping pathways not independent replicates'))
pd.DataFrame(rows).to_csv(A/'eligible31_vs_unavailable19.csv',index=False)
x=hm.query('cohort=="external138"').copy();x['coverage_tertile']=pd.qcut(x.coverage,3,labels=False)+1
x.groupby('coverage_tertile').agg(sets=('delta','size'),positive=('delta',lambda z:sum(z>0)),negative=('delta',lambda z:sum(z<0)),median_delta=('delta','median'),min_delta=('delta','min'),max_delta=('delta','max'),median_coverage=('coverage','median')).reset_index().to_csv(A/'pathway_coverage_heterogeneity.csv',index=False)
# M6. Reconstruct every displayed adjusted point estimate, and report exact point design columns and ranks.
# The pre-existing bootstrap intervals remain separately verified by their own original computations.
points=[];designs=[]
def analyze(co,model,m,frame,good,r,p,bootstrap_source=None):
 ids=m.index[good];X=np.column_stack([np.ones(len(ids)),frame.loc[good].to_numpy(dtype=float)]);const=np.ptp(X,axis=0)==0;const[0]=False;X=X[:,~const];X[:,1:]=(X[:,1:]-X[:,1:].mean(axis=0))/X[:,1:].std(axis=0)
 ranks=r.rank(pct=True);pz=p.loc[obs].sub(p.loc[obs].mean(axis=1),axis=0).div(p.loc[obs].std(axis=1),axis=0);Y=np.column_stack([signed(ranks,u,d),signed(ranks,uo,do),signed(pz,uo,do)]);Y=Y[m.index.get_indexer(ids)];Z=Y-X@np.linalg.lstsq(X,Y,rcond=None)[0];c1=rho(Z[:,0],Z[:,2]);c2=rho(Z[:,1],Z[:,2]);row=dict(cohort=co,model=model,n=len(ids),design_columns=X.shape[1],design_rank=int(np.linalg.matrix_rank(X)),rho_full=c1,rho_common=c2,delta=c2-c1);points.append(row);designs.append(row.copy())
# Native detailed models reproduce Patsy's categorical design, retaining zero columns where present.
for co,(r,p,m0) in cohorts.items():
 m=m0.copy()
 if co=='external138':m['age10']=pd.to_numeric(m.age,errors='coerce')/10;m['grade']=m.Histologic_grade;m['histology']=m.Histologic_Type;m['subtype']=m.Genomic_subtype;m['batch']=m.Batch
 forms={'clinical_subtype_batch':'age10+C(grade)+C(histology)+C(subtype)+C(batch)'}
 forms.update({'fully_adjusted_original':forms['clinical_subtype_batch']+'+immune_logratio+stroma_logratio+proliferation+inflammation+interferon+hypoxia+oxidative_stress+dissociation'} if co=='discovery95' else {'xCell_sensitivity':forms['clinical_subtype_batch']+'+xCell_immune_score+xCell_stroma_score'})
 ranks=r.rank(pct=True);pz=p.loc[obs].sub(p.loc[obs].mean(axis=1),axis=0).div(p.loc[obs].std(axis=1),axis=0);rz=r.loc[obs].sub(r.loc[obs].mean(axis=1),axis=0).div(r.loc[obs].std(axis=1),axis=0);Y=pd.DataFrame(np.column_stack([signed(ranks,u,d),signed(ranks,uo,do),signed(rz,uo,do),signed(pz,uo,do)]),index=m.index)
 for model,formula in forms.items():
  X=patsy.dmatrix(formula,m,return_type='dataframe');Z=Y.loc[X.index].values-X.values@np.linalg.lstsq(X.values,Y.loc[X.index].values,rcond=None)[0];c1=float(spearmanr(Z[:,0],Z[:,3]).statistic);c2=float(spearmanr(Z[:,1],Z[:,3]).statistic);ref=pd.read_csv(ROOT/'results/precision/precision_background_results.csv').query('cohort==@co and panel=="historical48" and background=="native" and model==@model').iloc[0]
  row=dict(cohort=co,model=model,n=len(X),design_columns=X.shape[1],design_rank=int(np.linalg.matrix_rank(X)),rho_full=c1,rho_common=c2,delta=c2-c1,bootstrap_draws=int(ref.requested),rank_deficient_draws=int(ref.rank_deficient_draws_retained),bootstrap_record='precision_background_results.csv');designs.append(row);points.append(row)
# Simplified external points using the EXACT original factor definitions; this is not the purity model's separate design.
m=EM.copy();m['age10']=pd.to_numeric(m.age,errors='coerce')/10;m['grade_binary']=m.Histologic_grade.map(lambda x:'low_intermediate' if str(x).startswith(('G1','G2')) else ('high' if str(x).startswith('G3') else np.nan));m['hist_binary']=np.where(m.Histologic_Type=='Endometrioid carcinoma','endometrioid','other');m['subtype']=m.Genomic_subtype;m['batch']=m.Batch
for model,form in [('minimal','age10+C(subtype)+C(batch)'),('binary_grade_histology','age10+C(subtype)+C(batch)+C(grade_binary)+C(hist_binary)')]:
 X=patsy.dmatrix(form,m,return_type='dataframe');ranks=ER.rank(pct=True);pz=EP.loc[obs].sub(EP.loc[obs].mean(axis=1),axis=0).div(EP.loc[obs].std(axis=1),axis=0);rz=r.loc[obs].sub(r.loc[obs].mean(axis=1),axis=0).div(r.loc[obs].std(axis=1),axis=0);Y=pd.DataFrame(np.column_stack([signed(ranks,u,d),signed(ranks,uo,do),signed(rz,uo,do),signed(pz,uo,do)]),index=m.index);Z=Y.loc[X.index].values-X.values@np.linalg.lstsq(X.values,Y.loc[X.index].values,rcond=None)[0];c1=float(spearmanr(Z[:,0],Z[:,3]).statistic);c2=float(spearmanr(Z[:,1],Z[:,3]).statistic);ref=pd.read_csv(ROOT/'results/precision/reduced_model_sensitivity.csv').query('model==@model').iloc[0];row=dict(cohort='external138',model=model,n=len(X),design_columns=X.shape[1],design_rank=int(np.linalg.matrix_rank(X)),rho_full=c1,rho_common=c2,delta=c2-c1,bootstrap_draws=int(ref.draws),rank_deficient_draws=int(ref.rank_deficient_draws_retained),bootstrap_record='reduced_model_sensitivity.csv');designs.append(row);points.append(row)
# Purity models use a fixed shared 128-patient intersection, and the same standardized frame as the re-run MS module.
m=EM;clinical=pd.DataFrame({'age10':pd.to_numeric(m.age,errors='coerce')/10,'grade':m.Histologic_grade.str.extract(r'^(G[123])',expand=False).map({'G1':'Low','G2':'Low','G3':'High'}),'histology':np.where(m.Histologic_Type.str.contains('Endometrioid',case=False,na=False),'Endometrioid','Other'),'subtype':m.Genomic_subtype,'batch':m.Batch},index=m.index);purity=pd.to_numeric(m.ABSOLUTE_tumor_purity,errors='coerce');purity=purity.where(purity.between(0,1));xc=m[['xCell_immune_score','xCell_stroma_score']].apply(pd.to_numeric,errors='coerce');good=clinical.notna().all(axis=1)&purity.notna()&xc.notna().all(axis=1);cl=pd.concat([clinical[['age10']],pd.get_dummies(clinical[['grade','histology','subtype','batch']].astype(str),drop_first=True,dtype=float)],axis=1);pur=pd.DataFrame({'ABSOLUTE_purity':purity},index=m.index)
for model,frame in [('clinical_batch',cl),('clinical_DNA_purity',pd.concat([cl,pur],axis=1)),('clinical_xCell',pd.concat([cl,xc],axis=1)),('clinical_DNA_xCell',pd.concat([cl,pur,xc],axis=1))]:
 analyze('external138',model,m,frame,good,ER,EP);ref=pd.read_csv(ROOT/'results/ms/DNA_purity_sensitivity.csv').query('model==@model').iloc[0];designs[-1].update(bootstrap_draws=int(ref.valid_bootstraps),rank_deficient_draws=int(ref.rank_deficient),bootstrap_record='rerun DNA_purity_sensitivity.csv')
pd.DataFrame(points).to_csv(A/'adjusted_point_correlations.csv',index=False);des=pd.DataFrame(designs);des['rank_deficient_fraction']=des.rank_deficient_draws/des.bootstrap_draws;des.to_csv(A/'all_adjustment_design_audit.csv',index=False)
# M1. Identical participants, explicitly different objects. No biologic truth or decision threshold is invented.
panel=pd.read_csv(ROOT/'results/masks/real_panel_factorial_metrics.csv');case=panel.query('cohort=="external138" and signature=="historical48" and member_state=="observed" and background=="measured_panel"');case.to_csv(A/'operational_real_data_example.csv',index=False)
# All point estimates compared to prior records are numeric consistency checks, not new confidence intervals.
checks=[]
for row in points:
 if row['model'] in ['binary_grade_histology','minimal']:
  ref=pd.read_csv(ROOT/'results/precision/reduced_model_sensitivity.csv').query('model==@row["model"]').iloc[0];error=abs(row['delta']-ref.membership_delta)
 elif row['model'] in ['clinical_batch','clinical_DNA_purity','clinical_xCell','clinical_DNA_xCell']:
  ref=pd.read_csv(ROOT/'results/ms/DNA_purity_sensitivity.csv').query('model==@row["model"]').iloc[0];error=abs(row['delta']-ref.delta)
 else:
  ref=pd.read_csv(ROOT/'results/precision/precision_background_results.csv').query('cohort==@row["cohort"] and panel=="historical48" and background=="native" and model==@row["model"]').iloc[0];error=abs(row['delta']-ref.membership_delta)
 checks.append(dict(cohort=row['cohort'],model=row['model'],error=error))
(A/'review_closure_validation.json').write_text(json.dumps({'point_checks':checks,'all_passed':all(x['error']<1e-10 for x in checks),'descriptive_qualification_comparisons':True,'prospective_or_causal_validation':False},indent=2));assert all(x['error']<1e-10 for x in checks)
print('Reviewer requested descriptive audits completed',flush=True)
