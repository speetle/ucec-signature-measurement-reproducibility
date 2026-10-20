from common import *
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform
import warnings,time
warnings.filterwarnings('ignore',category=RuntimeWarning)
rng=np.random.default_rng(20261020)
def zscale(p):
 return p.sub(p.mean(axis=1),axis=0).div(p.std(axis=1,ddof=1),axis=0)
def safe(a,b):
 a=np.asarray(a);b=np.asarray(b);g=np.isfinite(a)&np.isfinite(b)
 return rho(a[g],b[g]) if g.sum()>=3 and np.ptp(a[g])>0 and np.ptp(b[g])>0 else np.nan
def delt(a,b,t):
 g=np.isfinite(a)&np.isfinite(b)&np.isfinite(t)
 return safe(b[g],t[g])-safe(a[g],t[g])
def interval(draw):return dict(low=float(np.nanquantile(draw,.025)),high=float(np.nanquantile(draw,.975)),valid_bootstraps=int(np.isfinite(draw).sum()))
def endpoint(pp,up=uo,down=do):return signed(zscale(pp),up,down)
# All decisions frozen before result calculation in analysis_plan.json.
summary=[];batchout=[];joint=[];boots=[];member=[];hm=[];elig=[];hmsummary=[];coverage=[]
old=pd.read_csv(SNAP/'missingness_associations.csv')
remove=old.query('cohort=="discovery95" and comparison=="missingness_by_batch" and BH_within_cohort<0.01').gene.tolist()
(O/'jointly_removed_members.json').write_text(json.dumps({'rule':'discovery batch-associated missingness BH<.01; frozen without external effect selection','genes':remove,'up_retained':[g for g in uo if g not in remove],'down_retained':[g for g in do if g not in remove]},indent=2))
keys=sorted(sets);dist=np.array([[1-len(set(sets[a])&set(sets[b]))/len(set(sets[a])|set(sets[b])) for b in keys] for a in keys]);np.fill_diagonal(dist,0);labels=fcluster(linkage(squareform(dist),method='average'),t=.8,criterion='distance');clmap=dict(zip(keys,labels))
frozen={}
for name in keys:
 full=[g for g in sets[name] if g in R.index and np.nanstd(R.loc[g].values)>0]
 common=[g for g in full if g in P.index and P.loc[g].notna().mean()>=.8 and P.loc[g].std()>0]
 frozen[name]=(full,common)
 for g in full:member.append(dict(set=name,gene=g,protein_common=g in common,discovery_detection=float(P.loc[g].notna().mean()) if g in P.index else 0))
pd.DataFrame(member).to_csv(O/'hallmark_frozen_members.csv',index=False)
for co,(r,p,m) in cohorts.items():
 q=r.rank(pct=True);pp=p.loc[obs];a=signed(q,u,d);b=signed(q,uo,do);t=endpoint(pp);base=delt(a,b,t)
 for mode in ['patient','within_batch']:
  draws=[];ref=[];bf=m['batch' if co=='discovery95' else 'Batch'].astype(str).values;groups=[np.where(bf==v)[0] for v in sorted(set(bf))]
  for it in range(2000):
   j=rng.integers(len(m),size=len(m)) if mode=='patient' else np.concatenate([rng.choice(ix,len(ix),replace=True) for ix in groups])
   fixed=delt(a[j],b[j],t[j]);new=delt(a[j],b[j],endpoint(pp.iloc[:,j]));draws.append(fixed);ref.append(new);boots.append(dict(cohort=co,mode=mode,iteration=it,fixed=fixed,refitted=new))
  for scale,v in [('fixed',draws),('refitted',ref)]:summary.append(dict(cohort=co,mode=mode,scale=scale,n=len(m),delta=base,draws=2000,**interval(v)))
 for unit,col in [('batch','batch' if co=='discovery95' else 'Batch')]+([('plex','Plex')] if co=='external138' else []):
  labs=m[col].astype(str).values
  for v in sorted(set(labs)):
   keep=np.where(labs!=v)[0];inside=np.where(labs==v)[0]
   for scale,tt in [('fixed',t[keep]),('refitted',endpoint(pp.iloc[:,keep]))]:batchout.append(dict(cohort=co,unit=unit,omitted=v,n=len(keep),scale=scale,delta=delt(a[keep],b[keep],tt),within_omitted_n=len(inside),within_omitted_delta=delt(a[inside],b[inside],t[inside]) if len(inside)>=10 else np.nan))
 ou=[g for g in uo if g not in remove];od=[g for g in do if g not in remove];tj=endpoint(p.loc[ou+od],ou,od);bj=signed(q,ou,od)
 for label,bb in [('endpoint_only',b),('joint_RNA_protein',bj)]:
  dr=[]
  for it in range(2000):
   j=rng.integers(len(m),size=len(m));tt=endpoint(p.loc[ou+od].iloc[:,j],ou,od);dr.append(delt(a[j],bb[j],tt))
  joint.append(dict(cohort=co,estimand=label,removed=';'.join(remove),up=len(ou),down=len(od),n=len(m),delta=delt(a,bb,tj),**interval(dr)))
 # Hallmark actual paired contrasts: same bootstrap draws across pathways.
 n=len(m);jj=rng.integers(n,size=(1000,n));mat=[];names=[]
 for name in keys:
  full,common=frozen[name];reason=''
  if len(full)<15 or len(common)<5:reason='insufficient_discovery_members'
  elif not set(full).issubset(r.index) or not set(common).issubset(p.index):reason='frozen_members_absent_externally'
  elif (p.loc[common].std()<=0).any() or p.loc[common].std().isna().any():reason='nonvariable_external_endpoint_member'
  if reason:elig.append(dict(cohort=co,set=name,full=len(full),common=len(common),n=0,eligible=False,reason=reason));continue
  aa=q.loc[full].mean().values;bb=q.loc[common].mean().values;zz=zscale(p.loc[common]);tt=zz.mean().values;good=(zz.notna().sum().values>=int(np.ceil(.8*len(common))));tt[~good]=np.nan
  ng=int(np.isfinite(tt).sum());eligible=ng>=30;elig.append(dict(cohort=co,set=name,full=len(full),common=len(common),n=ng,eligible=eligible,reason='' if eligible else 'fewer_than30_patients'))
  if not eligible:continue
  dr=np.array([delt(aa[j],bb[j],tt[j]) for j in jj]);mat.append(dr);names.append(name)
  hm.append(dict(cohort=co,set=name,cluster=int(clmap[name]),full=len(full),common=len(common),coverage=len(common)/len(full),n=ng,rho_full=safe(aa,tt),rho_common=safe(bb,tt),delta=delt(aa,bb,tt),**interval(dr)))
  coverage.extend([dict(cohort=co,set=name,patient=m.index[i],observed_members=int(zz.notna().sum().iloc[i]),required=int(np.ceil(.8*len(common))),included=bool(good[i])) for i in range(n)])
 mat=np.array(mat);hd=pd.DataFrame([x for x in hm if x['cohort']==co]);cls=hd.cluster.values
 clusterpoint=hd.groupby('cluster').delta.median();clusterboot=np.array([np.nanmedian(mat[cls==c],axis=0) for c in sorted(set(cls))]);aggregates=[('median_delta',hd.delta.median(),np.nanmedian(mat,axis=0)),('mean_of_cluster_medians',clusterpoint.mean(),np.nanmean(clusterboot,axis=0))]
 for typ,point,dr in aggregates:hmsummary.append(dict(cohort=co,summary=typ,eligible=len(hd),positive=int((hd.delta>0).sum()),negative=int((hd.delta<0).sum()),clusters=len(clusterpoint),value=point,**interval(dr),coverage_delta_rho=safe(hd.coverage,hd.delta)))
 pd.DataFrame(mat.T,columns=names).to_csv(O/f'{co}_hallmark_bootstraps.csv.gz',index=False)
 print('completed benchmark/uncertainty',co,flush=True)
pd.DataFrame(summary).to_csv(O/'full_pipeline_bootstrap_summary.csv',index=False);pd.DataFrame(boots).to_csv(O/'full_pipeline_bootstrap_draws.csv.gz',index=False);pd.DataFrame(batchout).to_csv(O/'leave_batch_plex_out.csv',index=False);pd.DataFrame(joint).to_csv(O/'joint_missingness_deletion.csv',index=False);pd.DataFrame(hm).to_csv(O/'hallmark_actual_correspondence.csv',index=False);pd.DataFrame(elig).to_csv(O/'hallmark_eligibility.csv',index=False);pd.DataFrame(hmsummary).to_csv(O/'hallmark_dependence_summary.csv',index=False);pd.DataFrame(coverage).to_csv(O/'hallmark_patient_coverage.csv.gz',index=False)
# DNA-based purity sensitivity uses the same complete-case patients in every model.
m=EM;q=ER.rank(pct=True);a=signed(q,u,d);b=signed(q,uo,do);pp=EP.loc[obs];t=endpoint(pp)
clinical=pd.DataFrame({'age10':pd.to_numeric(m.age,errors='coerce')/10,'grade':m.Histologic_grade.str.extract(r'^(G[123])',expand=False).map({'G1':'Low','G2':'Low','G3':'High'}),'histology':np.where(m.Histologic_Type.str.contains('Endometrioid',case=False,na=False),'Endometrioid','Other'),'subtype':m.Genomic_subtype,'batch':m.Batch},index=m.index)
purity=pd.to_numeric(m.ABSOLUTE_tumor_purity,errors='coerce');purity=purity.where(purity.between(0,1));xc=m[['xCell_immune_score','xCell_stroma_score']].apply(pd.to_numeric,errors='coerce');good=clinical.notna().all(axis=1)&purity.notna()&xc.notna().all(axis=1);idx=np.where(good)[0]
cl=pd.concat([clinical[['age10']],pd.get_dummies(clinical[['grade','histology','subtype','batch']].astype(str),drop_first=True,dtype=float)],axis=1);pur=pd.DataFrame({'ABSOLUTE_purity':purity},index=m.index);adj=[];ys=np.column_stack([a,b,t])[idx];samples=[]
for i in range(len(m)):samples.append(dict(patient=m.index[i],ABSOLUTE_purity=purity.iloc[i],valid_purity=bool(purity.notna().iloc[i]),clinical_complete=bool(clinical.notna().all(axis=1).iloc[i]),xCell_complete=bool(xc.notna().all(axis=1).iloc[i]),shared_intersection=bool(good.iloc[i])))
pd.DataFrame(samples).to_csv(O/'purity_patient_eligibility.csv',index=False)
def resid(y,x):return y-x@np.linalg.lstsq(x,y,rcond=None)[0]
bootj=rng.integers(len(idx),size=(1000,len(idx)))
for label,frame in [('clinical_batch',cl),('clinical_DNA_purity',pd.concat([cl,pur],axis=1)),('clinical_xCell',pd.concat([cl,xc],axis=1)),('clinical_DNA_xCell',pd.concat([cl,pur,xc],axis=1))]:
 x=np.column_stack([np.ones(len(idx)),frame.loc[good].values]);const=np.ptp(x,axis=0)==0;const[0]=False;x=x[:,~const];x[:,1:]=(x[:,1:]-x[:,1:].mean(axis=0))/x[:,1:].std(axis=0)
 rr=resid(ys,x);point=delt(rr[:,0],rr[:,1],rr[:,2]);dr=[];rankdef=0
 for j in bootj:
  yy=ys[j].copy();yy[:,2]=endpoint(pp.iloc[:,idx[j]]);rr=resid(yy,x[j]);dr.append(delt(rr[:,0],rr[:,1],rr[:,2]));rankdef+=np.linalg.matrix_rank(x[j])<x.shape[1]
 adj.append(dict(cohort='external138',model=label,n=len(idx),delta=point,design_columns=x.shape[1],rank_deficient=rankdef,scale='protein refit + covariate refit',**interval(dr)))
pd.DataFrame(adj).to_csv(O/'DNA_purity_sensitivity.csv',index=False)
print('DONE',flush=True)
