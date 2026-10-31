from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,str(Path(sys.argv[1]).resolve()/'src'))
from audit_panel_alignment import require_matplotlib_panel_alignment
ROOT=Path(sys.argv[1]);F=ROOT/'frozen';OUT=Path(sys.argv[2]) if len(sys.argv)>2 else Path(__file__).parent.parent/'figures';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],'font.size':7,'axes.labelsize':6.5,'axes.titlesize':7,'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':5.8,'legend.frameon':False,'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none'})
BLUE='#35688A';ORANGE='#C97742';GRAY='#929AA1';GREEN='#508F86';PURPLE='#967EA1';OLIVE='#8C9667';sources={}
files={}
for folder,names in {
 'primary':['primary_results.csv','discovery95_scores.csv','external138_scores.csv'],
 'null':['weight_contrasts_bootstrap.csv','within_batch_permutations.csv','null_summary.csv'],
 'controls':['coexpression_null_results.csv','coexpression_null_summary.csv','protein29_jackknife.csv','LOXL3_assumption_stress.csv'],
 'ms':['hallmark_actual_correspondence.csv','hallmark_shared_scope_summary.csv','hallmark_dependence_summary.csv','hallmark_external_missing_members.csv'],
 'rppa':['RPPA_Hallmark_dependent_correlations.csv'],
 'masks':['real_panel_factorial_metrics.csv','frozen_transfer_summary_by_family.csv','signature_specific_summary.csv'],
 'arrays':['GSE17025_tumor_primary_scores.csv','GSE63678_tumor_primary_scores.csv'],
 'heldout':['heldout_score_models.csv'],
 'cells':['admixture_curves.csv','lineage_pseudobulk.csv','exact_donor_tests.csv'],
 'clinical':['external_survival_exploratory.csv','source48_increment_bootstrap.csv','source48_decision_summary.csv','TCGA_source48_PH_sensitivity.csv']}.items():
 for name in names:files[name]=ROOT/'results'/folder/name
def read(n):
 p=files[n];sources[n]=hashlib.sha256(p.read_bytes()).hexdigest();return pd.read_csv(p,float_precision='round_trip')
def new():
 sources.clear()
 f=plt.figure(figsize=(180/25.4,230/25.4));axs=[]
 for k,(x,y) in enumerate([(28,169),(116,169),(28,101),(116,101),(28,33),(116,33)]):
  a=f.add_axes([x/180,y/230,58/180,42/230]);a.set_gid(chr(97+k));a.text(-.16,1.13,chr(97+k),transform=a.transAxes,fontsize=9,fontweight='bold');axs.append(a)
 return f,axs
contracts={1:'Does restricting RNA members change correspondence with the same protein endpoint in independent patients?',2:'Do structural references, patient mismatch and member influence support a distinctive source-panel interpretation?',3:'How broadly does member alignment generalize across pathways and protein assays?',4:'How do measured members, signed branches and ranking background alter ordering?',5:'Which composition and outcome checks delimit interpretation without validating member changes?'}
def save(f,axs,n,note):
 f.text(.155,.025,note,fontsize=6,color='#454D53')
 require_matplotlib_panel_alignment(f,axes=axs,panel_ids=list('abcdef'),row_groups=([list('ad'),list('be'),list('cf')] if n==5 else [list('ab'),list('cd'),list('ef')]),column_groups=([list('abc'),list('def')] if n==5 else [list('ace'),list('bdf')]),require_panel_labels=False,json_out=OUT/f'Figure{n}_alignment.json',strict=True,tolerance_pt=1.5)
 f.savefig(OUT/f'Figure{n}.pdf');f.savefig(OUT/f'Figure{n}.svg');f.savefig(OUT/f'Figure{n}.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
 f.savefig(OUT/f'Figure{n}.png',dpi=600);f.savefig(OUT/f'Figure{n}_preview.png',dpi=300)
 (OUT/f'Figure{n}_sources.json').write_text(json.dumps({'question':contracts[n],'archetype':'quantitative grid','backend':'Python matplotlib','dimensions_mm':[180,230],'source_sha256':sources.copy(),'all_data_retained':True},indent=2));plt.close(f)
def forest(ax,labels,val,lo,hi,colors,zero=0):
 for i,(v,l,h,c) in enumerate(zip(val,lo,hi,colors)):ax.errorbar(v,i,xerr=[[v-l],[h-v]],fmt='o',capsize=2,color=c,ms=3,lw=.8)
 ax.set_yticks(range(len(labels)),labels);ax.invert_yaxis();ax.axvline(zero,color=GRAY,lw=.65,ls='--')
f,axs=new();a,b,c,d,e,ff=axs
from matplotlib.patches import FancyBboxPatch
# Figure 1 compares endpoint correlation with preservation of the original RNA order.
a.set_axis_off()
for yy,lab,col in [(0.82,'full48 RNA | 32 up, 16 down',BLUE),(0.49,'common29 RNA | 19 up, 10 down',ORANGE),(0.16,'fixed29 protein | partial endpoint',GREEN)]:
 a.add_patch(FancyBboxPatch((.02,yy-.085),.96,.18,boxstyle='round,pad=0.01',facecolor=col+'18',edgecolor=col,lw=.8,transform=a.transAxes));a.text(.5,yy,lab,ha='center',va='center',transform=a.transAxes,fontsize=6.5)
a.text(.5,.995,'Same patients; different measurement objects',ha='center',transform=a.transAxes,fontsize=6.5)
a.text(.5,-.16,'Biological construct accuracy is not measured',ha='center',transform=a.transAxes,fontsize=5.8,color='#454D53')
pr=read('primary_results.csv');forest(b,['Discovery 95','External 138'],pr.delta,pr.low,pr.high,[BLUE,ORANGE]);b.set(xlabel='Common-minus-full Spearman rho',title='Fixed protein endpoint; 20,000 draws',xlim=(-.05,.33))
newpath=ROOT/'results/v17_review/main_ordering_tradeoff.csv';tt=pd.read_csv(newpath).set_index('cohort');sources['v17/main_ordering_tradeoff.csv']=hashlib.sha256(newpath.read_bytes()).hexdigest()
for co,ax,col in [('discovery95',c,BLUE),('external138',d,ORANGE)]:
 q=read(f'{co}_scores.csv');z=tt.loc[co];ax.scatter(q.full48,q.common29,s=8,color=col,alpha=.7);ax.set(xlabel='Full48 RNA rank score',ylabel='Common29 RNA rank score',title=f'{int(z.n)} patients: RNA order rho = {z.full_common_rho:.3f}')
labels=['D95','E138'];colors=[BLUE,ORANGE];z=tt.loc[['discovery95','external138']]
e.bar(labels,z.discordance_percent,color=colors,width=.5);e.set(ylabel='Discordant untied patient pairs (%)',title='RNA ordering changes',ylim=(0,28))
for j,row in enumerate(z.itertuples()):e.text(j,row.discordance_percent+1,f'{int(row.discordant_pairs)}/{int(row.untied_pairs)}',ha='center',fontsize=6)
ff.bar(labels,z.retention_percent,color=colors,width=.5);ff.set(ylabel='Full48 top-quartile patients retained (%)',title='Highest-quartile retention',ylim=(0,108))
for j,row in enumerate(z.itertuples()):ff.text(j,row.retention_percent+3,f'{int(row.retained)}/{int(row.top_n)}',ha='center',fontsize=6)
save(f,axs,1,'Higher protein correlation and preservation of the full48 RNA order are separate empirical questions.')
f,axs=new();a,b,c,d,e,ff=axs
pr=read('primary_results.csv')
w=read('weight_contrasts_bootstrap.csv').query('model=="unadjusted" and metric in ["deletion_step","reweight_step"]');labs=[];vals=[];low=[];high=[];cols=[]
for co,col in [('discovery95',BLUE),('external138',ORANGE)]:
 for metric,lab in [('deletion_step','Remove'),('reweight_step','Reweight')]:
  z=w.query('cohort==@co and metric==@metric').iloc[0];labs.append(('D95' if col==BLUE else 'E138')+' '+lab);vals.append(z.value);low.append(z.low);high.append(z.high);cols.append(col)
forest(a,labs,vals,low,high,cols);a.set(xlabel='Sequential change in Spearman rho',title='Removal versus branch renormalization')
nr=read('coexpression_null_results.csv');ens=['feature_only','module_matched','covariance_constrained'];v=[nr.query('cohort=="external138" and ensemble==@en').delta for en in ens];b.boxplot(v,vert=False,widths=.5,showfliers=True,flierprops={'markersize':1.2,'alpha':.3},medianprops={'color':BLUE});b.set_yticks([1,2,3],['Features','+ Module','+ Covariance']);b.invert_yaxis();b.axvline(float(pr.query('cohort=="external138"').delta.iloc[0]),color=ORANGE,lw=1);b.set(xlabel='Common-minus-full rho',title='Ordinary panels; external patients')
pe=read('within_batch_permutations.csv');ns=read('null_summary.csv')
for ax,co in [(c,'discovery95'),(d,'external138')]:
 z=pe[pe.cohort.eq(co)];ax.hist(z.total_delta,bins=28,color=GRAY,edgecolor='white',lw=.3);ax.axvline(z.observed_delta.iloc[0],color=ORANGE,lw=1);p=ns.query('cohort==@co and null_type=="within_batch_permutation"').upper_tail_fraction_plus1.iloc[0];ax.set(xlabel='Common-minus-full rho',ylabel='Permutations',title=f'{co}: within-batch mismatch');ax.set_title(f'{co}: within-batch mismatch\nUpper-tail fraction {p:.4f}',fontsize=6.7,pad=7)
j=read('protein29_jackknife.csv').query('cohort=="external138"');e.scatter(j.endpoint_only_delta,j.joint_member_deletion_delta,color=BLUE,s=10);e.set(xlabel='Delete from protein endpoint only',ylabel='Delete from RNA and protein',title='All 29 external member deletions')
for ix in [j.endpoint_only_delta.idxmin(),j.endpoint_only_delta.idxmax()]:
 zz=j.loc[ix];e.annotate(zz.gene,(zz.endpoint_only_delta,zz.joint_member_deletion_delta),xytext=(4,3),textcoords='offset points',fontsize=5.5)
mp=j.copy().reset_index(drop=True);mp.insert(0,'point_id',range(1,len(mp)+1));mp.to_csv(OUT/'Figure2e_member_points.csv',index=False)
l=read('LOXL3_assumption_stress.csv').query('cohort=="external138"');ff.scatter(range(len(l)),l.delta,color=BLUE,s=12);ff.axhline(0,color=GRAY,lw=.6);ff.set(xlabel='Twelve assumed missing-z scenarios',ylabel='Common-minus-full rho',title='LOXL3 stress; values were not observed',xticks=[0,3,6,9,11],ylim=(-.01,.24))
lm=l.copy().reset_index(drop=True);lm.insert(0,'point_id',range(len(lm)));lm.to_csv(OUT/'Figure2f_scenario_points.csv',index=False)
save(f,axs,2,'Reference panels and stress scenarios are dependent computational comparisons, not new biological replicates.')
f,axs=new();a,b,c,d,e,ff=axs
hm=read('hallmark_actual_correspondence.csv')
for i,(co,col) in enumerate([('discovery95',BLUE),('external138',ORANGE)]):
 z=hm[hm.cohort.eq(co)].delta.to_numpy();a.boxplot([z],positions=[i],widths=.45,showfliers=False,medianprops={'color':col});a.scatter(i+np.linspace(-.13,.13,len(z)),z,s=8,color=col,alpha=.7)
a.set_xticks([0,1],['D95: 50 sets','E138: 31 sets']);a.axhline(0,color=GRAY,lw=.6,ls='--');a.set(ylabel='Common-minus-full rho',title='Actual mass-spectrometry benchmarks')
x=hm.pivot(index='set',columns='cohort',values='delta').dropna();assert len(hm)==81 and len(x)==31;b.scatter(x.discovery95,x.external138,s=10,color=BLUE);b.plot([-.15,.42],[-.15,.42],color=GRAY,lw=.6,ls='--');b.axhline(0,color=GRAY,lw=.5);b.axvline(0,color=GRAY,lw=.5);b.set(xlabel='Discovery contrast',ylabel='External contrast',title='Same 31 transferable Hallmarks')
ss=read('hallmark_shared_scope_summary.csv');hd=read('hallmark_dependence_summary.csv');z=hd.query('cohort=="discovery95" and summary=="median_delta"').iloc[0];forest(c,['D95: all 50','D95: same 31','E138: same 31'],[z.value,*ss.median_delta],[z.low,*ss.low],[z.high,*ss.high],[BLUE,BLUE,ORANGE]);c.set(xlabel='Median contrast; patient-bootstrap interval',title='Matched scope limits the general claim')
rp=read('RPPA_Hallmark_dependent_correlations.csv').query('method=="within_sample_gene_rank"');d.scatter(rp.full_RNA_protein_rho,rp.common_RNA_protein_rho,color=GREEN,s=11);d.plot([-.25,.85],[-.25,.85],ls='--',color=GRAY,lw=.6);d.set(xlabel='Full RNA versus same RPPA endpoint',ylabel='Common RNA versus same endpoint',title='RPPA: 422 patients; 25 other pathways')
e.scatter(100*rp.RPPA_common_total_genes/rp.RNA_full_measured,rp.paired_delta,color=GREEN,s=11);e.axhline(0,color=GRAY,lw=.6);e.set(xlabel='Measured RNA members with RPPA (%)',ylabel='Common-minus-full rho',title='Other Hallmarks: RPPA coverage and contrast')
miss=read('hallmark_external_missing_members.csv');cnt=miss.missing_external_protein.fillna('').apply(lambda s:len(s.split(';')) if s else 0);vv=cnt[cnt>0];ff.hist(vv,bins=np.arange(.5,vv.max()+1.5),color=GRAY,edgecolor='white');ff.set(xlabel='Absent frozen external protein rows / set',ylabel='Hallmark sets',title=f'{len(vv)} sets were non-transferable')
save(f,axs,3,'RPPA median delta 0.161 [0.119, 0.192]; 2,000 shared patient draws. Source48 has no qualifying RPPA targets.')
f,axs=new();a,b,c,d,e,ff=axs
rp=read('real_panel_factorial_metrics.csv');platforms=['GPL570','GPL571','GeoMx1779','RPPA_gene_background'];names=['GPL570','GPL571','GeoMx1779','RPPA genes'];q=rp.query('cohort=="external138" and signature=="historical48" and member_state=="observed" and background=="measured_panel"').set_index('panel').loc[platforms];y=np.arange(4);a.barh(y,q.up_observed,color=BLUE,label="Up");a.barh(y,q.down_observed,left=q.up_observed,color=ORANGE,label="Down");a.legend(loc="lower right",fontsize=5.8);a.set_yticks(y,names);a.invert_yaxis();a.axvline(48,color=GRAY,lw=.6,ls='--');a.set(xlabel='Observed source members; up / down',title='Real measurement-list coverage')
mat=[]
for pn in platforms:
 mat.append([rp.query('cohort=="external138" and signature=="historical48" and panel==@pn and member_state==@me and background==@bg').rho.iloc[0] for me,bg in [('full','full_transcriptome'),('observed','full_transcriptome'),('full','measured_panel'),('observed','measured_panel')]])
mat=np.array(mat);cm=plt.get_cmap('Blues').copy();cm.set_bad('#ECEEF0');b.imshow(mat,vmin=0,vmax=1,aspect='auto',cmap=cm);b.set_yticks(y,names);b.set_xticks(range(4),['Full / RNA','Obs / RNA','Full / panel','Obs / panel'],rotation=30,ha='right',rotation_mode='anchor');b.tick_params(labelsize=5.6,length=0);b.set_title('Same-patient member x background test')
for i in range(4):
 for k in range(4):
  v=mat[i,k];b.text(k,i,'NA' if np.isnan(v) else f'{v:.3f}',fontsize=6,ha='center',va='center',color='white' if v>.65 else 'black')
S=read('frozen_transfer_summary_by_family.csv');methods=['percentile_rank','discovery_frozen_z','top1500_rank_area'];mechs=['uniform','low_expression','high_variance'];m=np.array([[S.query('signature_family=="Hallmark50" and cohort=="external138" and method==@md and mechanism==@mk and fraction==0.5',local_dict={'md':md,'mk':mk}).median_rho.iloc[0] for mk in mechs] for md in methods]);c.imshow(m,aspect='auto',vmin=0,vmax=1,cmap=cm);c.set_yticks(range(3),['Percentile rank','Frozen z','Truncated rank area']);c.set_xticks(range(3),['Uniform','Low expr.','High var.']);c.tick_params(length=0);c.set_title('Fixed 50% loss; median rank correlation')
for i in range(3):
 for k in range(3):c.text(k,i,f'{m[i,k]:.3f}',ha='center',va='center',fontsize=6.5,color='white')
P=read('signature_specific_summary.csv');cols=[GRAY,OLIVE,PURPLE,BLUE,ORANGE]
for mk,col,lab in zip(mechs+['up_only','down_only'],cols,['Uniform','Low expr.','High var.','Up only','Down only']):
 z=P.query('signature=="source48" and cohort=="external138" and method=="percentile_rank" and mechanism==@mk').sort_values('fraction');d.plot(z.fraction*100,z.median_rho,'o-',color=col,ms=3,lw=.8,label=lab)
d.set(xlabel='Removed members (% of full 48)',ylabel='Median rank correlation',title='Direction-specific removal differs',ylim=(0,1.04),xticks=[10,30,50]);d.legend(loc='lower left',ncol=2,fontsize=5.5);d.text(50,.05,'NA',ha='center',color=ORANGE,fontsize=5.5)
for ax,geo,n in [(e,'GSE17025',91),(ff,'GSE63678',7)]:
 z=read(geo+'_tumor_primary_scores.csv');xx=z.historical48_within_sample_gene_rank;yy=z.fixed_discovery29_within_sample_gene_rank;ax.scatter(xx,yy,s=10,color=BLUE);ax.set(xlabel='Measured full48 RNA projection',ylabel='Measured common29 RNA projection',title=f'{geo}: {n} tumors; RNA-only\n'+('47/48 vs 29/29 measured' if n==91 else '34/48 vs 24/29 measured'))
save(f,axs,4,'Lists are RNA counterfactuals, not independent assays. Complete-member/panel-background cells require full RNA.')
f,physical=new();axs=[physical[0],physical[2],physical[4],physical[1],physical[3],physical[5]]
for label,ax in zip('abcdef',axs):ax.set_gid(label);ax.texts[0].set_text(label)
a,b,c,d,e,ff=axs
cv=read('heldout_score_models.csv');order=['clinical','composition','generic','both'];a.boxplot([cv.query('model==@m').pooled_outer_R2 for m in order],widths=.5,showfliers=False);a.set_xticks(range(1,5),['Clinical','+ Comp.','+ Generic','+ Both']);a.tick_params(axis='x',labelsize=5.8);a.axhline(0,color=GRAY,lw=.5);a.set(ylabel='Held-out RNA score prediction R²',title='91 tumors; 20 repeated outer evaluations')
mix=read('admixture_curves.csv')
for j,(sam,g) in enumerate(mix.query('panel=="response"').groupby('sample')):b.plot(g.epithelial_RNA_fraction,g.score,color=[BLUE,ORANGE,GREEN,PURPLE,OLIVE][j],lw=.8,label=f'W{j+1}')
b.set(xlabel='Epithelial RNA fraction',ylabel='48-member rank score',title='Constructed mixtures from five donors');b.legend(ncol=3,fontsize=5.5,loc='upper left')
sc=read('lineage_pseudobulk.csv').query('dataset=="GSE173682" and correction=="corrected" and n_cells>=30').pivot(index='sample',columns='lineage',values='response')[['Epithelial','Stromal']].dropna()
assert len(sc)==5
for _,z in sc.iterrows():c.plot([0,1],z.values,'o-',color=BLUE,ms=3,lw=.7)
test=read('exact_donor_tests.csv').query('panel=="response" and contrast=="Epithelial minus Stromal"').iloc[0];c.set_xticks([0,1],['Epithelial','Stromal']);c.set(ylabel='Count-summed pseudobulk rank score',title=f'Five donors; exact P = {test.p:.4f}; BH10 q = {test.BH10:.3f}')
ex=read('external_survival_exploratory.csv');tc=read('TCGA_source48_PH_sensitivity.csv').query('model=="stage_subtype_stratified"').iloc[0];labs=[f'TCGA adj.; stage/subtype\nstrat.; n={int(tc.n)}, deaths={int(tc.events)}',f'CPTAC univ. full48\nn={int(ex.iloc[0].n)}, deaths={int(ex.iloc[0].events)}',f'CPTAC univ. common29\nn={int(ex.iloc[1].n)}, deaths={int(ex.iloc[1].events)}'];forest(d,labs,[tc.HR,*ex.HR],[tc.low,*ex.low],[tc.high,*ex.high],[BLUE,ORANGE,ORANGE],zero=1);d.set(xlabel='HR per cohort SD; 95% interval',title='Separate OS associations; no pooled HR',xlim=(.4,1.9));d.tick_params(axis='y',labelsize=5.6)
bt=read('source48_increment_bootstrap.csv');z=bt.pivot(index=['iteration','metric'],columns='model',values='oob');v=(z.plus_source48-z.clinical_molecular).xs('C_index',level='metric');e.hist(v,bins=30,color=BLUE,edgecolor='white',lw=.3);e.axvline(0,color=GRAY,lw=.7,ls='--');e.set(xlabel='OOB C-index: plus48 minus baseline',ylabel='Bootstrap fits',title='503 patients; 82 events; 500 refits')
nb=read('source48_decision_summary.csv');ff.errorbar(nb.threshold,nb.median_incremental_NB,yerr=[nb.median_incremental_NB-nb.stability_low,nb.stability_high-nb.median_incremental_NB],fmt='o-',color=BLUE,ms=3,capsize=2,lw=.8);ff.axhline(0,color=GRAY,lw=.5);ff.set(xlabel='Illustrative 36-month risk threshold',ylabel='Incremental OOB net benefit',title='Median and stability range; not CI',xticks=[.05,.1,.15,.2,.25,.3]);ff.tick_params(axis='x',labelsize=5.6)
f.text(.155,.977,'Score and composition descriptions',fontsize=8,fontweight='bold');f.text(.64,.977,'Separate exploratory clinical checks',fontsize=8,fontweight='bold')
save(f,axs,5,'These are interpretation checks, not validation of member restriction, death activity or clinical utility.')
print('Five figures rebuilt from retained source data.',flush=True)
