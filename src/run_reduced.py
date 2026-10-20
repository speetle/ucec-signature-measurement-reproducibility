from pathlib import Path
import numpy as np,pandas as pd,patsy
from scipy.stats import rankdata
r=Path(__file__).resolve().parents[1];out=r/'results/precision';out.mkdir(parents=True,exist_ok=True);m=pd.read_csv(r/'input/external/metadata.csv',index_col=0);s=pd.read_csv(r/'frozen/external138_patient_scores.csv',index_col=0);m['age10']=pd.to_numeric(m.age,errors='coerce')/10;m['grade_binary']=m.Histologic_grade.map(lambda x:'low_intermediate' if str(x).startswith(('G1','G2')) else ('high' if str(x).startswith('G3') else np.nan));m['hist_binary']=np.where(m.Histologic_Type=='Endometrioid carcinoma','endometrioid','other');m['subtype']=m.Genomic_subtype;m['batch']=m.Batch;rows=[]
def corr(Y,X):
 z=Y-X@np.linalg.lstsq(X,Y,rcond=None)[0];z=rankdata(z,axis=0);z-=z.mean(axis=0);z/=np.sqrt((z*z).sum(axis=0));return z[:,:3].T@z[:,3]
for name,f in [('minimal','age10+C(subtype)+C(batch)'),('binary_grade_histology','age10+C(subtype)+C(batch)+C(grade_binary)+C(hist_binary)')]:
 X=patsy.dmatrix(f,m,return_type='dataframe');Y=s.loc[X.index].values;point=corr(Y,X.values);rng=np.random.default_rng(2026101802);b=[];rankdef=0
 for k in range(10000):
  ix=rng.integers(0,len(X),len(X));xx=X.values[ix];rankdef+=int(np.linalg.matrix_rank(xx)<xx.shape[1]);c=corr(Y[ix],xx);b.append(c[1]-c[0])
 lo,hi=np.quantile(b,[.025,.975]);rows.append({'model':name,'n':len(X),'columns':X.shape[1],'complete_rank_rho':point[0],'common_rank_rho':point[1],'common_z_rho':point[2],'membership_delta':point[1]-point[0],'CI_low':lo,'CI_high':hi,'draws':10000,'rank_deficient_draws_retained':rankdef})
pd.DataFrame(rows).to_csv(out/'reduced_model_sensitivity.csv',index=False)
missing=m[['age10','subtype','Histologic_Type','Histologic_grade']].isna();missing[missing.any(axis=1)].to_csv(out/'external_adjustment_exclusions.csv');print(pd.DataFrame(rows).to_string(index=False))
