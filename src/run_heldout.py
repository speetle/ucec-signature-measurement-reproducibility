from pathlib import Path
import json
import numpy as np,pandas as pd,patsy
from scipy.stats import spearmanr
from sklearn.model_selection import KFold,GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
r=Path(__file__).resolve().parents[1];i=r/'input/discovery';out=r/'results/heldout';out.mkdir(parents=True,exist_ok=True);d=pd.read_csv(i/'CPTAC_model.csv',index_col=0);a=pd.read_csv(i/'RNA.csv.gz',index_col=0)[d.index];p=json.loads((i/'plan.json').read_text())['formal_source_response_panel'];y=a.rank(pct=True).loc[p['up']].mean()-a.rank(pct=True).loc[p['down']].mean();base='age10+C(grade)+C(histology)+C(subtype)+C(batch)';forms={'clinical':base,'composition':base+'+immune_logratio+stroma_logratio','generic':base+'+proliferation+inflammation+interferon+hypoxia+oxidative_stress+dissociation','both':base+'+immune_logratio+stroma_logratio+proliferation+inflammation+interferon+hypoxia+oxidative_stress+dissociation'};ids=patsy.dmatrix(forms['both'],d,return_type='dataframe').index;rows=[]
for rep in range(20):
 splits=list(KFold(5,shuffle=True,random_state=20261018+rep).split(ids))
 for model,formula in forms.items():
  X=patsy.dmatrix(formula,d.loc[ids],return_type='dataframe').values[:,1:];Y=y.loc[ids].values;pred=np.empty(len(ids))
  for tr,te in splits:
   search=GridSearchCV(make_pipeline(StandardScaler(),Ridge()),{'ridge__alpha':[.1,1,10,100]},cv=3,scoring='neg_mean_squared_error');search.fit(X[tr],Y[tr]);pred[te]=search.predict(X[te])
  rows.append({'repeat':rep,'model':model,'n':len(ids),'pooled_outer_R2':r2_score(Y,pred),'RMSE':np.sqrt(np.mean((Y-pred)**2))})
pd.DataFrame(rows).to_csv(out/'heldout_score_models.csv',index=False)
