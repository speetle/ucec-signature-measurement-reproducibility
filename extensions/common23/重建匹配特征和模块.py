from pathlib import Path
import argparse,numpy as np,pandas as pd
from scipy.stats import rankdata
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits
ap=argparse.ArgumentParser();ap.add_argument('project',type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--release',type=Path);args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
i=(args.release if args.release else args.project/'UCEC_v15_逐条实质修订_20261020/public_release')/'input/discovery';D=pd.read_csv(i/'CPTAC_model.csv',index_col=0);R=pd.read_csv(i/'RNA.csv.gz',index_col=0)[D.index];P=pd.read_csv(i/'protein.csv.gz',index_col=0)[D.index]
genes=R.index[R.notna().all(axis=1)&(R.std(axis=1)>0)];det=P.notna().mean(axis=1).reindex(genes).fillna(0);features=pd.DataFrame({'mean':R.loc[genes].mean(axis=1),'log_sd':np.log(R.loc[genes].std(axis=1)),'RNA_detection':(R.loc[genes]>0).mean(axis=1),'protein_detection':det});rr=R.loc[genes].rank(axis=1).values;rr-=rr.mean(axis=1,keepdims=True);rr/=np.sqrt((rr*rr).sum(axis=1,keepdims=True))
for c in ['proliferation','inflammation','interferon','hypoxia','oxidative_stress','dissociation']:
 a=rankdata(D[c].values);a-=a.mean();a/=np.linalg.norm(a);features[c]=rr@a
features.to_csv(out/'discovery_matching_features.csv');threadpool_limits(2);embed=PCA(n_components=12,random_state=20261019).fit_transform(rr);mods=KMeans(n_clusters=12,random_state=20261019,n_init=10).fit_predict(embed);pd.DataFrame({'module':mods},index=genes).to_csv(out/'discovery_coexpression_modules.csv');print('重建完成',len(genes),flush=True)
