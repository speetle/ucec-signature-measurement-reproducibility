from controls_common import *
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
rng=np.random.default_rng(20261019);loo=[];missing=[];missingassoc=[];strat=[];adjust=[];sings=[]
for co,(r,p,m) in cohorts.items():
 q=r.rank(pct=True);z=p.sub(p.mean(axis=1),axis=0).div(p.std(axis=1,ddof=1),axis=0);a=signed(q,u,d);b=signed(q,uo,do);t=signed(z,uo,do);base=rho(b,t)-rho(a,t)
 for g in obs:
  ou=[x for x in uo if x!=g];od=[x for x in do if x!=g];tt=signed(z,ou,od);bb=signed(q,ou,od)
  loo.append(dict(cohort=co,gene=g,original_delta=base,endpoint_only_delta=rho(b,tt)-rho(a,tt),joint_member_deletion_delta=rho(bb,tt)-rho(a,tt)))

pd.DataFrame(loo).to_csv(O/'protein29_jackknife.csv',index=False)
