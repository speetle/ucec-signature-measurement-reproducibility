from pathlib import Path
import sys
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1];out=P/'figures';out.mkdir(exist_ok=True)
r=pd.read_csv(P/'results/GSE73024跨平台全部候选排名.csv');r=r[(r.scenario=='all92')&(r.case=='mean_drug_minus_control')&r.FGFR_positive];s=pd.read_csv(P/'results/GSE316734_AN3两个面板样本分数.csv');null=pd.read_csv(P/'results/GSE316734随机子面板5000次.csv')
assert len(r)==12 and len(s)==6 and len(null)==5000
plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
fig,axes=plt.subplots(1,3,figsize=(10.5,3.5));fig.subplots_adjust(left=.065,right=.985,bottom=.31,top=.86,wspace=.42)
ax=axes[0];names=['Infigratinib','Erdafitinib','Pemigatinib','Futibatinib'];colors=['#0072B2','#E69F00','#009E73','#CC79A7']
for name,c in zip(names,colors):
 z=r[r.drug==name].set_index('panel').loc[['source48','common29','full']];ax.plot([0,1,2],z['rank'],marker='o',ms=4,lw=1,color=c,label=name)
ax.axhline(10,color='0.45',ls='--',lw=.8);ax.invert_yaxis();ax.set_xticks([0,1,2],['48 genes','29 genes','All genes']);ax.set_ylabel('Rank among 92 conditions');ax.set_title('Cross-study class retrieval');ax.legend(fontsize=6,frameon=False,loc='upper left',bbox_to_anchor=(0,-.19),ncol=2,borderaxespad=0)
ax=axes[1]
for condition,color in [('DMSO','#0072B2'),('ATX968','#D55E00')]:
 d=s[s.condition==condition];ax.scatter(d.source48,d.common29,s=26,color=color,label=condition)
ax.set_xlabel('48-gene score');ax.set_ylabel('29-gene score');ax.set_title('Six samples; nested panels');ax.text(.04,.95,'Spearman ρ = 1.00',transform=ax.transAxes,va='top');ax.legend(frameon=False,fontsize=7,loc='lower right')
ax=axes[2];vals=null.spearman.to_numpy();bins=np.linspace(-1.05,1.05,15);ax.hist(vals,bins=bins,color='#8c9da8',edgecolor='white');ax.set_xlabel('Spearman ρ versus 48-gene score');ax.set_ylabel('Number of random subsets');ax.set_title('Post hoc subset reference');ax.text(.04,.95,'23.68% identical ordering',transform=ax.transAxes,va='top',fontsize=7)
for ax,label in zip(axes,'abc'):ax.annotate(label,xy=(0,1),xycoords='axes fraction',xytext=(-16,9),textcoords='offset points',fontweight='bold',fontsize=10)
fig.savefig(out/'独立扰动整合图.pdf');fig.savefig(out/'独立扰动整合图.png',dpi=300);print('saved')
