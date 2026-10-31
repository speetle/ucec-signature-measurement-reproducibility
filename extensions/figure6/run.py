from pathlib import Path
import subprocess,sys,hashlib,json
p=Path(__file__).resolve().parent
assets=json.loads((p/'MANIFEST.json').read_text())
for name, info in assets.items():
 f=p/'assets'/name
 if not f.is_file() or hashlib.sha256(f.read_bytes()).hexdigest()!=info['sha256']:raise ValueError('Frozen asset mismatch: '+name)
for script in ['prepare_input.py','scripts/分析GSE316734随机子面板.py','scripts/GSE73024跨平台检索.py','scripts/draw_figure6.py']:
 subprocess.run([sys.executable,str(p/script)],check=True)
a=json.loads((p/'results/GSE316734_AN3覆盖与扰动摘要.json').read_text());b=json.loads((p/'results/GSE316734随机子面板摘要.json').read_text());assert a['rank_spearman']==1 and a['discordant_pairs']==0 and abs(b['proportion_perfect_order']-.2368)<1e-12
import pandas as pd
r=pd.read_csv(p/'results/GSE73024跨平台全部候选排名.csv');r=r[(r.scenario=='all92')&(r.case=='mean_drug_minus_control')&r.FGFR_positive];v={k:int(r[r.panel==k]['rank'].min()) for k in ['source48','common29','full']};assert v=={'source48':11,'common29':14,'full':19}
print('PASS: Figure 6 data assertions',v,'; rho=1; random perfect order=0.2368')
