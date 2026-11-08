from pathlib import Path
import argparse,subprocess,sys,shutil
p=argparse.ArgumentParser();p.add_argument('project',type=Path);p.add_argument('--work',type=Path,required=True,help='新的空工作目录');p.add_argument('--release',type=Path);p.add_argument('--cervical',type=Path);a=p.parse_args();project=a.project.expanduser().resolve();work=a.work.expanduser().resolve();here=Path(__file__).resolve().parent
if work.exists() and any(work.iterdir()):raise SystemExit('工作目录非空，为避免覆盖请指定新的空目录')
work.mkdir(parents=True,exist_ok=True);inp=work/'input';inp.mkdir();
def run(name,*args):subprocess.run([sys.executable,str(here/name),*map(str,args)],check=True)
release=(a.release.expanduser().resolve() if a.release else project/'UCEC_v15_逐条实质修订_20261020/public_release');opts=['--release',release];cerv=['--cervical',a.cervical.expanduser().resolve()] if a.cervical else []
run('重建共同成员与背景.py',project,'--output',inp,*opts,*cerv);run('重建匹配特征和模块.py',project,'--output',inp,*opts)
for cohort in ['discovery','external']:
 target=inp/cohort;target.mkdir();shutil.copy2(release/'input'/cohort/'protein.csv.gz',target/'protein.csv.gz')
for n in ['生成参照.py','汇总与审计.py']:shutil.copy2(here/n,work/n)
subprocess.run([sys.executable,str(work/'生成参照.py')],check=True);subprocess.run([sys.executable,str(work/'汇总与审计.py')],check=True)
subprocess.run(['Rscript',str(here/'三队列配对统计.R'),str(inp),str(work/'output')],check=True)
print('完成：',work/'output',flush=True)
