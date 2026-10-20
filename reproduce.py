#!/usr/bin/env python3
"""Processed-input reproduction. Public-source retrieval is explicit and separately logged."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,shutil,platform,time,os,datetime
R=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--figures-only',action='store_true');p.add_argument('--verify-only',action='store_true');p.add_argument('--fetch-missing',action='store_true');a=p.parse_args()
mode='verify_only' if a.verify_only else ('figures_only' if a.figures_only else 'full_analysis')
run_id=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
run_dir=R/'results/runs'/run_id;log=run_dir/'logs';log.mkdir(parents=True,exist_ok=True)
try:commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,stderr=subprocess.DEVNULL,text=True).strip()
except Exception:commit=None
try:dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=R,stderr=subprocess.DEVNULL,text=True).strip())
except Exception:dirty=None
source_hashes={str(x.relative_to(R)):hashlib.sha256(x.read_bytes()).hexdigest() for x in [R/'reproduce.py',R/'requirements.txt',*sorted((R/'src').glob('*.py')),*sorted((R/'src').glob('*.R'))]}
record={'source_tree_dirty':dirty,'analysis_source_sha256':source_hashes,'run_id':run_id,'release':'v1.0.1','mode':mode,'commit':commit,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'running','stages':[],'input_manifest_verified':False}
def persist():
 (run_dir/'run_record.json').write_text(json.dumps(record,indent=2))
 status='stage_status.json' if mode=='full_analysis' else mode+'_stage_status.json'
 (R/'results'/status).write_text(json.dumps(record['stages'],indent=2))
persist()
env=os.environ.copy()
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:env[key]='2'
stages=record['stages']
def run(name,cmd):
 start=time.time()
 with (log/(name+'.log')).open('w') as f:result=subprocess.run(cmd,cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT)
 rec=dict(stage=name,exit_code=result.returncode,seconds=round(time.time()-start,2));stages.append(rec);print(json.dumps(rec),flush=True)
 record['status']='failed' if result.returncode else 'running';persist()
 if result.returncode:raise SystemExit('STOP. See '+str(log/(name+'.log')))
if a.fetch_missing:
 if any(not (R/'input/external'/n).is_file() for n in ['RNA.csv.gz','protein.csv.gz','metadata.csv']):run('fetch_external',[sys.executable,str(R/'src/prepare_external_inputs.py')])
 if any(not (R/x['file']).is_file() for x in json.loads((R/'optional_sources.json').read_text())):run('fetch_rppa',[sys.executable,str(R/'src/fetch_optional_inputs.py')])
items=json.loads((R/'input_manifest.json').read_text());bad=[]
for item in items:
 f=R/item['file']
 if not f.is_file() or hashlib.sha256(f.read_bytes()).hexdigest()!=item['sha256']:bad.append(item['file'])
if bad:raise SystemExit('Input verification failed. Exact versions required; use --fetch-missing for source-retrieved inputs: '+str(bad))
record['input_manifest_verified']=True;record['input_count']=len(items);persist()
print('Manifest input verification PASS',len(items),flush=True)
if a.verify_only:
 record['status']='verified_inputs';persist();sys.exit(0)
if not a.figures_only:
 for name,script in [('primary','recompute_primary'),('patient_coverage','run_patient_coverage'),('actual_ms','run_ms_benchmarks'),('ms_verification','verify_ms_points'),('weight_and_mismatch','run_null'),('structure','run_coexpression_null'),('jackknife','run_jackknife'),('missing_assumptions','run_missingness_bounds'),('measurement_masks','run_masks'),('arrays','run_arrays'),('rppa','run_rppa'),('heldout','run_heldout'),('cells','run_cells'),('precision','run_precision'),('reduced_model','run_reduced')]:run(name,[sys.executable,str(R/'src'/f'{script}.py')])
 if True:
  rs=shutil.which('Rscript')
  if not rs:raise SystemExit('Rscript and survival required for clinical analyses')
  out=R/'results/clinical';out.mkdir(exist_ok=True);shutil.copy2(R/'input/clinical/TCGA48_input.csv',out/'TCGA48_input.csv');run('clinical',[rs,str(R/'src/clinical_increment.R'),str(out)]);run('associations',[rs,str(R/'src/associations.R'),str(R)])
 run('review_audits',[sys.executable,str(R/'src/review_closure.py')])
 if True:
  run('original_numeric_validation',[sys.executable,str(R/'src/validate_outputs.py')]);run('extended_validation',[sys.executable,str(R/'src/validate_extended.py')])
run('figures',[sys.executable,str(R/'src/build_figures.py')])
(R/'results/run_environment.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'figures_only':a.figures_only,'source_fetch_requested':a.fetch_missing,'raw_assay_reprocessing':False,'scope':'All displayed panel statistical computations from explicit processed matrices, phenotype designs, platform definitions, lineage count aggregates and clinical analysis inputs; upstream cell QC and original assay processing are separate.'},indent=2))
record['status']='completed';record['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();persist()
print('Completed. Run record: '+str(run_dir/'run_record.json'),flush=True)
