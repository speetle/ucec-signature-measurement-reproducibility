from pathlib import Path
import json,hashlib,numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];records=[]
for stage in ['ms','clinical']:
 for p in (R/'results'/stage).glob('*.csv*'):
  e=R/'expected'/stage/p.name
  if not e.exists():continue
  a=pd.read_csv(p,float_precision='round_trip');b=pd.read_csv(e,float_precision='round_trip');cols=[c for c in a.columns if c in b and pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c]) and a[c].dtype!='bool'];same=a.shape==b.shape;error=0.
  if same:
   for c in cols:
    x=a[c].to_numpy(float);y=b[c].to_numpy(float)
    if not np.array_equal(np.isnan(x),np.isnan(y)):same=False;break
    finite=np.isfinite(x)&np.isfinite(y)
    if finite.any():error=max(error,float(np.max(abs(x[finite]-y[finite]))))
   strings=[c for c in a.columns if c not in cols];same=same and all(a[c].fillna('<NA>').astype(str).tolist()==b[c].fillna('<NA>').astype(str).tolist() for c in strings)
  records.append({'module':stage,'file':p.name,'shape_equal':a.shape==b.shape,'labels_and_missingness_equal':bool(same),'max_numeric_error':error,'passed':bool(same and error<1e-10)})
primary=json.loads((R/'results/primary/validation.json').read_text());rec={'comparisons':records,'primary':primary,'all_passed':all(z['passed'] for z in records) and primary['passed'],'scope':'Primary/MS and clinical tables in this check; extended_validation.json separately covers recalculated controls, masks, arrays, RPPA, heldout and processed-lineage statistics. Original raw assay and cell QC are outside this processed-input scope.'};(R/'results/numerical_validation.json').write_text(json.dumps(rec,indent=2));print(json.dumps({'tables_compared':len(records),'passed':rec['all_passed'],'max_numeric_error':max(z['max_numeric_error'] for z in records)},indent=2));assert rec['all_passed'],[z for z in records if not z['passed']]
