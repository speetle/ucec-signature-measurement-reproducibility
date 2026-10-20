from pathlib import Path
import json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];F=R/'frozen';O=R/'results';rows=[]
def compare(folder,name,reference=None,columns=None):
 a=pd.read_csv(O/folder/name,float_precision='round_trip');b=pd.read_csv(F/(reference or name),float_precision='round_trip')
 if columns:a=a[columns];b=b[columns]
 assert a.shape==b.shape,(folder,name,a.shape,b.shape)
 errors=[];textgood=True
 for c in a:
  if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c]):
   x=a[c].to_numpy(dtype=float);y=b[c].to_numpy(dtype=float);assert np.array_equal(np.isnan(x),np.isnan(y)),(name,c,'missing');valid=np.isfinite(x)&np.isfinite(y);errors.append(float(np.max(np.abs(x[valid]-y[valid]))) if valid.any() else 0.)
  else:textgood &= a[c].fillna('<NA>').astype(str).equals(b[c].fillna('<NA>').astype(str))
 error=max(errors or [0]);record=dict(module=folder,file=name,rows=len(a),columns=len(a.columns),max_absolute_error=error,labels_and_missingness_equal=bool(textgood),passed=bool(error<1e-10 and textgood));rows.append(record)
for module,names in {
 'null':['weight_contrasts_bootstrap.csv','within_batch_permutations.csv','null_summary.csv'],
 'controls':['coexpression_null_results.csv','coexpression_null_summary.csv','protein29_jackknife.csv','LOXL3_assumption_stress.csv'],
 'masks':['real_panel_factorial_metrics.csv','frozen_transfer_scenarios.csv.gz','frozen_transfer_summary_by_family.csv','signature_specific_summary.csv'],
 'arrays':['GSE17025_tumor_primary_scores.csv','GSE63678_tumor_primary_scores.csv','array_panel_sensitivity.csv'],
 'rppa':['RPPA_Hallmark_dependent_correlations.csv','RPPA_Hallmark_coverage.csv'],
 'heldout':['heldout_score_models.csv'],
 'clinical':['external_survival_exploratory.csv'],
 'precision':['precision_background_results.csv','reduced_model_sensitivity.csv']}.items():
 for name in names:compare(module,name)
compare('cells','lineage_pseudobulk.csv',columns=['sample','dataset','lineage','n_cells','correction','legacy','response']);compare('cells','admixture_curves.csv')
check=json.loads((O/'review/review_closure_validation.json').read_text());summary=dict(comparisons=rows,all_passed=all(x['passed'] for x in rows) and check['all_passed'],tolerance=1e-10,scope='fresh recomputation from processed RNA/protein, antibody measurements, mapped array gene matrices, lineage count aggregates, and clinical analysis inputs; not raw-cell reprocessing or raw-assay re-quantification')
(O/'extended_validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));assert summary['all_passed']
