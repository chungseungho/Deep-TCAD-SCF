"""Append explicit false validation flags for failed actual bias attempts."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[2]
status=json.loads((ROOT/'STATUS.json').read_text());failed=[]
for record in status['actual_run_records']:
 if record['exit_code']==0:continue
 cmd=record['command'];get=lambda key:cmd[cmd.index(key)+1] if key in cmd else None
 log=ROOT/'results'/record['run']/'run/execution.log';text=log.read_text()
 trial=ROOT/'results'/record['run']
 bound=any('BOUND_PREPARATION_REQUIRED' in p.read_text() for p in trial.glob('*newton_trials.jsonl'))
 failed.append({'run':record['run'],'vg_v':float(get('--vg')),'vd_v':float(get('--vd')),'id_a':None,
  'SCF_CONVERGED':False,'TRANSPORT_CONSERVED':False,'NUMERICAL_CONVERGENCE_VERIFIED':False,
  'PHYSICAL_MODEL_VALID':False,'QUANTITATIVE_VALIDATED':False,'TRAINING_READY':False,
  'status':'BOUND_PREPARATION_REJECTIONS_AND_TIMEOUT' if bound else 'TIMEOUT' if record['exit_code']==124 else 'RUN_FAILED',
  'scope':'No final physical state was accepted. Timeout is not proof that no solution exists.',
  'execution_record':f"results/{record['run']}/run/execution.json"})
status['failed_bias_attempts']=failed
(ROOT/'STATUS.json').write_text(json.dumps(status,indent=2)+'\n')
print('failed actual attempts recorded',len(failed))
