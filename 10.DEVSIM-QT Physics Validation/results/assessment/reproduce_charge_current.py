"""Check saved SCF charge residuals and full-device bond-current endpoints."""
from pathlib import Path
import json,numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).parent
rows=[]
for folder in sorted((ROOT/'results').glob('final_*')):
 for p in sorted((folder/'states').glob('*.npz')):
  a=np.load(p);metadata=json.loads(p.with_suffix('.json').read_text())
  nq=float(a['electron_density_cm3'].sum()*a['cell_volume_nm3']*1e-21)
  np_=float(np.dot(a['poisson_density_cm3'],a['node_volume_nm3'])*1e-21)
  b=a['bond_current_a'];current=float(a['id_a']);scale=max(abs(current),1e-15)
  mapping=json.loads((folder/'mapping.json').read_text())
  rows.append({'run':folder.name,'quantum_electrons':nq,'saved_poisson_electrons':np_,
    'saved_charge_difference_relative':abs(nq-np_)/max(abs(nq),1e-30),
    'charge_difference_scope':'finite SCF residual between saved input charge and reevaluated quantum output; not a mapping-operator defect',
    'map_operator_random_field_charge_error':mapping['charge_relative_error'],
    'map_cell_marginal_error':mapping.get('cell_marginal_error'),
    'source_plus_drain_bond_current_a':float(b[0]-b[-1]),
    'all_bond_current_relative_defect':float(max(abs(b-current))/scale),
    'gauss_relative_error':metadata['gauss_relative_error']})
(OUT/'charge_current_checks.json').write_text(json.dumps({'states':rows,'training_ready':False},indent=2)+'\n')
print('saved state charge/current checks',len(rows))
