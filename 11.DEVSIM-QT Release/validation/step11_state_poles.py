import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config import Config
from modespace.hamiltonian import Grid
from modespace.junction_basis import JunctionCompatibleBuilder
from device.definition import DeviceDefinition
from transport.contact_model import solve_reference
from transport.bound_states import audit
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'results/state_poles';out.mkdir(exist_ok=True);rows=[]
for name in ('raised_diagnostic','raised_eq_final','raised_refined_eq','raised_fine_eq','buffer_eq'):
 s=np.load(ROOT/f'results/{name}/states/vg0.400000_vd0.000000.npz');c=Config(**json.loads(str(s['config_json'])));g=Grid.build(c);valleys=DeviceDefinition().frame.valleys();r=solve_reference(c,g,valleys);audits=[]
 for v in valleys:
  b=JunctionCompatibleBuilder(c,g,v);b.block_only=True;b.contacts=r;system=b.build(s['phi_v']);audits.append(audit(system,c,r.mu-c.potential_reference_v,0.))
 row={'state':name,'modes':c.modes,'sd_modes':c.sd_modes,'audits':audits};rows.append(row);(out/'summary.json').write_text(json.dumps(rows,indent=2));print(name,[(a['valley'],[p['energy_ev'] for p in a['bound_poles_below_both_leads']]) for a in audits],flush=True)
