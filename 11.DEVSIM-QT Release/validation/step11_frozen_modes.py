import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from dataclasses import replace
from config import Config
from modespace.hamiltonian import Grid
from modespace.junction_basis import JunctionCompatibleBuilder
from device.definition import DeviceDefinition
from transport.contact_model import solve_reference
from transport.negf_reference import solve_negf,blocks_density
from transport.bound_states import audit

def run(state,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);data=np.load(state);base=Config(**json.loads(str(data['config_json'])));grid=Grid.build(base);valleys=DeviceDefinition().frame.valleys()
    ref=solve_reference(base,grid,valleys);mu=ref.mu-base.potential_reference_v;rows=[];old=None
    for modes,sd_modes in [(4,16),(8,32),(12,48),(16,64),(24,100)]:
        c=replace(base,modes=modes,sd_modes=sd_modes);density=np.zeros(len(grid.xyz));dims=[];audits=[];tic=time.time()
        for v in valleys:
            b=JunctionCompatibleBuilder(c,grid,v);b.contacts=ref;b.block_only=True;s=b.build(data['phi_v']);dims.append([len(h) for h in s.onsite])
            blocks,flux,diag,parts=solve_negf(s,c,mu,mu,0);density+=blocks_density(s,blocks,2*v.multiplicity)
            audits.append(audit(s,c,mu,0))
        row={'requested_modes':modes,'requested_sd_modes':sd_modes,'actual_slice_dimensions':dims,'fixed_lead_mu_ev':mu,'seconds':time.time()-tic,'bound_audits':audits}
        if old is not None:
            x=grid.xyz[:,0]
            for name,mask in [('all',x>=0),('channel',(x>=c.gate_start)&(x<=c.gate_end)),('junction_extension',((x>=c.ls)&(x<c.gate_start))|((x>c.gate_end)&(x<=c.drain_start))),('reservoir',(x<c.ls)|(x>c.drain_start))]:
                row[name+'_density_relative_l2_change']=float(np.linalg.norm((old-density)[mask])/max(np.linalg.norm(density[mask]),1e-30))
        old=density;rows.append(row);np.savez_compressed(out/f'modes_{modes}_{sd_modes}.npz',density_cm3=density,xyz_nm=grid.xyz,mu_ev=mu)
        (out/'summary.json').write_text(json.dumps({'experiment':'A: fixed device phi and full-basis SP lead mu; not SCF observable convergence','rows':rows},indent=2));print(modes,sd_modes,row.get('all_density_relative_l2_change'),flush=True)
if __name__=='__main__':run(*sys.argv[1:])
