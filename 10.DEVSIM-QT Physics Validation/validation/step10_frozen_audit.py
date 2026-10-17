"""Mode/contact/pole audit at a saved potential; never claims a new SCF state."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import replace
import numpy as np
from config import Config
from device.definition import Frame
from modespace.hamiltonian import Grid,ModeBuilder
from modespace.junction_basis import JunctionCompatibleBuilder
from transport.registry import create_transport_backend
from transport.negf_reference import device_blocks,neutral_lead_mu
from transport.bound_states import audit


def run(state,out):
    state=np.load(state);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    base=Config(**json.loads(str(state['config_json'])));grid=Grid.build(base);phi=state['phi_v'];rows=[];densities=[]
    for m,sd in ((2,4),(4,8),(8,16),(12,24)):
        c=replace(base,modes=m,sd_modes=sd)
        cls=JunctionCompatibleBuilder if c.junction_basis=='nested' else ModeBuilder
        builders=[cls(c,grid,v) for v in Frame(c.transport_axis).valleys()]
        mu,contact=neutral_lead_mu(builders,c,grid)
        backend=create_transport_backend('negf',c)
        q=backend.solve(phi,0.,{'builders':builders,'grid':grid,'fixed_mu':mu,'bias':0.,'detail':False})
        densities.append(q['n']);pole_reports=[];mismatch=[]
        for s in q['systems']:
            pole_reports.append(audit(s,c,mu,0.))
            onsite,hop,t,leads=device_blocks(s,c,(c.potential_reference_v,c.potential_reference_v))
            mismatch.append(max(float(np.linalg.norm(onsite[i]-(u*e)@u.T,2))
                for i,(e,u) in zip((0,-1),leads)))
        row={'channel_modes':m,'requested_sd_modes':sd,'dimensions':[s.dimension for s in q['systems']],
            'mu_ev':mu,'electrons':float(q['n'].sum()*grid.cell_volume*1e-21),
            'contact_neutrality_relative_error':contact['neutrality_relative_error'],
            'endpoint_to_assumed_lead_onsite_norm_ev':mismatch,'bound_poles':pole_reports,
            'scope':'frozen potential; recalibrated reservoir mu; not a new SCF solution'}
        rows.append(row)
        np.savez_compressed(out/f'm{m}_sd{sd}.npz',xyz_nm=grid.xyz,phi_v=phi,density_cm3=q['n'],mu_ev=mu)
        print(m,sd,row['electrons'],mu,flush=True)
    reference=densities[-1]
    for i,row in enumerate(rows):
        row['density_relative_l2_to_largest']=float(np.linalg.norm(densities[i]-reference)/np.linalg.norm(reference))
        row['density_relative_l2_to_previous']=None if i==0 else float(np.linalg.norm(densities[i]-densities[i-1])/np.linalg.norm(densities[i]))
    summary={'cases':rows,'thresholds_predeclared':{'density_relative_l2':.01,'current_relative':.01,'chemical_potential_ev':1e-4},
        'mode_plateau_observed_at_frozen_potential':bool(rows[-1]['density_relative_l2_to_previous']<.01 and abs(rows[-1]['mu_ev']-rows[-2]['mu_ev'])<1e-4),
        'physical_vs_projection_bound_state_classification':'UNRESOLVED: retained-model poles must survive full mode/contact/grid convergence',
        'device_numerical_convergence_verified':False}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary

if __name__=='__main__':run(sys.argv[1],sys.argv[2])
