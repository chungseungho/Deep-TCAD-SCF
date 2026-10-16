"""Controlled frozen-field bias/mode/partition sensitivity, not biased SCF."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import replace
import numpy as np
from scipy.linalg import block_diag
from config import Config,VALLEYS
from modespace.hamiltonian import Grid
from modespace.junction_basis import JunctionCompatibleBuilder
from multispace.occupation import solve_states,occupations
from multispace.constraints import number_bounds
from multispace.tracking import match_states
from device.regions import mode_regions
from solver import write_json,atomic_npz
from provenance import source_digest


def run(source,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);saved=np.load(source)
    base=Config(**json.loads(str(saved['config_json'])));g=Grid.build(base);mu=float(saved['mu_l_ev']);rows=[]
    for m in (2,4,10):
        c=replace(base,modes=m,sd_modes=2*m);builders=[JunctionCompatibleBuilder(c,g,v) for v in VALLEYS];last=None
        for vd in (0.,.001,.002,.005,.01,.05,.1):
            phi=saved['phi_v']+vd*g.xyz[:,0]/c.length
            states=[];systems=[];details=[]
            for j,b in enumerate(builders):
                s=b.build(phi);labels=mode_regions(s,c);alts=[mode_regions(s,c,sign*c.dx) for sign in (-1,1)]
                state=occupations(solve_states(s.H,labels,c,alts),mu,mu-vd,c.temperature,strict=False)
                active=state['active'];ids=np.flatnonzero(state['ambiguous']);profiles=np.array([[sum(abs(state['u'][s.sl(i),a])**2) for i in range(len(g.x))] for a in ids])
                d={'valley':s.valley.name,'active':int(active.sum()),'ambiguous_active':int(state['ambiguous'].sum()),
                    'boundary_unstable_active':int(sum(state['boundary_unstable']&active)),
                    'occupation_uncertainty_electrons':float(4*state['occupation_uncertainty_per_state'].sum())}
                if last is not None:
                    match=match_states(last[j],state,c.ms_degenerate_ev);a,z=match['old_indices'],match['new_indices'];keep=last[j]['active'][a]|state['active'][z]
                    d.update(matched_active_pairs=int(keep.sum()),matched_origin_changes=int(sum(last[j]['origin'][a][keep]!=state['origin'][z][keep])),
                             median_active_overlap=float(np.median(match['overlap_squared'][keep])) if keep.any() else 1.)
                atomic_npz(out/f'M{m}_vd{vd:g}_v{j}.npz',energies_ev=state['e'],weights=state['weights'],origin=state['origin'],
                    active=active,ambiguous=state['ambiguous'],origin_contrast=state['contrast'],reservoir_reach=state['reach'],
                    active_ambiguous_ids=ids,longitudinal_probability=profiles,x_nm=g.x,mu_l_ev=mu,mu_r_ev=mu-vd)
                states.append(state);details.append(d)
            lo,hi=number_bounds(states,[4]*3,mu,vd,c.temperature)
            rows.append({'channel_modes':m,'sd_modes':2*m,'vd_v':vd,'mu_l_ev':mu,'number_interval':[lo,hi],'valleys':details})
            last=states
        print('origin sensitivity mode',m,flush=True)
    write_json(out/'summary.json',{'source_sha256':source_digest(),'cases':rows,
        'scope':'phi(vd)=saved equilibrium phi + vd*x/L controlled frozen family; NOT self-consistent bias continuation',
        'unknown_occupations_applied_to_density':False,'source_state':str(source)})

if __name__=='__main__':run(sys.argv[1],sys.argv[2])
