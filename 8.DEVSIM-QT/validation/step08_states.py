"""Audit every recorded active ambiguous state without changing its origin."""
import sys,json,csv
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config import Config,VALLEYS
from modespace.hamiltonian import Grid
from modespace.junction_basis import JunctionCompatibleBuilder
from multispace.tracking import match_states
from solver import write_json,atomic_npz
from provenance import source_digest


def run(folder,out):
    folder=Path(folder);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    c=Config(**json.loads((folder/'settings.json').read_text())['config']);grid=Grid.build(c)
    rows=[];tracking=[]
    for valley,v in enumerate(VALLEYS):
        system=JunctionCompatibleBuilder(c,grid,v).build(np.zeros(len(grid.xyz)));last=None;lastvd=None
        for f in sorted((folder/'failures').glob(f'*_valley{valley}_occupation.npz')):
            s=dict(np.load(f));vd=float(f.stem.split('_vd')[1].split('_')[0]);idx=np.flatnonzero(s['ambiguous']);profiles=[]
            matched=match_states(last,s,c.ms_degenerate_ev) if last is not None else None
            if matched is not None:
                i,j=matched['old_indices'],matched['new_indices'];active=last['active'][i]|s['active'][j]
                tracking.append({'valley':v.name,'from_vd':lastvd,'to_vd':vd,
                    'active_pairs':int(active.sum()),'matched_origin_changes':int(sum(last['origin'][i][active]!=s['origin'][j][active])),
                    'active_overlap_min':float(min(matched['overlap_squared'][active])),
                    'does_not_determine_injection_origin':True})
            for j in idx:
                prob=np.array([sum(abs(s['u'][system.sl(i),j])**2) for i in range(len(grid.x))]);profiles.append(prob)
                nearest=float(min(abs(np.delete(s['e'],j)-s['e'][j])))
                rows.append({'vd_v':vd,'valley':v.name,'state_index':int(j),'energy_ev':float(s['e'][j]),
                    'wL':float(s['weights'][j,0]),'wC':float(s['weights'][j,1]),'wR':float(s['weights'][j,2]),
                    'occupation_difference':float(s['occupation_uncertainty_per_state'][j]),
                    'reach':float(s['reach'][j]),'contrast':float(s['contrast'][j]),
                    'boundary_unstable':bool(s['boundary_unstable'][j]),'cluster_unstable':bool(s['cluster_unstable'][j]),
                    'nearest_energy_gap_ev':nearest,'longitudinal_ipr':float(sum(prob**2)),
                    'center_of_probability_nm':float(prob@grid.x),'largest_region':int(s['largest_region'][j]),
                    'origin':int(s['origin'][j])})
            atomic_npz(out/f'{f.stem}_profiles.npz',x_nm=grid.x,state_indices=idx,
                       longitudinal_probability=np.asarray(profiles),energy_ev=s['e'][idx])
            last,lastvd=s,vd
    with (out/'active_states.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write_json(out/'summary.json',{'source_sha256':source_digest(),'active_ambiguous_state_records':len(rows),
        'tracking':tracking,'scope':'same fixed mode basis, biased failed-iterate profiles; not converged biased wavefunctions',
        'classification_and_tracking_are_independent':True})
    print('state audit',len(rows),flush=True)

if __name__=='__main__':run(sys.argv[1],sys.argv[2])
