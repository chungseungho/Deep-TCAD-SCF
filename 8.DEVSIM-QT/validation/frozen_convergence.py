"""Frozen-field GAA subspace study; not a self-consistent plateau claim."""
import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import replace
import numpy as np
from config import Config,VALLEYS
from modespace.hamiltonian import Grid,ModeBuilder,density_from_matrix
from modespace.junction_basis import JunctionCompatibleBuilder
from device.regions import mode_regions
from multispace.occupation import solve_states,occupations
from transport.ms_postprocess import transmission
from solver import write_json,atomic_npz
from provenance import source_digest


def run(source,out):
    s=np.load(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    base=Config(**json.loads(str(s['config_json'])));grid=Grid.build(base)
    # Fixed μ isolates the retained-subspace effect from finite-number μ recalibration.
    mu=float(s['mu_l_ev']);energies=np.linspace(mu-.1,mu+.3,161);rows=[];datasets={}
    settings=[('nested_2_4','fixed','nested',2,4),('nested_4_8','fixed','nested',4,8),
              ('nested_6_12','fixed','nested',6,12),('plain_2_4','fixed','plain',2,4),
              ('local_nested_2_4','local','nested',2,4)]
    for name,basis,junction,m,sd in settings:
        c=replace(base,basis=basis,junction_basis=junction,modes=m,sd_modes=sd)
        cls=JunctionCompatibleBuilder if junction=='nested' else ModeBuilder
        n=np.zeros(len(grid.xyz));T=np.zeros(len(energies));ortherr=0.;missing=0.;ranks=[];tail=0.;tic=time.time()
        for v in VALLEYS:
            system=cls(c,grid,v).build(s['phi_v']);g=2*v.multiplicity
            state=occupations(solve_states(system.H,mode_regions(system,c),c),mu,mu,c.temperature)
            rho=(state['u']*state['occupation'])@state['u'].T
            n+=density_from_matrix(system,rho,g)
            T+=g*transmission(system,c,energies,0)
            for U,mask in zip(system.bases,grid.mask):
                u=U[mask];ortherr=max(ortherr,float(np.max(abs(u.T@u-np.eye(u.shape[1])))))
            i=np.where(grid.x<c.ls)[0][-1]
            a=system.bases[i].reshape(-1,system.bases[i].shape[-1]);b=system.bases[i+1].reshape(-1,system.bases[i+1].shape[-1])
            missing=max(missing,float(np.linalg.norm(np.eye(b.shape[1])-(a.T@b).T@(a.T@b))))
            # Population above μ+8 kBT is a measured tail, not all added-mode charge.
            tail+=g*state['occupation'][state['e']>mu+8*8.617333262145e-5*c.temperature].sum()
            ranks.append(system.dimension)
        data={'density_cm3':n,'energies_ev':energies,'weighted_transmission':T,'source_sha256':source_digest(),
              'config_json':json.dumps(c.to_dict()),'mu_l_ev':mu,'phi_v':s['phi_v'],'xyz_nm':grid.xyz}
        atomic_npz(out/f'{name}.npz',**data);datasets[name]=data
        row={'name':name,'config':c.to_dict(),'electrons':float(n.sum()*grid.cell_volume*1e-21),
            'orthogonality_error':ortherr,'junction_missing_trace_norm':missing,'dimensions':ranks,
            'electrons_above_mu_plus_8kbt':float(tail),'runtime_s':time.time()-tic}
        rows.append(row);print(json.dumps(row),flush=True)
    reference=datasets['nested_6_12']
    for row in rows:
        q=datasets[row['name']]
        row['density_relative_to_largest_frozen_subspace']=float(np.linalg.norm(q['density_cm3']-reference['density_cm3'])/np.linalg.norm(reference['density_cm3']))
        row['transmission_l2_relative_to_largest_frozen_subspace']=float(np.linalg.norm(q['weighted_transmission']-reference['weighted_transmission'])/max(np.linalg.norm(reference['weighted_transmission']),1e-300))
    result={'source_sha256':source_digest(),'source_state':str(source),'mu_prescription':'same fixed μ for all frozen subspaces',
        'scf_recomputed':False,'plateau_certified':False,'cases':rows,
        'tail_interpretation':'changed subspace changes eigenstates and DOS across the entire Fermi window; added charge is not assumed to be an unoccupied high-energy tail'}
    write_json(out/'summary.json',result);return result


if __name__=='__main__':run(sys.argv[1],sys.argv[2])
