"""Compare separately converged coupled solutions, with ensemble checks."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config import Config
from device.regions import five_regions,NAMES
from solver import write_json,atomic_npz
from provenance import source_digest


def compare(a_path,b_path,out):
    a=np.load(a_path);b=np.load(b_path);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if not bool(a['converged']) or not bool(b['converged']):raise ValueError('Need two numerically converged states')
    if not np.array_equal(a['xyz_nm'],b['xyz_nm']):raise ValueError('Different grids')
    if not np.array_equal(a['node_xyz_nm'],b['node_xyz_nm']):raise ValueError('Different Poisson nodes')
    ca=Config(**json.loads(str(a['config_json'])));cb=Config(**json.loads(str(b['config_json'])))
    n1=a['electron_density_cm3'];n2=b['electron_density_cm3'];p1=a['phi_v'];p2=b['phi_v'];xyz=a['xyz_nm']
    channel=(xyz[:,0]>=ca.gate_start)&(xyz[:,0]<=ca.gate_end)
    regions=five_regions(xyz[:,0],ca);dv=float(a['cell_volume_nm3'])*1e-21
    matches=abs(float(a['mu_l_ev'])-float(b['mu_l_ev']))<1e-12 and abs(float(a['mu_r_ev'])-float(b['mu_r_ev']))<1e-12
    regional={}
    for i,name in enumerate(NAMES):
        use=regions==i;regional[name]={'ms_electrons':float(n1[use].sum()*dv),'negf_electrons':float(n2[use].sum()*dv)}
    longitudinal=(a['x_nm']>=ca.ls)&(a['x_nm']<=ca.drain_start)
    barriers=[]
    for j in range(3):
        ba=float(max(a[f'valley{j}_subbands_ev'][longitudinal,0]))
        bb=float(max(b[f'valley{j}_subbands_ev'][longitudinal,0]))
        barriers.append({'valley':j,'ms_barrier_ev':ba,'negf_barrier_ev':bb,'difference_ev':ba-bb})
    result={'source_sha256':source_digest(),'a':str(a_path),'b':str(b_path),'coupled_fixed_points':True,
        'matched_chemical_potentials':bool(matches),'a_mu_l_ev':float(a['mu_l_ev']),'b_mu_l_ev':float(b['mu_l_ev']),
        'density_relative_volume_l2':float(np.linalg.norm(n1-n2)/np.linalg.norm(n2)),
        'maximum_quantum_phi_difference_v':float(max(abs(p1-p2))),
        'maximum_node_phi_difference_v':float(max(abs(a['node_phi_v']-b['node_phi_v']))),
        'channel_mean_barrier_difference_ev':float(np.mean(-p1[channel]+p2[channel])),
        'a_current_a':float(a['id_a']),'b_current_a':float(b['id_a']),'regions':regional,'barriers':barriers,
        'a_config':ca.to_dict(),'b_config':cb.to_dict(),'quantitative_validation_pass':False,
        'interpretation':'same-contact comparison' if matches else 'different finite-N/open-reservoir ensembles; not same-contact equivalence'}
    write_json(out/'summary.json',result)
    atomic_npz(out/'comparison.npz',xyz_nm=xyz,n_ms_cm3=n1,n_negf_cm3=n2,phi_ms_v=p1,phi_negf_v=p2,
        source_sha256=source_digest(),a_source_sha256=a['source_sha256'],b_source_sha256=b['source_sha256'])
    print(json.dumps(result,indent=2));return result


if __name__=='__main__':compare(sys.argv[1],sys.argv[2],sys.argv[3])
