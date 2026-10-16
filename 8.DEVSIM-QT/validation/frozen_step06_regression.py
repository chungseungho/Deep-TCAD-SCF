"""Repair regression at the previously difficult frozen STEP06 GAA field."""
import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config import Config,VALLEYS
from modespace.hamiltonian import Grid,ModeBuilder
from solver import electronic,write_json,atomic_npz
from provenance import source_digest


def run(source,out):
    source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    s=np.load(source);c=Config(**json.loads(str(s['config_json'])))
    c.backend='negf';c.junction_basis='plain';c.energy_max_nodes=180000
    grid=Grid.build(c);builders=[ModeBuilder(c,grid,v) for v in VALLEYS]
    all_results=[]
    for step in (.02,.01):
        c.energy_step=step;tic=time.time()
        q=electronic(c,builders,grid,s['phi_v'],float(s['mu_l_ev']),float(s['vd_v']),True)
        all_results.append(q)
        atomic_npz(out/f'panel_{step:g}.npz',density_cm3=q['n'],bonds_a=q['bonds'],phi_v=s['phi_v'],
            xyz_nm=grid.xyz,mu_l_ev=q['mu'],mu_r_ev=q['mu']-float(s['vd_v']),
            source_sha256=source_digest(),input_source_sha256=str(s['source_sha256']),config_json=json.dumps(c.to_dict()))
        write_json(out/f'panel_{step:g}.json',{'current_a':q['id_a'],'runtime_s':time.time()-tic,
            'valleys':q['valleys'],'source_sha256':source_digest(),'config':c.to_dict()})
        print('panel',step,'current',q['id_a'],'runtime',time.time()-tic,flush=True)
    a,b=all_results
    report={'comparison':'semi-infinite NEGF integration on frozen STEP06 finite-spectral field; not an MS or NEGF SCF',
        'source_field':str(source),'source_sha256':source_digest(),
        'density_relative_volume_l2':float(np.linalg.norm(a['n']-b['n'])/np.linalg.norm(b['n'])),
        'current_relative_difference':float(abs(a['id_a']-b['id_a'])/abs(b['id_a'])),
        'currents_a':[a['id_a'],b['id_a']],'eta_added_ev':0.,
        'pass':bool(np.linalg.norm(a['n']-b['n'])/np.linalg.norm(b['n'])<5e-5 and abs(a['id_a']/b['id_a']-1)<5e-5)}
    write_json(out/'summary.json',report);print(json.dumps(report,indent=2));return report


if __name__=='__main__':run(sys.argv[1],sys.argv[2])
