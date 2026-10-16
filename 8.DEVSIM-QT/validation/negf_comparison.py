"""Frozen-field density comparison: never relabel it a coupled NEGF SCF."""
import json
from pathlib import Path
from dataclasses import replace
import numpy as np
from config import Config
from modespace.hamiltonian import Grid,ModeBuilder
from modespace.junction_basis import JunctionCompatibleBuilder
from device.definition import Frame
from device.regions import five_regions,NAMES
from solver import electronic,write_json,atomic_npz
from provenance import source_digest


def frozen(path,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);s=np.load(path)
    c=Config(**json.loads(str(s['config_json'])));g=Grid.build(c)
    cls=JunctionCompatibleBuilder if c.junction_basis=='nested' else ModeBuilder
    builders=[cls(c,g,v) for v in Frame(c.transport_axis).valleys()]
    mu=float(s['mu_l_ev']);vd=float(s['vd_v'])
    b=electronic(replace(c,backend='negf'),builders,g,s['phi_v'],mu,vd,True)
    a=s['electron_density_cm3'];d=b['n'];dv=g.cell_volume*1e-21
    regions=five_regions(g.xyz[:,0],c);regional={}
    for i,name in enumerate(NAMES):
        select=regions==i
        regional[name]={'ms_electrons':float(a[select].sum()*dv),
            'negf_electrons':float(d[select].sum()*dv),
            'density_relative_l2':float(np.linalg.norm(a[select]-d[select])/max(np.linalg.norm(d[select]),1e-300)),
            'charge_difference_c':float((a[select]-d[select]).sum()*dv*1.602176634e-19)}
    from transport.ms_postprocess import transmission,direct_transmission
    energies=np.linspace(min(mu,mu-vd)-.2,max(mu,mu-vd)+.3,101)
    T=sum(2*sys.valley.multiplicity*transmission(sys,c,energies,vd) for sys in b['systems'])
    dense_indices=[20,40,60,80];dense_defect=0.
    for sys in b['systems']:
        r=transmission(sys,c,energies[dense_indices],vd)
        direct=direct_transmission(sys,c,energies[dense_indices],vd)
        dense_defect=max(dense_defect,float(max(abs(r-direct))))
    result={'comparison':'n_MS versus n_NEGF[phi_MS], frozen field, same μ and lead bands',
        'coupled_negf_fixed_point':False,'ms_state':str(path),'source_sha256':source_digest(),
        'input_source_sha256':str(s['source_sha256']),'config':c.to_dict(),
        'density_relative_volume_l2':float(np.linalg.norm(a-d)/max(np.linalg.norm(d),1e-300)),
        'ms_electrons':float(a.sum()*dv),'negf_electrons':float(d.sum()*dv),
        'ms_current_a':float(s['id_a']),'negf_current_a':b['id_a'],
        'regions':regional,'dense_vs_rgf_max_transmission_defect':dense_defect,
        'negf_diagnostics':b['valleys'],'quantitative_validation_pass':False}
    atomic_npz(out/'comparison.npz',xyz_nm=g.xyz,n_ms_cm3=a,n_negf_cm3=d,phi_ms_v=s['phi_v'],
               energies_ev=energies,weighted_transmission=T,source_sha256=source_digest(),
               config_json=json.dumps(c.to_dict()),mu_l_ev=mu,mu_r_ev=mu-vd)
    write_json(out/'comparison.json',result);return result
