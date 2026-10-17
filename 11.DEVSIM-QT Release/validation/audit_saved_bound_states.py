import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config import Config
from device.definition import Frame
from modespace.hamiltonian import Grid,ModeBuilder
from modespace.junction_basis import JunctionCompatibleBuilder
from validation.bound_states import audit
from solver import write_json
from provenance import source_digest


def run(path,out):
    s=np.load(path);c=Config(**json.loads(str(s['config_json'])));grid=Grid.build(c)
    cls=JunctionCompatibleBuilder if c.junction_basis=='nested' else ModeBuilder
    results=[]
    for v in Frame(c.transport_axis).valleys():
        system=cls(c,grid,v).build(s['phi_v']);result=audit(system,c,float(s['mu_l_ev']),float(s['vd_v']))
        results.append(result);print(v.name,result,flush=True)
    report={'source_sha256':source_digest(),'state_source_sha256':str(s['source_sha256']),
        'state':str(path),'config':c.to_dict(),'valleys':results,
        'needs_preparation_model':any(x['bias_dependent_bound_occupation_unresolved'] for x in results),
        'scope':'below all lead bands only; embedded dark states not excluded'}
    write_json(out,report);return report


if __name__=='__main__':run(sys.argv[1],sys.argv[2])
