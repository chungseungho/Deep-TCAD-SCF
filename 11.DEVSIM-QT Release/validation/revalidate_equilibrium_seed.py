"""Explicit equilibrium seed migration after a reviewed numerical/metadata source change.

Refuses biased states and unlisted exact old/new production-file hash pairs.
Recomputes the entire equilibrium quantum density using current source before
creating a NON-certified starting density. Original state remains untouched.
"""
import sys,json,hashlib
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config import Config
from modespace.hamiltonian import Grid,ModeBuilder
from modespace.junction_basis import JunctionCompatibleBuilder
from device.definition import DeviceDefinition
from transport.registry import create_transport_backend
from solver import electronic
from provenance import source_digest

def migrate(state,manifest,out):
    root=Path(__file__).resolve().parents[1];s=np.load(state,allow_pickle=False);old=json.loads(Path(manifest).read_text());changes=[]
    if float(s['vd_v'])!=0:raise ValueError('Only zero-bias equilibrium seeds can be migrated')
    for rel,sha in old.items():
        p=Path(rel)
        if not rel.endswith('.py') or any(x in ('validation','tests','results','__pycache__') for x in p.parts):continue
        if hashlib.sha256((root/rel).read_bytes()).hexdigest()!=sha:changes.append(rel)
    allow=json.loads((root/'validation/equilibrium_migration_pairs.json').read_text())
    for rel in changes:
        current=hashlib.sha256((root/rel).read_bytes()).hexdigest()
        if [old[rel],current] not in allow.get(rel,[]):raise ValueError('Unreviewed source version pair: '+rel)
    c=Config(**json.loads(str(s['config_json'])));g=Grid.build(c);cls=JunctionCompatibleBuilder if c.junction_basis=='nested' else ModeBuilder
    builders=[cls(c,g,v) for v in DeviceDefinition(transport_axis=c.transport_axis).frame.valleys()];backend=create_transport_backend(c.backend,c);mu,cal=backend.calibrate_contacts({'builders':builders,'grid':g,'config':c})
    q=electronic(backend,builders,g,s['phi_v'],mu,0);error=float(np.linalg.norm(q['n']-s['electron_density_cm3'])/np.linalg.norm(q['n']))
    if error>1e-10:
        Path(out).with_suffix('.rejected.json').write_text(json.dumps({'status':'FAIL','density_relative_error':error,'maximum_density_difference_cm3':float(max(abs(q['n']-s['electron_density_cm3']))),'input':str(state)},indent=2))
        raise ValueError(f'Equilibrium density not reproduced: relative error {error:.12g}')
    report={'input_state':str(state),'input_source_sha256':str(s['source_sha256']),'current_source_sha256':source_digest(),'changed_files':changes,'equilibrium_density_relative_error':error,'meaning':'revalidated initial guess only; next SCF must converge and reevaluate','TEACHER_READY':False}
    np.savez_compressed(out,config_json=s['config_json'],source_sha256=source_digest(),node_xyz_nm=s['node_xyz_nm'],poisson_density_cm3=s['poisson_density_cm3'],numerically_converged=False,origin_json=json.dumps(report))
    Path(out).with_suffix('.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':migrate(*sys.argv[1:])
