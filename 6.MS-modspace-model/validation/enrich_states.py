"""Add static material/mass metadata; never change a computed physical field.

Run after solver processes finish. Dataset eligibility remains false. This is
structured state packaging, not AI training or interpolation of failed data.
"""
import hashlib,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config import VALLEYS
from poisson.discretization import mesh_regions
from solver import atomic_npz


EXCLUDE={'development_point','negf_pilot','spectral_pilot'}


def run(root):
    root=Path(root);catalog=[];enrichment_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for settings in sorted(root.rglob('settings.json')):
        directory=settings.parent
        if any(p in EXCLUDE for p in directory.relative_to(root).parts):continue
        statefiles=list((directory/'states').glob('vg*.npz'))
        statefiles=[p for p in statefiles if not p.name.endswith('_spectra.npz')]
        if not statefiles:continue
        record=json.loads(settings.read_text());mesh=directory/'mesh/gaa_nm.msh'
        if mesh.exists():
            xyz,tets=mesh_regions(mesh)
            static={'mesh_vertices_nm':xyz,'material_names':np.array(list(tets)),
                    'geometry_config_json':json.dumps(record['config'],sort_keys=True)}
            for material,indices in tets.items():static['tetra_'+material]=indices
            atomic_npz(directory/'mesh_static.npz',**static)
        for file in sorted(statefiles):
            s=np.load(file);data={name:s[name] for name in s.files}
            data.update(valley_names=np.array([v.name for v in VALLEYS]),
                        effective_masses_m0=np.array([[v.mx,v.my,v.mz] for v in VALLEYS]),
                        valley_pair_multiplicity=np.array([v.multiplicity for v in VALLEYS]),
                        spin_degeneracy=2,density_already_degeneracy_weighted=True,
                        material_distribution_reference='../mesh_static.npz',ai_teacher_eligible=False,
                        state_packaging_source_sha256=enrichment_sha)
            atomic_npz(file,**data)
            catalog.append({'state':str(file.relative_to(root)),
                'vg_v':float(data['vg_v']),'vd_v':float(data['vd_v']),
                'scf_converged':bool(data['converged']),'ai_teacher_eligible':False,
                'quantitative_validation_pass':bool(data['quantitative_validation_pass']),
                'source_sha256':str(data['source_sha256']),
                'configuration_fingerprint':str(data['configuration_fingerprint']),
                'extension_source_sha256':str(data.get('extension_source_sha256','')),
                'material_distribution':str((directory/'mesh_static.npz').relative_to(root))})
    result={'states':catalog,'teacher_eligible_states':0,
            'catalog_count_scope':'File entries include copied sweep partitions and shared convergence baselines, not independent solves or unique training samples.',
            'reason':'contact, joint discretization/basis and independent biased-lead validation incomplete; no learning was performed'}
    (root/'dataset_catalog.json').write_text(json.dumps(result,indent=2));return result


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--root',default='results');a=p.parse_args()
    result=run(a.root);print(json.dumps({'packaged_states':len(result['states']),'teacher_eligible_states':0}))
