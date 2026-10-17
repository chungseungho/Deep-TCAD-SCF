from pathlib import Path
import numpy as np,json,hashlib,shutil
root=Path(__file__).resolve().parents[2];out=Path(__file__).parent
old=root/'results/final_uniform_vd005';ref=root/'results/final_uniform_equilibrium'
a=(old/'mesh/gaa_cm.msh').read_bytes();b=(ref/'mesh/gaa_cm.msh').read_bytes()
assert hashlib.sha256(a).digest()==hashlib.sha256(b).digest()
checkpoint=np.load(old/'checkpoints/vg0.400000_vd0.050000.npz');state=np.load(ref/'states/vg0.400000_vd0.000000.npz')
assert len(checkpoint['n'])==len(state['node_xyz_nm'])
np.savez_compressed(out/'restart_seed.npz',node_xyz_nm=state['node_xyz_nm'],poisson_density_cm3=checkpoint['n'])
shutil.copytree(old/'mesh',out/'mesh',dirs_exist_ok=True)
(out/'restart_conversion.json').write_text(json.dumps({'legacy_checkpoint':'../final_uniform_vd005/checkpoints/vg0.400000_vd0.050000.npz','iteration':int(checkpoint['iteration']),'identical_reference_mesh_sha256':hashlib.sha256(a).hexdigest(),'SCF_CONVERGED':False,'operation':'Coordinate metadata copied from an identical mesh; actual saved density reused without alteration.'},indent=2)+'\n')
