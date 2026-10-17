"""Source/config identity excludes generated output and tests by definition."""
import hashlib,json
from pathlib import Path


def fingerprint(c):
    return hashlib.sha256(json.dumps(c.to_dict(),sort_keys=True).encode()).hexdigest()


def source_digest():
    root=Path(__file__).resolve().parent;h=hashlib.sha256()
    for f in sorted(root.rglob('*.py')):
        if any(p in ('output','results','tests','validation','__pycache__','.venv','venv','.git','.pytest_cache') for p in f.relative_to(root).parts): continue
        h.update(str(f.relative_to(root)).encode());h.update(f.read_bytes())
    return h.hexdigest()


def validate_restart(path,c,xyz):
    import numpy as np
    saved=np.load(path,allow_pickle=False)
    required={'config_json','source_sha256','node_xyz_nm','poisson_density_cm3'}
    if not required.issubset(saved.files):raise ValueError('Restart missing provenance or coordinates')
    if str(saved['source_sha256'])!=source_digest():raise ValueError('Restart source mismatch')
    old=json.loads(str(saved['config_json']))
    if old!=c.to_dict():raise ValueError('Restart model/geometry/config mismatch')
    if saved['node_xyz_nm'].shape!=xyz.shape or not np.allclose(saved['node_xyz_nm'],xyz,rtol=0,atol=1e-9):
        raise ValueError('Restart mesh coordinate mismatch')
    density=saved['poisson_density_cm3']
    if density.shape!=(len(xyz),) or not np.isfinite(density).all() or np.min(density)<0:raise ValueError('Invalid restart density')
    return saved
