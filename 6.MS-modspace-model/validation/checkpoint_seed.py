#!/usr/bin/env python3
"""Convert an incomplete checkpoint into an explicitly unaccepted initial guess."""
import hashlib
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config import Config
from poisson.device import Poisson


def seed(directory,destination):
    directory=Path(directory);settings=json.loads((directory/'settings.json').read_text())
    files=sorted((directory/'checkpoints').glob('*.npz'))
    if len(files)!=1:raise ValueError('Use a single-bias incomplete run with exactly one checkpoint')
    checkpoint=files[0];data=np.load(checkpoint);c=Config(**settings['config'])
    p=Poisson(c,directory/'mesh/gaa_cm.msh')
    try:
        if data['n'].shape!=p.donors.shape:raise ValueError('Checkpoint and Poisson node dimensions differ')
        destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(destination,node_xyz_nm=p.xyz['silicon'],poisson_density_cm3=data['n'],
            vg_v=data['vg'],vd_v=data['vd'],seed_only=True,converged=False,
            quantitative_validation_pass=False,ai_teacher_eligible=False,
            checkpoint_iteration=data['iteration'],configuration_fingerprint=settings['fingerprint'],
            checkpoint_source_sha256=settings['source_sha256'],
            checkpoint_file_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest())
        return {'seed':str(destination),'seed_only':True,'accepted_result':False,
                'iteration':int(data['iteration']),'note':'Re-run the full quantum-Poisson SCF on the new output; no current/density acceptance is implied by this seed.'}
    finally:p.close()


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('directory');parser.add_argument('--out',required=True);args=parser.parse_args()
    print(json.dumps(seed(args.directory,args.out),indent=2))
