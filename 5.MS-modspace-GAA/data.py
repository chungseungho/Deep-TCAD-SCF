#!/usr/bin/env python3
"""Export whole-bias snapshots for future DeepSCF training; preserve validation flags.

Input is the donor-neutral, non-SCF Poisson field and boundary biases. Target
is converged CHANNEL charge. Converged phi is a diagnostic, never a one-shot
input. Full bias points, not randomly shuffled spatial cells, form the split.
"""
import argparse,json
from pathlib import Path
import numpy as np


def export_dataset(run_dir,destination,include_prototypes=False):
    out=Path(run_dir);settings=json.loads((out/'settings.json').read_text());c=settings['config']
    files=[];seen=set()
    for f in sorted((out/'states').glob('vg*.npz')):
        if f.name.endswith(('_spectra.npz','.tmp.npz')) or f.name in seen: continue
        s=np.load(f)
        if not bool(s['converged']): continue
        if not bool(s['quantitative_validation_pass']) and not include_prototypes: continue
        if 'initial_neutral_poisson_phi_v' not in s: continue
        seen.add(f.name);files.append(f)
    if not files: raise ValueError('No eligible snapshots; prototype data requires explicit --include-prototypes')
    first=np.load(files[0]);xyz=first['xyz_nm']
    # Channel includes the gate-controlled bridge and the two nanosheet extensions.
    region=(xyz[:,0]>=c['ls'])&(xyz[:,0]<=c['ls']+2*c['spacer']+c['lg'])
    features=[];target=[];bias=[];valid=[];phi=[]
    for f in files:
        s=np.load(f)
        if not np.array_equal(s['xyz_nm'],xyz): raise ValueError('Grid mismatch; separate datasets by geometry')
        features.append(s['initial_neutral_poisson_phi_v']);target.append(s['electron_density_cm3'][region])
        phi.append(s['phi_v'][region]);bias.append([float(s['vg_v']),float(s['vd_v'])]);valid.append(bool(s['quantitative_validation_pass']))
    # Deterministic whole-bias test set; never split cells from the same solution.
    test=np.arange(len(files))%5==4
    destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(destination,input_poisson_phi_v=np.array(features),bias_v=np.array(bias),
                        target_channel_density_cm3=np.array(target),target_channel_phi_v=np.array(phi),
                        channel_xyz_nm=xyz[region],full_xyz_nm=xyz,channel_mask=region,
                        test_mask=test,quantitative_validation_pass=np.array(valid),
                        config_json=json.dumps(c),fingerprint=settings['fingerprint'],
                        source_files=np.array([f.name for f in files]))
    return {'snapshots':len(files),'train_points':int((~test).sum()),'heldout_bias_points':int(test.sum()),
            'quantitatively_validated_points':int(sum(valid)),'dataset':str(destination)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run_dir');p.add_argument('destination')
    p.add_argument('--include-prototypes',action='store_true');a=p.parse_args()
    print(json.dumps(export_dataset(a.run_dir,a.destination,a.include_prototypes),indent=2))
