#!/usr/bin/env python3
"""Read-only audit reruns the preserved STEP05 Poisson at a saved SCF density.

This optional audit needs STEP05 next to STEP06. The STEP06 simulator itself
never imports or needs any previous stage.
"""
import argparse,json,subprocess,sys
from pathlib import Path


def run(old,out):
    old=Path(old).resolve()
    program='''
import numpy as np,json
from pathlib import Path
from scipy.spatial import cKDTree
from config import Config
from poisson.device import Poisson
c=Config(**json.loads(Path('results/transfer_settings.json').read_text())['config'])
s=np.load('results/states/vg0.400000_vd0.100000.npz')
p=Poisson(c,Path('results/mesh_iv/gaa_cm.msh'))
try:
 p.set_charge(s['poisson_density_cm3']);p.solve(.4,.1);report={}
 for a,b in [('silicon','oxide'),('silicon','spacer'),('oxide','spacer')]:
  distance,index=cKDTree(p.xyz[a]).query(p.xyz[b]);keep=distance<1e-6
  error=abs(p.values('Potential',a)[index[keep]]-p.values('Potential',b)[keep])
  report[a+'_'+b]={'shared_nodes':int(keep.sum()),'maximum_potential_jump_v':float(error.max()),'rms_jump_v':float(np.sqrt(np.mean(error**2)))}
 print('AUDIT_JSON='+json.dumps(report))
finally:p.close()
'''
    process=subprocess.run([sys.executable,'-c',program],cwd=old,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    out.with_suffix('.log').write_text(process.stdout)
    if process.returncode: raise RuntimeError('STEP05 audit failed; see '+str(out.with_suffix('.log')))
    result=json.loads(next(v[len('AUDIT_JSON='):] for v in process.stdout.splitlines() if v.startswith('AUDIT_JSON=')))
    result.update(saved_bias={'vg_v':.4,'vd_v':.1},source='preserved STEP05 original Poisson and saved SCF density',
                  fix='STEP06 uses one shared potential DOF for all touching materials')
    out.write_text(json.dumps(result,indent=2));return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--old',default='../5.MS-modspace-GAA');p.add_argument('--out',default='results/step05_interface_audit.json');a=p.parse_args()
    print(json.dumps(run(a.old,a.out),indent=2))
