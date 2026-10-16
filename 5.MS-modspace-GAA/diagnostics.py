#!/usr/bin/env python3
"""Actual one-bias convergence campaign, not manufactured device output.

Every case calls run.py. Failure is retained; comparisons use only converged
states. Pairwise agreement at one bias is not a whole-device certification.
"""
import argparse,csv,json,subprocess,sys
from pathlib import Path
import numpy as np
from config import Config


def cases():
    base=Config(dx=1,modes=3,sd_modes=6,mesh_h=3.5,channel_h=1.8,junction_h=1.3).to_dict()
    changes={
      'baseline':{},'mesh_medium':dict(mesh_h=2.5,channel_h=1.2,junction_h=.9),
      'mesh_fine':dict(mesh_h=1.8,channel_h=.8,junction_h=.6),
      'transverse_fine':dict(dy=.5,dz=.5),'longitudinal_fine':dict(dx=.5),
      'modes_low':dict(modes=2,sd_modes=4),'modes_high':dict(modes=5,sd_modes=12),
      'source_drain_24nm':dict(ls=24,ld=24),'source_drain_36nm':dict(ls=36,ld=36),
      'gamma_5mev':dict(gamma=.005),'gamma_20mev':dict(gamma=.02),
      'basis_local':dict(basis='local'),'mixing_unscreened':dict(preconditioner='none'),
    }
    return {k:{**base,**v} for k,v in changes.items()}


def summarize(out):
    rows=[]
    for directory in sorted(Path(out).glob('*')):
        f=directory/'point.csv'
        if not f.exists(): continue
        r=list(csv.DictReader(f.open()))[-1]
        row={k:r[k] for k in ['vg_v','vd_v','id_a','status','iterations','phi_residual_v','charge_residual',
                              'sd_neutrality_relative_error','current_conservation_relative','runtime_s']}
        row['case']=directory.name
        files=list((directory/'states').glob('vg*.npz'));files=[f for f in files if not f.name.endswith(('_spectra.npz','.tmp.npz'))]
        if files:
            s=np.load(files[-1]);row['electron_number']=float(s['electron_density_cm3'].sum()*float(s['cell_volume_nm3'])*1e-21)
            c=json.loads((directory/'settings.json').read_text())['config']
            mask=(s['x_nm']>=c['ls']+c['spacer'])&(s['x_nm']<=c['ls']+c['spacer']+c['lg'])
            row['barrier_ev']=float(np.min([s[f'valley{v}_subbands_ev'][mask,0].max() for v in range(3)]))
            row['channel_electrons']=float(s['electron_density_cm3'][(s['xyz_nm'][:,0]>=c['ls'])&(s['xyz_nm'][:,0]<=c['ls']+2*c['spacer']+c['lg'])].sum()*float(s['cell_volume_nm3'])*1e-21)
        rows.append(row)
    base=next((r for r in rows if r['case']=='baseline' and r['status']=='SCF_CONVERGED_PROTOTYPE'),None)
    if base:
        for r in rows:
            r['relative_current_vs_baseline']=float(r['id_a'])/float(base['id_a'])-1
            r['barrier_change_mev']=(r.get('barrier_ev',np.nan)-base['barrier_ev'])*1000
    if rows:
        keys=list(dict.fromkeys(k for r in rows for k in r))
        with (Path(out)/'convergence.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    (Path(out)/'summary.json').write_text(json.dumps({'cases':rows,'quantitative_validation_pass':False,
        'note':'Each variation needs successive refinement; one-bias sensitivity is not global convergence.'},indent=2))
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',default='output/convergence')
    p.add_argument('--cases',default='baseline,mesh_medium,modes_low,modes_high,gamma_5mev,gamma_20mev')
    p.add_argument('--vg',type=float,default=.8);p.add_argument('--vd',type=float,default=.1)
    p.add_argument('--summary-only',action='store_true');p.add_argument('--resume',action='store_true');a=p.parse_args()
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    if not a.summary_only:
        available=cases()
        for name in a.cases.split(','):
            if name not in available: raise ValueError('Unknown convergence case '+name)
            conf=out/(name+'.json');conf.write_text(json.dumps(available[name],indent=2))
            command=[sys.executable,'-u',str(Path(__file__).with_name('run.py')),'--phase','point','--vg',str(a.vg),
                     '--vd',str(a.vd),'--out',str(out/name),'--config',str(conf)]
            if a.resume: command.append('--resume')
            with (out/(name+'.log')).open('a' if a.resume else 'w') as f:
                ret=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT)
            print(name,'exit',ret.returncode,flush=True)
    print(json.dumps(summarize(out),indent=2))
