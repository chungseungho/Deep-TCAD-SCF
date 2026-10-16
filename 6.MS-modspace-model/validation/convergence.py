#!/usr/bin/env python3
"""Real SCF campaign. Plan first; --run executes local Python, never SLURM.

Three levels per axis are required. Frozen-potential comparisons are separate
and cannot certify coupled SCF convergence. Failed/missing levels fail plateau.
"""
import argparse,csv,json,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config import Config
from resources import estimate


def campaign(backend='spectral'):
    base=Config(backend=backend,dx=.5,modes=3,sd_modes=6,mesh_h=3.5,channel_h=1.8,junction_h=1.3).to_dict()
    groups={
        'mesh':[dict(mesh_h=3.5,channel_h=1.8,junction_h=1.3),dict(mesh_h=2.5,channel_h=1.2,junction_h=.9),dict(mesh_h=1.8,channel_h=.8,junction_h=.6)],
        'dx':[dict(dx=.5),dict(dx=.25),dict(dx=.125)],
        'transverse':[dict(dy=1.,dz=1.),dict(dy=.5,dz=.5),dict(dy=.25,dz=.25)],
        'modes':[dict(modes=3,sd_modes=6),dict(modes=6,sd_modes=12),dict(modes=9,sd_modes=18)],
        'length':[dict(ls=12.,ld=12.),dict(ls=24.,ld=24.),dict(ls=36.,ld=36.)],
        'gamma':[dict(gamma=.005),dict(gamma=.01),dict(gamma=.02),dict(gamma=.04)],
        'scf':[dict(phi_tol=1e-4,charge_tol=2e-3,current_tol=2e-3),dict(phi_tol=5e-5,charge_tol=1e-3,current_tol=1e-3),dict(phi_tol=2.5e-5,charge_tol=5e-4,current_tol=5e-4)],
        'energy':[dict(energy_step=.02),dict(energy_step=.01),dict(energy_step=.005)],
        'poles':[dict(fermi_poles=32),dict(fermi_poles=48),dict(fermi_poles=64)],
        'basis':[dict(basis='fixed',modes=3,sd_modes=6),dict(basis='local',modes=3,sd_modes=6),dict(basis='fixed',modes=9,sd_modes=18),dict(basis='local',modes=9,sd_modes=18)]}
    if backend=='negf': groups.pop('gamma')
    else:
        groups.pop('energy');groups.pop('poles')
    return {name:[{**base,**change} for change in changes] for name,changes in groups.items()}


def summarize(out):
    out=Path(out);rows=[];plateaus={};planned={};incomplete=[]
    for conf in sorted(out.glob('*.json')):
        if conf.name in ('summary.json','plan.json'): continue
        configuration=json.loads(conf.read_text())
        if not isinstance(configuration,dict) or 'backend' not in configuration or 'dx' not in configuration: continue
        name=conf.stem;table=out/name/'point.csv'
        planned.setdefault(name.rsplit('_',1)[0],[]).append(name)
        if not table.exists():
            history=out/name/'scf_history.csv'
            attempts=list(csv.DictReader(history.open())) if history.exists() else []
            incomplete.append({'case':name,
                'execution_status':'NO_ACCEPTED_FINAL_STATE' if attempts else 'NOT_EXECUTED_OR_NO_SCF_ITERATE',
                'last_saved_iteration':attempts[-1] if attempts else None,
                'checkpoint_exists':bool(list((out/name/'checkpoints').glob('*.npz')))})
            continue
        for row in csv.DictReader(table.open()):
            row={**row,'case':name};state=out/name/'states'/(row['key']+'.npz')
            if state.exists():
                s=np.load(state);c=json.loads(conf.read_text())
                row['electron_number']=float(s['electron_density_cm3'].sum()*s['cell_volume_nm3']*1e-21)
                channel=(s['xyz_nm'][:,0]>=c['ls'])&(s['xyz_nm'][:,0]<=c['ls']+2*c['spacer']+c['lg'])
                row['channel_electrons']=float(s['electron_density_cm3'][channel].sum()*s['cell_volume_nm3']*1e-21)
                gate=(s['x_nm']>=c['ls']+c['spacer'])&(s['x_nm']<=c['ls']+c['spacer']+c['lg'])
                row['barrier_ev']=float(min(s[f'valley{j}_subbands_ev'][gate,0].max() for j in range(3)))
            rows.append(row)
    for group in sorted(planned):
        r=sorted([v for v in rows if v['case'].rsplit('_',1)[0]==group],key=lambda v:v['case'])
        missing=sorted(set(planned[group])-{v['case'] for v in r})
        good=len(r)>=3 and not missing and all(v['status']=='SCF_CONVERGED_PROTOTYPE' for v in r)
        changes=[]
        if good:
            for a,b in zip(r,r[1:]):
                changes.append({'current':abs(float(a['id_a'])-float(b['id_a']))/max(abs(float(b['id_a'])),1e-18),
                                'channel_charge':abs(a['channel_electrons']/b['channel_electrons']-1) if b['channel_electrons'] else float('inf'),
                                'barrier_v':abs(a['barrier_ev']-b['barrier_ev'])})
        plateau=good and all(v['current']<.02 and v['channel_charge']<.02 and v['barrier_v']<.001 for v in changes)
        plateaus[group]={'planned_levels':len(planned[group]),'completed_levels':len(r),'missing_cases':missing,
                         'three_or_more_converged_levels':good,'successive_changes':changes,'plateau_pass':bool(plateau)}
    if rows:
        keys=list(dict.fromkeys(k for r in rows for k in r))
        with (out/'convergence.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    result={'cases':rows,'incomplete_cases':incomplete,'plateaus':plateaus,'global_reference_certified':False,
            'criteria':'All successive current/charge changes <2%, barrier <1mV, >=3 levels; single bias never certifies full domain'}
    (out/'summary.json').write_text(json.dumps(result,indent=2));return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--backend',choices=['spectral','negf'],default='spectral')
    p.add_argument('--groups',default='mesh,dx,transverse,modes,length,gamma,basis,scf')
    p.add_argument('--vg',type=float,default=.8);p.add_argument('--vd',type=float,default=.1)
    p.add_argument('--out',default='results/convergence');p.add_argument('--run',action='store_true')
    p.add_argument('--resume',action='store_true');p.add_argument('--summary-only',action='store_true')
    a=p.parse_args();out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    if a.summary_only: print(json.dumps(summarize(out),indent=2));return
    available=campaign(a.backend);plan=[]
    for group in a.groups.split(','):
        if group not in available: raise ValueError('Not applicable to selected backend: '+group)
        for i,conf in enumerate(available[group]):
            name=f'{group}_{i}';f=out/(name+'.json');f.write_text(json.dumps(conf,indent=2))
            command=[sys.executable,'-u',str(Path(__file__).resolve().parents[1]/'run.py'),'--phase','point','--vg',str(a.vg),'--vd',str(a.vd),'--config',str(f),'--out',str(out/name)]
            if a.resume: command.append('--resume')
            plan.append({'case':name,'config':conf,'resources':estimate(Config(**conf)),'command':command})
    (out/'plan.json').write_text(json.dumps(plan,indent=2))
    print(json.dumps([{'case':v['case'],'resources':v['resources']} for v in plan],indent=2),flush=True)
    if a.run:
        for item in plan:
            log=out/(item['case']+'.log')
            with log.open('a' if a.resume else 'w') as f: result=subprocess.run(item['command'],stdout=f,stderr=subprocess.STDOUT)
            print(item['case'],'exit',result.returncode,flush=True)
        summarize(out)


if __name__=='__main__': main()
