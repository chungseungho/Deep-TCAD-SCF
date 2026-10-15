#!/usr/bin/env python3
"""STEP03 production-oriented Id-Vgs -> Vth -> three Id-Vds curves.

Uses STEP01 DEVSIM Poisson and geometry unchanged, STEP01 valley EMA inputs
unchanged; STEP03 replaces initial demonstrator transport with a vectorized,
fixed-basis coupled-mode-space ballistic NEGF. Never writes fabricated data.

Must run from 3.modsspace-iv-test/ inside a DEVSIM+gmsh venv.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import replace,asdict
import json
import os
import shutil
from pathlib import Path
import sys
import time
import traceback
import numpy as np
from scipy.interpolate import RegularGridInterpolator

HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'1.example_build'
if not BASE.joinpath('config.py').exists():
    raise SystemExit('Missing ../1.example_build/config.py; restore STEP01 reference directory.')
sys.path.insert(0,str(BASE))
from config import Config
from poisson.devsim_poisson import DevsimPoisson
from geometry.make_mesh import build_mesh
from quantum import make_bases,solve_quantum,min_bandwidth_ev
from analysis import voltages,summarize_iv,choose_vgs
from plotting import make_modes_figure,make_idvgs,make_idvds,make_device_maps,make_scf_history
from quality import verify_snapshot


def axis_grid(cfg):
    return (np.linspace(0,cfg.length_nm,cfg.nx),
            np.linspace(-cfg.width_nm/2,cfg.width_nm/2,cfg.ny+2)[1:-1],
            np.linspace(-cfg.thickness_nm/2,cfg.thickness_nm/2,cfg.nz+2)[1:-1])


def charge_to_nodes(p,axes,n):
    # Clamp source/drain ends and semiconductor surfaces to the nearest valid
    # interior quantum-site value rather than silently zeroing surface charge.
    q=p.si_nodes_nm.copy()
    for k,a in enumerate(axes):q[:,k]=np.clip(q[:,k],a[0],a[-1])
    ans=RegularGridInterpolator(axes,n,method='linear',bounds_error=True)(q)
    if np.any(~np.isfinite(ans)) or np.any(ans<0):
        raise ArithmeticError('Charge interpolation produced invalid densities')
    return ans


def append_dict(path,row,columns):
    path=Path(path);exists=path.exists()
    with path.open('a',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns)
        if not exists:w.writeheader()
        w.writerow({k:row.get(k,'') for k in columns})


def already(path,key):
    if not Path(path).exists():return False
    with open(path,newline='') as f:
        return any(r['key']==key and r['status']=='CONVERGED' for r in csv.DictReader(f))


def state_path(out,key):return out/'states'/f'{key}.npz'


def run_bias(p,cfg,bases,vg,vd,args,axes,out,key):
    tic=time.perf_counter()
    # Retain the previous converged carrier density during continuous bias sweep.
    # On process restart, restore the nearest saved converged state if available.
    p.solve(vg-args.gate_flatband,vd)
    phi=p.sample_si_potential(*axes)
    old_current=None
    converged=False
    hist=[]
    for it in range(1,args.max_iterations+1):
        current,n,vc,_,_=solve_quantum(cfg,phi,bases,vd,args.energy_step,density=True)
        nn=charge_to_nodes(p,axes,n)
        old_nodes=np.asarray(p.d.get_node_model_values(device=p.device,region=p.si,name='InjectedElectrons'),float)
        # fixed-point residual BEFORE damping, not the artificially small mixed update.
        nres=float(np.linalg.norm(nn-old_nodes)/max(np.linalg.norm(nn),1e-30))
        new_nodes=(1-args.mix)*old_nodes+args.mix*nn
        p.set_charge(new_nodes)
        p.solve(vg-args.gate_flatband,vd)
        new_phi=p.sample_si_potential(*axes)
        dphi=float(np.max(np.abs(new_phi-phi)))
        icurrent=abs(current-old_current)/max(abs(current),1e-22) if old_current is not None else float('inf')
        row={'key':key,'iteration':it,'id_a':current,'delta_phi_v':dphi,
             'charge_residual':nres,'relative_current_change':icurrent,
             'elapsed_s':time.perf_counter()-tic}
        hist.append(row)
        print(f'[{key}] SCF {it:02d}/{args.max_iterations:02d} Id={current:.5e} A  '
              f'dPhi={dphi:.3e} V  chargeR={nres:.3e}  currentR={icurrent:.3e}',flush=True)
        phi=new_phi
        if (it>=3 and dphi < args.phi_tol and nres < args.charge_tol
                and icurrent < args.current_tol):
            converged=True
            break
        old_current=current
    # Terminal NEGF on final potential, not outdated previous Poisson field.
    final_I,final_n,vc,E,T=solve_quantum(cfg,phi,bases,vd,args.energy_step,density=True)
    nres_end=float(np.linalg.norm(charge_to_nodes(p,axes,final_n)-np.asarray(
         p.d.get_node_model_values(device=p.device,region=p.si,name='InjectedElectrons'),float))
         /max(np.linalg.norm(charge_to_nodes(p,axes,final_n)),1e-30))
    # Only completed SCF iterates with a consistent final charge are accepted.
    converged=bool(converged and nres_end<args.charge_tol*1.2)
    status='CONVERGED' if converged else 'NOT_CONVERGED'
    sp=state_path(out,key)
    np.savez_compressed(sp,x_nm=axes[0],y_nm=axes[1],z_nm=axes[2],
        phi_v=phi,electron_density_cm3=final_n,id_a=final_I,vg_v=vg,vd_v=vd,
        converged=int(converged),injected_si_cm3=np.asarray(
         p.d.get_node_model_values(device=p.device,region=p.si,name='InjectedElectrons')))
    for r in hist:append_dict(out/'scf_history.csv',r,
       ['key','iteration','id_a','delta_phi_v','charge_residual','relative_current_change','elapsed_s'])
    summary={'key':key,'vg_v':float(vg),'vd_v':float(vd),'id_a':float(final_I),
        'status':status,'iterations':len(hist),'delta_phi_v':float(hist[-1]['delta_phi_v']),
        'charge_residual':nres_end,'relative_current_change':float(hist[-1]['relative_current_change']),
        'runtime_s':time.perf_counter()-tic,'state_file':str(sp.relative_to(out)),
        'valley_currents_a':vc}
    print(f'[POINT {status}] {key}: Id={final_I:.6e} A, '
          f'SCF={len(hist)}, final charge residual={nres_end:.3e}, t={summary["runtime_s"]:.1f}s',flush=True)
    return summary


def upsert_point(path,row,columns):
    """Atomic per-key replacement prevents duplicate rows after failed/restarted bias."""
    path=Path(path)
    previous=get_rows(path) if path.exists() else []
    records=[r for r in previous if r['key']!=row['key']]+[row]
    records.sort(key=lambda r:(float(r['vg_v']),float(r['vd_v'])))
    tmp=path.with_suffix('.tmp')
    with tmp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore')
        w.writeheader()
        for r in records:w.writerow(r)
    os.replace(tmp,path)


def fill_state(p,out,step,key):
    """Allow restart without repeating earlier converged biases, via last state."""
    f=state_path(out,key)
    if f.exists():
        saved=np.load(f)
        nodes=saved['injected_si_cm3']
        if nodes.shape == (len(p.si_nodes_nm),):
            p.set_charge(nodes)
            return True
    return False


def get_rows(path):
    if not path.exists():return []
    with path.open(newline='') as f:return list(csv.DictReader(f))


def write_report(out,cfg,args,summary=None,triplet=None):
    rows=get_rows(out/'idvgs.csv')
    vds=get_rows(out/'idvds.csv')
    report=['# STEP03 — Real-device Mode-Space NEGF Validation', '',
       '> All plots are calculated from actual DEVSIM coupled SCF runs. This is an uncalibrated ballistic EMA prototype, NOT measurement-grade silicon TCAD.',
       '',f'- Id–Vgs converged: {sum(r["status"]=="CONVERGED" for r in rows)}/{len(rows)}',
       f'- Id–Vds converged: {sum(r["status"]=="CONVERGED" for r in vds)}/{len(vds)}',
       f'- Transverse: {cfg.ny} × {cfg.nz}, longitudinal: {cfg.nx} nodes, modes/valley: {cfg.modes}',
       f'- Minimum discrete lead bandwidth: {min_bandwidth_ev(cfg):.3f} eV',
       f'- Gate flatband adjustment: {args.gate_flatband:.3f} V (assumption, NOT calibrated)',
       '- HIGH-BIAS WARNING: 0–6 V with a 2 nm SiO2 is outside typical reliable oxide-operation conditions; breakdown, leakage and scattering are omitted.',
       '- The infinite hard-wall confinement and ideal source/drain contact models are assumptions.',
       '- Off-current may be sensitive to energy integration and contact/band alignment.',
       '', '## Plots', '']
    for filename in ['01_transverse_modes.png','02_idvgs.png','03_idvds.png','04_device_maps.png','05_scf_convergence.png',
                     '06_energy_convergence.png','07_mode_convergence.png','08_subbands_under_bias.png']:
        if (out/'figures'/filename).exists(): report.append(f'![{filename}](figures/{filename})')
    if summary is not None:
        report.extend(['','## Extracted Id–Vgs metrics','','```json',json.dumps(summary,indent=2),'```'])
    if triplet is not None:report.extend(['',f'Chosen [OFF, near ON, strong ON] Vgs: {triplet} V'])
    (out/'STATUS.md').write_text('\n'.join(report)+'\n')
    # Small reviewable images and Markdown automatically included in GitHub commits;
    # the bulky NPZ checkpoints stay in output/ and are gitignored.
    root=HERE/'figures';root.mkdir(exist_ok=True)
    for png in (out/'figures').glob('*.png'):
        shutil.copy2(png,root/png.name)
    (HERE/'STATUS.md').write_text('\n'.join(report)+'\n')
    metrics=out/'iv_metrics.json'
    if metrics.exists():shutil.copy2(metrics,root/'iv_metrics.json')
    quality=out/'quality_checks.json'
    if quality.exists():shutil.copy2(quality,root/'quality_checks.json')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',default='output')
    ap.add_argument('--phase',choices=['all','idvgs','idvds'],default='all')
    ap.add_argument('--vg-max',type=float,default=6.0)
    ap.add_argument('--vd-max',type=float,default=6.0)
    ap.add_argument('--step',type=float,default=0.2)
    ap.add_argument('--idvgs-vds',type=float,default=0.1)
    ap.add_argument('--ny',type=int,default=31)
    ap.add_argument('--nz',type=int,default=25)
    ap.add_argument('--nx',type=int,default=345)
    ap.add_argument('--modes',type=int,default=6)
    ap.add_argument('--energy-step',type=float,default=0.005)
    ap.add_argument('--max-iterations',type=int,default=100)
    ap.add_argument('--mix',type=float,default=0.20)
    ap.add_argument('--phi-tol',type=float,default=2e-3)
    ap.add_argument('--charge-tol',type=float,default=0.08)
    ap.add_argument('--current-tol',type=float,default=0.02)
    ap.add_argument('--gate-flatband',type=float,default=0.5,
        help='Assumed effective gate flatband offset V. Calibrate before interpreting Vth.')
    ap.add_argument('--manual-vth',type=float,default=None)
    ap.add_argument('--max-points',type=int,default=0,help='Stop after N NEW successful points (smoke tests)')
    ap.add_argument('--mesh-h',type=float,default=1.0,
        help='Gmsh mesh nm; use --regenerate-mesh when changing from STEP01')
    ap.add_argument('--regenerate-mesh',action='store_true')
    ap.add_argument('--allow-unconverged',action='store_true',help='Diagnostic only; do NOT claim physical validation')
    args=ap.parse_args()
    out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
    (out/'figures').mkdir(exist_ok=True);(out/'states').mkdir(exist_ok=True)
    cfg=replace(Config(),nx=args.nx,ny=args.ny,nz=args.nz,modes=args.modes,mesh_h_nm=args.mesh_h)
    cfg.validate()
    identity={'cfg':asdict(cfg),'voltage':{'vg_max':args.vg_max,
        'vd_max':args.vd_max,'step':args.step,'idvgs_vds':args.idvgs_vds},
        'numerics':{'energy_step':args.energy_step,'max_iterations':args.max_iterations,
         'mix':args.mix,'gate_flatband':args.gate_flatband,
         'phi_tol':args.phi_tol,'charge_tol':args.charge_tol,'current_tol':args.current_tol}}
    identity=json.loads(json.dumps(identity))
    meta=out/'run_settings.json'
    if meta.exists() and json.loads(meta.read_text())!=identity:
        raise SystemExit(f'Run settings changed relative to {meta}. Use a NEW --out directory '
                         'or exactly restore the original options to resume.')
    meta.write_text(json.dumps(identity,indent=2))
    if args.energy_step<=0 or not 0<args.mix<=1:ap.error('Energy step/mix must be positive')
    bandwidth=min_bandwidth_ev(cfg)
    print('=================================');print('STEP03 REAL DEVICE IV - NEGF SCF')
    print(f'DEVSIM 3D Poisson / Fixed-basis coupled MS-NEGF')
    print(f'EMA from STEP01, ny*nz={cfg.ny}*{cfg.nz}, nx={cfg.nx}, modes={cfg.modes}')
    print(f'dx={cfg.length_nm/(cfg.nx-1):.4f}nm; MIN lead bandwidth={bandwidth:.3f}eV')
    print(f'IdVgs 0..{args.vg_max:g} V step {args.step:g} V (Vds={args.idvgs_vds:g} V)')
    print(f'IdVds 0..{args.vd_max:g} V step {args.step:g} V (three Vgs picked after Vth)')
    print('=================================',flush=True)
    if bandwidth<args.vd_max+1.0:
        raise SystemExit(f'Discrete longitudinal lead bandwidth {bandwidth:.3f} eV too narrow for '
              f'{args.vd_max:g} V drain bias + safety margin. Increase --nx (smaller dx).')
    if args.vg_max>2 or args.vd_max>2:
        print('[PHYSICAL WARNING] 2nm oxide at 6V exceeds realistic operation; '
              'NO gate leakage/breakdown, inelastic scattering, or contact calibration.',flush=True)
    t0=time.perf_counter()
    bases=make_bases(cfg)
    make_modes_figure(cfg,bases,out/'figures'/'01_transverse_modes.png')
    x,y,z=axis_grid(cfg); axes=(x,y,z)
    # Default: generate a new *finer* 3D Poisson mesh for this physics stage.
    # Never overwrite the frozen STEP01 mesh.
    mesh=out/'mesh'/'mono_nanosheet_cm.msh'
    if args.regenerate_mesh or not mesh.exists():
        mesh=build_mesh(cfg,build_dir=out/'mesh')
    print('Mesh:',mesh,flush=True)
    p=DevsimPoisson(cfg,mesh)
    # initial charge approximately neutral source/drain; no DD invoked.
    xnodes=p.si_nodes_nm[:,0]
    p.set_charge(np.where((xnodes<cfg.ls_nm)|(xnodes>cfg.ls_nm+cfg.lg_nm),
                           cfg.donor_sd_cm3,cfg.donor_channel_cm3))
    new_points=0
    if args.phase in ('all','idvgs'):
        for vg in voltages(args.vg_max,args.step):
            key=f'vg_{vg:.2f}_vd_{args.idvgs_vds:.2f}'
            if already(out/'idvgs.csv',key):
                print(f'[RESUME] {key} already converged',flush=True)
                fill_state(p,out,'idvgs',key)
                continue
            fill_state(p,out,'idvgs',key)  # continue prior partial SCF if available
            try: res=run_bias(p,cfg,bases,float(vg),args.idvgs_vds,args,axes,out,key)
            except Exception as ex:
                print(f'[FAIL] {key}: {ex}',file=sys.stderr,flush=True)
                traceback.print_exc(); raise
            upsert_point(out/'idvgs.csv',res,
                ['key','vg_v','vd_v','id_a','status','iterations','delta_phi_v','charge_residual',
                 'relative_current_change','runtime_s','state_file'])
            if res['status']!='CONVERGED' and not args.allow_unconverged:
                write_report(out,cfg,args)
                raise SystemExit(f'STOP: SCF did not converge at {key}. Resume after correcting settings.')
            new_points+=1
            make_idvgs(out);make_scf_history(out);write_report(out,cfg,args)
            if args.max_points and new_points>=args.max_points:
                print('[PARTIAL] --max-points hit; no Vth extracted from incomplete Id-Vgs.');return
    rows=get_rows(out/'idvgs.csv')
    scan=voltages(args.vg_max,args.step)
    complete={round(float(r['vg_v']),7):r for r in rows if r['status']=='CONVERGED'}
    if not all(round(float(v),7) in complete for v in scan):
        raise SystemExit('Need fully converged Id-Vgs first; run --phase idvgs/--resume to complete.')
    ordered=[complete[round(float(v),7)] for v in scan]
    summary=summarize_iv(scan,[float(r['id_a']) for r in ordered],cfg.width_nm)
    (out/'iv_metrics.json').write_text(json.dumps(summary,indent=2))
    print('[VTH EXTRACTION]',summary['vth'],flush=True)
    if args.manual_vth is not None:
        print('[WARNING] manually supplied Vth override; NOT extracted from simulated Id-Vgs.',flush=True)
        vth=args.manual_vth
    else:vth=summary['vth']['vth_v']
    try:triplet=choose_vgs(vth,args.step,args.vg_max)
    except ValueError as exc:
        write_report(out,cfg,args,summary)
        raise SystemExit(f'Cannot select below/above/strong-on Vgs: {exc}')
    (out/'chosen_vgs.json').write_text(json.dumps({'extracted_vth_v':summary['vth']['vth_v'],
         'used_vth_v':vth,'manual_override':args.manual_vth is not None,
         'selected_vgs_v':triplet,'labels':['BELOW_VTH','JUST_ABOVE_VTH','STRONG_ON']},indent=2))
    print(f'[Id-Vds selected Vgs] {triplet} V',flush=True)
    # Fixed-potential transport convergence from the ACTUAL finished Id-Vgs point.
    vpick=min(scan,key=lambda v:abs(v-(triplet[1])))
    keypick=f'vg_{vpick:.2f}_vd_{args.idvgs_vds:.2f}'
    snap=state_path(out,keypick)
    if snap.exists():
        print(f'[QUALITY] Real SCF snapshot: {snap}',flush=True)
        data=np.load(snap)
        qstat=verify_snapshot(cfg,axes,data['phi_v'],args.idvgs_vds,out/'figures')
        print('[QUALITY] fixed-potential checks:',qstat,flush=True)
        if not qstat['energy_check_lt_3_percent'] or not qstat['mode_check_lt_5_percent']:
            print('[QUALITY WARNING] Current discretization not sufficiently converged. '
                  'Interpret extracted IV with caution.',flush=True)
    if args.phase=='idvgs':
        write_report(out,cfg,args,summary,triplet);return
    for vg in triplet:
        # New Vg family: reset to neutral donors rather than carrying prior ON charge.
        xnodes=p.si_nodes_nm[:,0]
        p.set_charge(np.where((xnodes<cfg.ls_nm)|(xnodes>cfg.ls_nm+cfg.lg_nm),
                               cfg.donor_sd_cm3,cfg.donor_channel_cm3))
        for vd in voltages(args.vd_max,args.step):
            key=f'vg_{vg:.2f}_vd_{vd:.2f}'
            if already(out/'idvds.csv',key):
                print(f'[RESUME] {key}',flush=True)
                fill_state(p,out,'idvds',key)
                continue
            fill_state(p,out,'idvds',key)  # continue prior partial SCF if available
            try:res=run_bias(p,cfg,bases,float(vg),float(vd),args,axes,out,key)
            except Exception as ex:
                print(f'[FAIL] {key}: {ex}',file=sys.stderr,flush=True)
                traceback.print_exc();raise
            upsert_point(out/'idvds.csv',res,
                ['key','vg_v','vd_v','id_a','status','iterations','delta_phi_v','charge_residual',
                 'relative_current_change','runtime_s','state_file'])
            if res['status']!='CONVERGED' and not args.allow_unconverged:
                write_report(out,cfg,args,summary,triplet)
                raise SystemExit(f'STOP: SCF did not converge at {key}. Resume after correcting settings.')
            new_points+=1
            make_idvds(out);make_scf_history(out)
            if vd in (0,round(args.vd_max/2,7),args.vd_max):make_device_maps(out,state_path(out,key))
            write_report(out,cfg,args,summary,triplet)
            if args.max_points and new_points>=args.max_points:
                print('[PARTIAL] --max-points hit. Restart with same --out.');return
    write_report(out,cfg,args,summary,triplet)
    print(f'[DONE] IV data and plots saved to {out}; elapsed {time.perf_counter()-t0:.1f}s',flush=True)

if __name__=='__main__':main()
