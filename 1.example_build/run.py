#!/usr/bin/env python3
"""Mono-Nanosheet GAA: Gmsh + DEVSIM Poisson + ballistic MS-NEGF reference.

Examples:
  python run.py --mode mesh
  python run.py --mode poisson --vg 0.4 --vd 0.05
  python run.py --mode negf-only --vg 0.4 --vd 0.05
  python run.py --mode coupled --vg 0.4 --vd 0.05 --iterations 8

The 'negf-only' mode uses an explicitly approximate analytical potential,
intended only for smoke tests of the transport engine without DEVSIM.
"""
import argparse
import csv
from pathlib import Path
import sys
import time

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from config import Config
from geometry.make_mesh import build_mesh
from negf.mode_space import evaluate_negf


def transport_grid(cfg):
    x = np.linspace(0, cfg.length_nm, cfg.nx)
    y = np.linspace(-cfg.width_nm/2, cfg.width_nm/2, cfg.ny+2)[1:-1]
    z = np.linspace(-cfg.thickness_nm/2, cfg.thickness_nm/2, cfg.nz+2)[1:-1]
    return x,y,z


def analytical_start_phi(cfg, vg, vd):
    """Uncalibrated test profile: gate controls center, contact potentials at ends."""
    x,_,_ = transport_grid(cfg)
    src_drn = vd*(x/cfg.length_nm)
    gate_window = np.clip((x-cfg.ls_nm)/3.0,0,1) * np.clip((cfg.ls_nm+cfg.lg_nm-x)/3.0,0,1)
    return np.broadcast_to((src_drn + vg*gate_window)[:,None,None],
                           (cfg.nx,cfg.ny,cfg.nz)).copy()


def initial_donor_field(cfg, nodes_nm):
    x = nodes_nm[:,0]
    return np.where((x < cfg.ls_nm) | (x > cfg.ls_nm+cfg.lg_nm),
                    cfg.donor_sd_cm3, cfg.donor_channel_cm3)


def save_outputs(outdir, cfg, x, y, z, phi, result, history, vg, vd, label):
    outdir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(outdir/'solution.npz', x_nm=x, y_nm=y, z_nm=z,
                        phi_v=phi, electron_density_cm3=result.electron_density_cm3,
                        energy_ev=result.spectrum_ev, transmission=result.transmission,
                        gate_voltage_v=vg, drain_voltage_v=vd)
    with (outdir/'history.csv').open('w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['iteration','Id_A','max_delta_phi_V','mean_electron_cm-3','elapsed_s'])
        w.writerows(history)
    with (outdir/'result.txt').open('w') as f:
        f.write(f'MODE={label}\nVgs={vg:.6g} V\nVds={vd:.6g} V\n')
        f.write(f'Id={result.current_a:.8e} A\n')
        f.write(f'NEGF per valley={result.valley_currents_a}\n')
        f.write('WARNING: Prototype reference calculation, not calibrated or validated.\n')
        f.write('No inelastic scattering; Dirichlet hard-wall confinement; source/drain contacts ideal.\n')
        f.write('See README for fixed electrochemical reference and limitations.\n')
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axs=plt.subplots(1,2,figsize=(10,4))
        midz=cfg.nz//2
        im=axs[0].pcolormesh(x,y,phi[:,:,midz].T,shading='auto')
        axs[0].set(xlabel='x (nm)',ylabel='y (nm)',title='Potential (V) at z~0')
        fig.colorbar(im,ax=axs[0])
        axs[1].plot(result.spectrum_ev,result.transmission)
        axs[1].set(xlabel='Energy (eV)',ylabel='Transmission (per spin)',title='Mode-space NEGF')
        axs[1].grid(alpha=.2)
        fig.tight_layout()
        fig.savefig(outdir/'overview.png',dpi=170)
        plt.close(fig)
    except ImportError:
        pass
    print(f'[saved] {outdir}/solution.npz, result.txt, history.csv, overview.png')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=['mesh','poisson','negf-only','coupled'],default='coupled')
    parser.add_argument('--vg',type=float,default=0.4)
    parser.add_argument('--vd',type=float,default=0.05)
    parser.add_argument('--iterations',type=int,default=8)
    parser.add_argument('--mix',type=float,default=0.15)
    parser.add_argument('--mesh-h',type=float,default=2.0)
    parser.add_argument('--modes',type=int,default=2)
    parser.add_argument('--nx',type=int,default=49)
    parser.add_argument('--energy-points',type=int,default=220)
    parser.add_argument('--reuse-mesh',action='store_true')
    parser.add_argument('--output',default='build/result')
    args=parser.parse_args()
    cfg=Config(mix=args.mix,mesh_h_nm=args.mesh_h,nx=args.nx,modes=args.modes,
               energy_points=args.energy_points)
    cfg.validate()
    x,y,z=transport_grid(cfg)
    if args.mode=='negf-only':
        phi=analytical_start_phi(cfg,args.vg,args.vd)
        result=evaluate_negf(cfg,phi,args.vg,args.vd)
        save_outputs(Path(args.output),cfg,x,y,z,phi,result,[],args.vg,args.vd,
                     'TEST ONLY: analytical electrostatics, no Poisson')
        print(f'NEGF-only Id = {result.current_a:.8e} A')
        return
    mesh_file=Path('build/mono_nanosheet_cm.msh')
    if not args.reuse_mesh or not mesh_file.exists():
        mesh_file=build_mesh(cfg,build_dir='build')
    if args.mode=='mesh':
        print('Use: gmsh build/mono_nanosheet_nm.msh')
        return
    from poisson.devsim_poisson import DevsimPoisson
    p=DevsimPoisson(cfg,mesh_file)
    # Zero-space-charge initial condition, matching fixed donors initially.
    p.set_charge(initial_donor_field(cfg,p.si_nodes_nm))
    start=time.perf_counter()
    p.solve(args.vg,args.vd)
    phi=p.sample_si_potential(x,y,z)
    if args.mode=='poisson':
        from devsim import write_devices
        Path(args.output).mkdir(parents=True, exist_ok=True)
        np.savez_compressed(Path(args.output)/'poisson_only.npz',phi_v=phi,x_nm=x,y_nm=y,z_nm=z)
        p.export_vtk(Path(args.output)/'poisson_only.vtk')
        print('Poisson solved; no NEGF run. Saved:',Path(args.output)/'poisson_only.npz')
        return
    history=[]
    old_current=None
    result=None
    for i in range(1,args.iterations+1):
        result=evaluate_negf(cfg,phi,args.vg,args.vd, density=True)
        if not np.all(np.isfinite(result.electron_density_cm3)):
            raise RuntimeError('Nonfinite NEGF electron density')
        # Mixed density, not mixed potential: external density fed to Poisson.
        new_nodes = RegularGridInterpolator(
            (x,y,z),result.electron_density_cm3,method='linear',
            bounds_error=False,fill_value=0)(p.si_nodes_nm)
        old_nodes=np.asarray(p.d.get_node_model_values(device=p.device,
                                                      region=p.si,name='InjectedElectrons'))
        p.set_charge((1-cfg.mix)*old_nodes + cfg.mix*new_nodes)
        p.solve(args.vg,args.vd)
        new_phi=p.sample_si_potential(x,y,z)
        dv=float(np.max(np.abs(new_phi-phi)))
        elapsed=time.perf_counter()-start
        history.append((i,result.current_a,dv,float(np.mean(result.electron_density_cm3)),elapsed))
        print(f'[SCF {i:02d}] Id={result.current_a:.5e} A  |Δphi|max={dv:.4e} V  t={elapsed:.1f}s',flush=True)
        phi=new_phi
        if dv<cfg.tol_v and old_current is not None:
            print('[converged] electrostatic potential tolerance reached')
            break
        old_current=result.current_a
    # Final transport on FINAL Poisson potential, not on the previous iteration.
    result=evaluate_negf(cfg,phi,args.vg,args.vd)
    last_dv=history[-1][2]
    if last_dv>=cfg.tol_v:
        print(f'[not converged] last Δphi={last_dv:.3g} V; increase iterations and/or reduce --mix',
              file=sys.stderr)
    save_outputs(Path(args.output),cfg,x,y,z,phi,result,history,args.vg,args.vd,
                 '3D DEVSIM Poisson + 2D modes + coupled MS-NEGF')
    try:
        p.export_vtk(Path(args.output)/'potential_3d.vtk')
    except Exception as ex:
        print('Optional VTK export failed:',ex,file=sys.stderr)
    print(f'Final Id={result.current_a:.8e} A; walltime={time.perf_counter()-start:.1f}s')


if __name__=='__main__':
    main()
