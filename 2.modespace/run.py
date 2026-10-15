#!/usr/bin/env python3
"""STEP 02 | Mode-space transverse Schrödinger verification.

Reuse original EMA implementation under ../1.example_build.  This step does NOT
modify its effective masses, Hamiltonian, or NEGF code.

Modes:
  python run.py                           # baseline (flat phi, unless Poisson found)
  python run.py --poisson ../1.example_build/build/result/poisson_only.npz
  python run.py --poisson none             # deliberately flat reference
  python -m unittest discover -s tests -v
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import json
import shutil
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).resolve().parent
DEFAULT_BASE = HERE.parent / '1.example_build'


def load_reference(base: Path):
    """Load the SAME EMA config/Schrodinger implementation from STEP 01."""
    base = base.resolve()
    if not (base / 'config.py').exists() or not (base / 'negf' / 'mode_space.py').exists():
        raise FileNotFoundError(
            f'Reference not found: {base}. Store the previous code under '
            '1.example_build/ (config.py and negf/mode_space.py).')
    sys.path.insert(0, str(base))
    from config import Config  # noqa: E402
    from negf.mode_space import (HB2_OVER_2M0, build_mode_slices,
                                 transverse_hamiltonian)  # noqa: E402
    return Config, HB2_OVER_2M0, build_mode_slices, transverse_hamiltonian


def axes(cfg):
    x = np.linspace(0, cfg.length_nm, cfg.nx)
    y = np.linspace(-cfg.width_nm / 2, cfg.width_nm / 2, cfg.ny + 2)[1:-1]
    z = np.linspace(-cfg.thickness_nm / 2, cfg.thickness_nm / 2, cfg.nz + 2)[1:-1]
    return x, y, z


def load_phi(path: Path | None, cfg):
    x, y, z = axes(cfg)
    if path is None:
        return np.zeros((len(x), len(y), len(z))), 'FLAT analytical reference (NOT a device potential)'
    data = np.load(path)
    expected = {'phi_v', 'x_nm', 'y_nm', 'z_nm'}
    if not expected.issubset(data.files):
        raise ValueError(f'{path}: missing arrays {sorted(expected.difference(data.files))}')
    grid = tuple(np.asarray(data[key], dtype=float) for key in ('x_nm', 'y_nm', 'z_nm'))
    phi = np.asarray(data['phi_v'], dtype=float)
    if phi.shape != tuple(len(a) for a in grid):
        raise ValueError(f'Poisson array shape {phi.shape} incompatible with axes')
    for a in grid:
        if len(a) < 2 or not np.all(np.diff(a) > 0):
            raise ValueError('Poisson axes must be strictly increasing')
    xi = np.stack(np.meshgrid(x, y, z, indexing='ij'), -1)
    interp = RegularGridInterpolator(grid, phi, bounds_error=True)
    out = interp(xi)
    if not np.isfinite(out).all():
        raise ValueError('Nonfinite Poisson potential')
    return out, f'DEVSIM Poisson: {path}'


def sorted_analytic_levels(cfg, valley, const):
    """Exact eigenenergies of the *discrete* 2D Dirichlet FD operator.

    This is a stronger implementation check than comparing to the continuum
    solution, which contains finite-difference truncation error.
    """
    dy = cfg.width_nm / (cfg.ny+1)
    dz = cfg.thickness_nm / (cfg.nz+1)
    ty = const / (valley.my * dy * dy)
    tz = const / (valley.mz * dz * dz)
    vals = []
    for p in range(1, cfg.ny+1):
        for q in range(1, cfg.nz+1):
            d = 2*ty*(1-np.cos(np.pi*p/(cfg.ny+1)))
            d += 2*tz*(1-np.cos(np.pi*q/(cfg.nz+1)))
            vals.append(d + cfg.ec0_ev)
    return np.sort(np.array(vals))


def continuum_ground(cfg, valley, const):
    return cfg.ec0_ev + np.pi**2 * const * (
        1.0/(valley.my*cfg.width_nm**2) +
        1.0/(valley.mz*cfg.thickness_nm**2))


def lowest_eigenpairs(h, nm):
    n = h.shape[0]
    if n <= 140:
        es, us = eigh(h, subset_by_index=(0, nm-1), check_finite=False)
    else:
        es, us = eigsh(h, k=nm, which='SA', tol=1e-11, maxiter=20000)
        idx = np.argsort(es)
        es, us = es[idx], us[:,idx]
    return es, us


def benchmark_grid(cfg, transverse_hamiltonian, const, levels=4):
    rows = []
    # Converged 2D discretization needed to make the 5nm confinement credible.
    sizes = [(7,5), (11,9), (19,15), (31,25), (47,39)]
    for ny, nz in sizes:
        c = replace(cfg, ny=ny, nz=nz)
        for valley in c.valleys:
            h = transverse_hamiltonian(c, valley, np.zeros((ny,nz)))
            e, u = lowest_eigenpairs(h, min(levels, ny*nz-1))
            e_exact_fd = sorted_analytic_levels(c, valley, const)
            e0_cont = continuum_ground(c, valley, const)
            err_fd = np.max(np.abs(e - e_exact_fd[:len(e)]))
            err_cont = abs(e[0]-e0_cont)/max(abs(e0_cont-c.ec0_ev), 1e-20)
            rows.append(dict(ny=ny,nz=nz,valley=valley.name,
                             e1_ev=float(e[0]), continuum_e1_ev=float(e0_cont),
                             relative_ground_error=float(err_cont),
                             discrete_analytic_max_error_ev=float(err_fd),
                             orthogonality_error=float(np.max(np.abs(u.conj().T@u-np.eye(len(e)))))))
    return rows


def save_rows(path, rows):
    with path.open('w', newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def style_axis(ax, xlabel=None, ylabel=None):
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    ax.grid(alpha=.22,linewidth=.6)


def render(cfg, phi, source, results, bench, out, const):
    out.mkdir(parents=True,exist_ok=True)
    x, y, z = axes(cfg)
    dy = cfg.width_nm/(cfg.ny+1)
    dz = cfg.thickness_nm/(cfg.nz+1)
    mid = cfg.nx // 2
    valley = cfg.valleys[0]
    e, u, _, _, _ = results[valley.name]

    # Plot A: 2D |psi|^2 at the channel center, full 2D transverse confinement.
    count = min(4, cfg.modes)
    fig, axs = plt.subplots(1,count,figsize=(3.7*count,3.5),squeeze=False)
    for k, ax in enumerate(axs[0]):
        prob = np.abs(u[mid,:,k].reshape(cfg.ny,cfg.nz))**2/(dy*dz)
        im=ax.imshow(prob.T, origin='lower', aspect='auto', interpolation='nearest',
                     extent=[-cfg.width_nm/2+dy/2,cfg.width_nm/2-dy/2,
                             -cfg.thickness_nm/2+dz/2,cfg.thickness_nm/2-dz/2])
        ax.set_title(f'Mode {k+1}: E={e[mid,k]:.4f} eV')
        ax.set_xlabel('y (nm)')
        if k==0: ax.set_ylabel('z (nm)')
        fig.colorbar(im,ax=ax,shrink=.82,label=r'$|\chi|^2$ (nm$^{-2}$)')
    fig.suptitle(f'Transverse confinement: {valley.name}  |  x={x[mid]:g} nm\n{source}',fontsize=11)
    fig.tight_layout()
    fig.savefig(out/'01_transverse_modes.png',dpi=170)
    plt.close(fig)

    # Plot B: local subbands along transport axis. Demonstrates gate-sensitive Ec(x).
    fig,ax=plt.subplots(figsize=(9,4.7))
    for name,(evals,_,_,_,_) in results.items():
        for m in range(min(3,cfg.modes)):
            ax.plot(x,evals[:,m],lw=1.65,label=f'{name} / mode {m+1}')
    for xp in (cfg.ls_nm,cfg.ls_nm+cfg.lg_nm):
        ax.axvline(xp,color='gray',linestyle=':',lw=1)
    ax.axvspan(cfg.ls_nm,cfg.ls_nm+cfg.lg_nm,color='gray',alpha=.06)
    ax.set_title('Local transverse subband energies along transport  (NOT current)')
    style_axis(ax,'x, transport direction (nm)','Subband energy (eV)')
    ax.legend(loc='best',fontsize=8,ncol=3)
    fig.tight_layout()
    fig.savefig(out/'02_subband_vs_x.png',dpi=170)
    plt.close(fig)

    # Plot C: rigorous 2D grid convergence vs continuum infinite well.
    fig,ax=plt.subplots(figsize=(8.5,4.8))
    for valley in cfg.valleys:
        sub=[r for r in bench if r['valley']==valley.name]
        ax.plot([r['nz'] for r in sub],[100*r['relative_ground_error'] for r in sub],
                'o-',label=f'{valley.name}')
    ax.axhline(1.0,color='gray',ls='--',lw=1,label='1% reference target')
    ax.set_yscale('log')
    ax.set_title('Grid convergence of ground-state confinement energy')
    style_axis(ax,'Number of interior grid points across thickness (z)','Error vs continuum, %')
    ax.legend()
    fig.tight_layout()
    fig.savefig(out/'03_grid_convergence.png',dpi=170)
    plt.close(fig)

    # Plot D: onsite orthogonality and adjacent-slice overlaps.
    u0=u[mid]
    S=u0.conj().T@u0
    if mid+1 < cfg.nx: Uchange=u0.conj().T@u[mid+1]
    else: Uchange=u[mid-1].conj().T@u0
    fig,axs=plt.subplots(1,2,figsize=(8.5,3.9))
    for ax,mat,tit in [(axs[0],np.abs(S),r'$|U_x^\dagger U_x|$: onsite'),
                       (axs[1],np.abs(Uchange),r'$|U_x^\dagger U_{x+1}|$: adjacent slices')]:
        im=ax.imshow(mat,vmin=0,vmax=1,cmap='viridis')
        ax.set_title(tit)
        ax.set_xlabel('Mode index')
        ax.set_ylabel('Mode index')
        ax.set_xticks(range(cfg.modes),range(1,cfg.modes+1))
        ax.set_yticks(range(cfg.modes),range(1,cfg.modes+1))
        fig.colorbar(im,ax=ax,shrink=.78)
    fig.suptitle(f'Mode orthogonality and local basis overlap: {valley.name}')
    fig.tight_layout()
    fig.savefig(out/'04_mode_overlap.png',dpi=170)
    plt.close(fig)

    # Plot E: gate potential profile on the centerline, to distinguish physics from test.
    fig,ax=plt.subplots(figsize=(8.5,4.4))
    ax.plot(x,phi[:,cfg.ny//2,cfg.nz//2],lw=2,label='Si centerline')
    ax.plot(x,phi[:,0,0],lw=1.5,label='Si near edge')
    for xp in (cfg.ls_nm,cfg.ls_nm+cfg.lg_nm):
        ax.axvline(xp,color='gray',ls=':',lw=1)
    ax.set_title('Electrostatic input to the transverse Hamiltonian')
    style_axis(ax,'x (nm)','Potential (V)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(out/'05_input_potential.png',dpi=170)
    plt.close(fig)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',type=Path,default=DEFAULT_BASE,
                        help='Path to untouched STEP 01 code (config.py + negf/)')
    parser.add_argument('--poisson',default='auto',
                        help='Poisson *.npz from STEP01; auto=use if available, none=flat reference')
    parser.add_argument('--modes',type=int,default=4)
    parser.add_argument('--out',type=Path,default=HERE/'results')
    args=parser.parse_args(argv)
    Config, const, builder, transverse_hamiltonian=load_reference(args.base)
    cfg=Config()
    cfg.modes=min(args.modes,cfg.ny*cfg.nz-1)
    cfg.validate()

    default_npz=args.base/'build/result/poisson_only.npz'
    if args.poisson=='auto':
        path=default_npz if default_npz.is_file() else None
    elif args.poisson.lower()=='none':
        path=None
    else:
        path=Path(args.poisson).resolve()
        if not path.exists(): raise FileNotFoundError(path)
    phi,source=load_phi(path,cfg)
    results={}
    for valley in cfg.valleys:
        results[valley.name]=builder(cfg,valley,phi)

    bench=benchmark_grid(cfg,transverse_hamiltonian,const)
    out=args.out.resolve()
    out.mkdir(parents=True,exist_ok=True)
    save_rows(out/'grid_convergence.csv',bench)
    np.savez_compressed(out/'modes.npz', x_nm=axes(cfg)[0],y_nm=axes(cfg)[1],z_nm=axes(cfg)[2],
                        phi_v=phi, **{f'E_{name}_ev':a[0] for name,a in results.items()},
                        **{f'U_{name}':a[1] for name,a in results.items()},
                        **{f'overlap_hop_{name}':a[4] for name,a in results.items()})
    render(cfg,phi,source,results,bench,out,const)

    # Validation metrics: *numerical* validity, not device calibration.
    check={}
    for name,(es,us,t,onsite,hop) in results.items():
        norms=np.array([np.max(np.abs(us[i].T@us[i]-np.eye(cfg.modes))) for i in range(cfg.nx)])
        herm=max(float(np.max(np.abs(mat-mat.T.conj()))) for mat in onsite)
        check[name]={'E1_source_ev':float(es[0,0]),
                     'E1_channel_ev':float(es[cfg.nx//2,0]),
                     'onsite_orthogonality_max':float(np.max(norms)),
                     'onsite_hermiticity_error':herm,
                     'modes':cfg.modes,
                     'hop_shape':list(hop.shape)}
    best=[r for r in bench if r['nz']==39]
    worst_discrete=max(r['discrete_analytic_max_error_ev'] for r in bench)
    worst_orth=max(v['onsite_orthogonality_max'] for v in check.values())
    worst_cont=max(r['relative_ground_error'] for r in best)
    report={'step':'02 Transverse Schrödinger / mode-space basis',
            'input_source':source,'potential_is_self_consistent':False,
            'no_ema_parameters_changed':True,
            'note':'STEP02 tests basis properties only. NEGF currents and SCF are STEP03/04.',
            'grid_analytic_error_max_ev':worst_discrete,
            'onsite_overlap_max_error':worst_orth,
            'finest_grid_continuum_relative_error_max':worst_cont,
            'checks':check,
            'acceptance':{
                'exact_2d_fd_eigenvalues':worst_discrete < 2e-8,
                'orthonormal_modes':worst_orth < 1e-10,
                'continuum_ground_under_1_percent':worst_cont < .01,
                'finite_values':bool(np.isfinite(phi).all() and all(np.isfinite(v[0]).all() for v in results.values()))
            }}
    report['step02_pass']=all(report['acceptance'].values())
    (out/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')

    # GitHub-reviewable figures. Keep bulky numerical arrays in ignored results/.
    figures = HERE / 'figures'
    figures.mkdir(parents=True, exist_ok=True)
    for png in sorted(out.glob('*.png')):
        shutil.copy2(png, figures/png.name)
    for name in ('validation.json', 'grid_convergence.csv'):
        shutil.copy2(out/name, figures/name)
    state = 'PASS' if report['step02_pass'] else 'FAIL'
    (figures/'STATUS.md').write_text(
        '# STEP 02 — ' + state + '\n\n'
        + '**Input potential:** ' + source + '\n\n'
        + '**Scope:** 2D transverse mode-space numerical validation only; '
          'no device NEGF current is calculated.\n\n'
        + '| Check | Outcome |\n|---|---|\n'
        + ''.join('| ' + k + ' | ' + ('PASS' if v else 'FAIL') + ' |\n'
                  for k,v in report['acceptance'].items())
        + '\n'
        + ''.join('![%s](%s)\n\n' % (png.stem,png.name)
                  for png in sorted(figures.glob('*.png')))
        ,encoding='utf-8')
    print('\nSTEP 02  |  TRANSVERSE SCHRODINGER / MODE BASIS')
    print('Potential source:', source)
    for k,v in check.items():
        print(f"{k:>8}  E1(source)={v['E1_source_ev']:.6f} eV  "
              f"E1(channel)={v['E1_channel_ev']:.6f} eV  "
              f"orthogonality={v['onsite_orthogonality_max']:.2e}")
    for k,v in report['acceptance'].items():
        print(f"{'PASS' if v else 'FAIL'}  {k}")
    print('STEP02 RESULT:', 'PASS' if report['step02_pass'] else 'FAIL')
    print('Plots and CSV:',out)
    print('GitHub figures:',figures)
    if not report['step02_pass']: raise SystemExit(2)


if __name__=='__main__': main()
