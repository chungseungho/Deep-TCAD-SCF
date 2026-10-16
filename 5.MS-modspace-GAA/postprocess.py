#!/usr/bin/env python3
"""Actual converged-device transmission quadrature convergence, after SCF only."""
import argparse,json,csv
from pathlib import Path
import numpy as np
from config import Config,VALLEYS,Q,HBAR_EV_S
from modespace.hamiltonian import Grid,ModeBuilder
from transport.spectral import finite_transmission
from multispace.occupation import fermi


def analyze(directory,key):
    out=Path(directory);c=Config(**json.loads((out/'settings.json').read_text())['config'])
    s=np.load(out/'states'/f'{key}.npz')
    if not bool(s['converged']): raise ValueError('Cannot postprocess a failed SCF state as accepted')
    if c.backend!='spectral': raise ValueError('This postprocessor checks the non-Markovian finite spectral backend')
    grid=Grid.build(c);mu=float(s['mu_l_ev']);mr=float(s['mu_r_ev'])
    systems=[ModeBuilder(c,grid,v).build(s['phi_v']) for v in VALLEYS]
    currents=s['bond_current_a'];select=(grid.x[:-1]>=c.ls-c.dx/2)&(grid.x[1:]<=c.drain_start+c.dx/2)
    reference=float(np.mean(currents[select]));rows=[]
    for step in [.002,.001,.0005]:
        E=np.arange(min(mu,mr)-.5,max(mu,mr)+.8+step*.1,step);total=np.zeros(len(E))
        for sys in systems:
            total+=2*sys.valley.multiplicity*finite_transmission(sys.H,sys.regions,c.gamma,E)
        J=Q/HBAR_EV_S/(2*np.pi)*total*(fermi(E,mu,c.temperature)-fermi(E,mr,c.temperature))
        integrated=float(np.trapezoid(J,E));row={'step_ev':step,'current_a':integrated,'bond_current_a':reference,
            'relative_mismatch':abs(integrated-reference)/max(abs(reference),1e-20)};rows.append(row)
        print(row,flush=True)
    # Independent full matrix inversion at five energies (not this pole expansion).
    direct_error=0.
    samples=np.linspace(mr-.05,mu+.15,5)
    for sys in systems:
        pole=finite_transmission(sys.H,sys.regions,c.gamma,samples)
        i=np.where(sys.regions==0)[0];j=np.where(sys.regions==2)[0]
        for k,E0 in enumerate(samples):
            G=np.linalg.solve(E0*np.eye(len(sys.H))-sys.H+.5j*np.diag((sys.regions!=1)*c.gamma),np.eye(len(sys.H)))
            T=c.gamma**2*np.sum(abs(G[np.ix_(i,j)])**2)
            direct_error=max(direct_error,abs(T-pole[k]))
    np.savez_compressed(out/'states'/f'{key}_spectra.npz',energy_ev=E,transmission=total,
                        spectral_current_a_per_ev=J,integrated_current_a=integrated)
    report={'key':key,'energy_quadrature':rows,'direct_vs_pole_transmission_max_abs_error':float(direct_error),
            'energy_window_ev':[float(E[0]),float(E[-1])],
            'note':'Spectral convergence at this saved bias only; not device mesh/contact convergence'}
    (out/'spectral_convergence.json').write_text(json.dumps(report,indent=2))
    from plotting import plot_run
    plot_run(out);return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory');p.add_argument('key');a=p.parse_args()
    analyze(a.directory,a.key)
