#!/usr/bin/env python3
"""Actual GAA frozen-field contact comparison, NOT a semi-infinite SCF result."""
import hashlib
import json
import sys
import time
import traceback
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from config import Config,VALLEYS,KB_EV,Q,HBAR_EV_S
from modespace.hamiltonian import Grid,ModeBuilder
from multispace.occupation import fermi
from transport.negf_reference import device_blocks,lead_matrices,rgf,neutral_lead_mu
from transport.quadrature import integrate_vector
from validation.post_scf_current import finite_rgf_transmission


def run(directory,key,out):
    directory=Path(directory);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    settings=json.loads((directory/'settings.json').read_text());c=Config(**settings['config'])
    with np.load(directory/'states'/(key+'.npz')) as saved:
        if not bool(saved['converged']):raise ValueError('Use an accepted finite-bath state')
        phi=saved['phi_v'].copy();mu=float(saved['mu_l_ev']);vd=float(saved['vd_v'])
        grid=Grid.build(c);central=(grid.x[:-1]>=c.ls-c.dx/2)&(grid.x[1:]<=c.drain_start+c.dx/2)
        finite_current=float(saved['bond_current_a'][central].mean())
    builders=[ModeBuilder(c,grid,valley) for valley in VALLEYS]
    systems=[builder.build(phi) for builder in builders]
    opened=[device_blocks(system,c,(0.,vd)) for system in systems]
    neutral_mu,calibration=neutral_lead_mu(builders,c,grid)
    def transmission(energies):
        total=np.zeros(len(energies))
        for system,(onsite,hop,tx,leads) in zip(systems,opened):
            sl,sr=lead_matrices(np.asarray(energies,complex),tx,leads)
            _,amplitude=rgf(energies,onsite,hop,sl,sr,False,0)
            gamma=1j*(sr-sr.conj().transpose(0,2,1));a=amplitude[-1]
            t=np.einsum('eij,ejk,eik->e',gamma,a,a.conj(),optimize=True).real
            if t.min()<-1e-10:raise RuntimeError('Negative semi-infinite transmission')
            total+=2*system.valley.multiplicity*np.maximum(t,0)
        return total
    chemical_potentials=np.array([mu,neutral_mu]);kt=KB_EV*c.temperature
    low=min(mu,neutral_mu,mu-vd,neutral_mu-vd)-20*kt
    high=max(mu,neutral_mu,mu-vd,neutral_mu-vd)+20*kt
    edges=np.concatenate([np.r_[e-2*tx,e+2*tx] for onsite,hop,tx,leads in opened for e,u in leads])
    edges=edges[(edges>low)&(edges<high)]
    def integrand(energy):
        df=fermi(energy[:,None],chemical_potentials[None,:],c.temperature)-fermi(energy[:,None],chemical_potentials[None,:]-vd,c.temperature)
        return transmission(energy)[:,None]*df
    checks=[];start=time.perf_counter()
    for panel in [.02,.01,.005]:
        breaks=np.unique(np.r_[np.arange(low,high,panel),high,edges])
        integral,diagnostics=integrate_vector(integrand,breaks,8,2e-7,1e-16,100000)
        current=integral.real*Q/HBAR_EV_S/(2*np.pi)
        checks.append({'initial_panel_ev':panel,'same_mu_current_a':float(current[0]),
                       'neutral_lead_mu_current_a':float(current[1]),**diagnostics})
    energy=np.linspace(low,high,800);semi_t=transmission(energy)
    finite_t=sum(2*system.valley.multiplicity*finite_rgf_transmission(system,c.gamma,energy) for system in systems)
    np.savez_compressed(out/'transmission.npz',energy_ev=energy,
        semiinfinite_spin_valley_summed_transmission=semi_t,finite_bath_spin_valley_summed_transmission=finite_t)
    result={'status':'COMPLETED_FROZEN_CURRENT_ONLY','frozen_potential_only':True,
            'semiinfinite_quantum_poisson_scf':False,'density_comparison_performed':False,
            'quantitative_validation_pass':False,'saved_state':str(directory/'states'/(key+'.npz')),
            'saved_state_source_sha256':settings['source_sha256'],
            'comparison_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'finite_bath_mu_ev':mu,'infinite_neutral_lead_mu_ev':neutral_mu,
            'lead_calibration':calibration,'finite_bath_bond_current_a':finite_current,
            'checks':checks,'runtime_s':time.perf_counter()-start,
            'model_difference':'Same saved potential, masses, grid, retained transverse space and interior FD projection. Semi-infinite continuation replaces only terminal x hard walls and finite distributed baths. The second chemical potential comparison is separately calibrated, never a fitted current.',
            'limitation':'Current quadrature on a frozen finite-bath SCF field does not establish semi-infinite density, charge neutrality or coupled SCF convergence; a quadrature estimator can miss exceptionally narrow peaks.'}
    (out/'summary.json').write_text(json.dumps(result,indent=2));return result


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('directory');parser.add_argument('key')
    parser.add_argument('--out',default='results/frozen_contact_comparison');args=parser.parse_args()
    try:print(json.dumps(run(args.directory,args.key,args.out),indent=2))
    except Exception:
        out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
        (out/'failure.txt').write_text(traceback.format_exc())
        (out/'summary.json').write_text(json.dumps({'status':'ERROR','frozen_potential_only':True,'quantitative_validation_pass':False,'error_file':'failure.txt'},indent=2))
        raise
