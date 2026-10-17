"""Number intervals before any unresolved origin is assigned a population.

For each unclassified eigenstate, min(fL,fR) <= f <= max(fL,fR).
The resulting monotone N_low/high define a feasible mu interval. A midpoint
is used ONLY to audit uncertainty; materially ambiguous cases return no rho.
Tail-only accepted populations retain an explicitly bounded count residual.
This is a finite electron constraint, not a variational DFT functional.
"""
import numpy as np
from scipy.optimize import brentq
from multispace.occupation import fermi,occupations,AmbiguousOccupation


def number_bounds(states,degeneracies,mu,vd,temperature):
    lo=hi=0.
    for s,g in zip(states,degeneracies):
        fl=fermi(s['e'],mu,temperature);fr=fermi(s['e'],mu-vd,temperature)
        known=np.where(s['origin']==2,fr,fl);unknown=s['origin']==-1
        lo+=g*np.where(unknown,np.minimum(fl,fr),known).sum()
        hi+=g*np.where(unknown,np.maximum(fl,fr),known).sum()
    return float(lo),float(hi)


def solve_number(states,degeneracies,target,vd,c,fixed_mu):
    bounds=lambda mu:number_bounds(states,degeneracies,mu,vd,c.temperature)
    if c.ms_number_constraint=='finite_total':
        if not 0<target<sum(len(s['e'])*g for s,g in zip(states,degeneracies)):
            raise ValueError('Not enough states for electron count')
        lower=min(s['e'].min() for s in states)-abs(vd)-5
        upper=max(s['e'].max() for s in states)+abs(vd)+5
        mu_lo=brentq(lambda m:bounds(m)[1]-target,lower,upper,xtol=1e-13)
        mu_hi=brentq(lambda m:bounds(m)[0]-target,lower,upper,xtol=1e-13)
        mu=.5*(mu_lo+mu_hi)
    else:mu_lo=mu_hi=mu=float(fixed_mu)
    nlo,nhi=bounds(mu)
    diagnostic={'constraint':c.ms_number_constraint,'target_electrons':float(target),
        'chemical_potential_interval_ev':[float(mu_lo),float(mu_hi)],
        'chemical_potential_interval_width_ev':float(mu_hi-mu_lo),
        'electron_number_interval_at_audit_mu':[nlo,nhi],
        'electron_number_uncertainty_bound':float(nhi-nlo),
        'audit_mu_ev':float(mu),'unresolved_assignment_in_number_equation':False,
        'bias_constraint_error_ev':float(mu-(mu-vd)-vd),'variational_functional_implemented':False}
    filled=[occupations(s,mu,mu-vd,c.temperature,c.ms_occupation_uncertainty,strict=False) for s in states]
    # Reject all material unknown states before returning a density or selected mu.
    for state in filled:
        if state['ambiguous'].any():
            error=AmbiguousOccupation(state);error.number_diagnostics=diagnostic
            error.all_states=filled
            raise error
    electrons=float(sum(g*s['occupation'].sum() for g,s in zip(degeneracies,filled)))
    return mu,filled,{**diagnostic,'electrons':electrons,
        'number_residual_electrons':electrons-target,
        'ambiguous_tail_electron_uncertainty_bound':float(nhi-nlo),
        'number_prescription':'exact when interval collapses; otherwise tail-bounded representative, not exact neutrality'}
