"""Finite electron-number constraint, not a variational DFT functional."""
import numpy as np
from scipy.optimize import brentq
from multispace.occupation import fermi,occupations


def solve_number(states,degeneracies,target,vd,c,fixed_mu):
    def count(mu):
        total=0.
        for s,g in zip(states,degeneracies):
            f=np.where(s['origin']==2,fermi(s['e'],mu-vd,c.temperature),fermi(s['e'],mu,c.temperature))
            total+=g*f.sum()
        return total
    if c.ms_number_constraint=='finite_total':
        if not 0<target<sum(len(s['e'])*g for s,g in zip(states,degeneracies)):
            raise ValueError('Not enough states for electron count')
        lower=min(s['e'].min() for s in states)-abs(vd)-5
        upper=max(s['e'].max() for s in states)+abs(vd)+5
        mu=brentq(lambda m:count(m)-target,lower,upper,xtol=1e-13)
    else:mu=float(fixed_mu)
    filled=[occupations(s,mu,mu-vd,c.temperature,c.ms_occupation_uncertainty) for s in states]
    uncertainty=sum(g*s['occupation_uncertainty_per_state'].sum() for g,s in zip(degeneracies,filled))
    return mu,filled,{'constraint':c.ms_number_constraint,'target_electrons':float(target),
        'electrons':float(count(mu)),'number_residual_electrons':float(count(mu)-target),
        'ambiguous_tail_electron_uncertainty_bound':float(uncertainty),
        'bias_constraint_error_ev':float(mu-(mu-vd)-vd),'variational_functional_implemented':False}
