"""Real-axis bound-pole audit below every propagating lead band.

This diagnoses preparation-dependent populations, rather than filling them
from an arbitrary average chemical potential. Embedded dark states require
an additional symmetry/rank analysis and are not excluded by this audit.
"""
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from transport.negf_reference import device_blocks,lead_matrices
from multispace.occupation import fermi


class BoundOccupationError(RuntimeError):
    def __init__(self,diagnostics):
        self.bound_diagnostics=diagnostics
        super().__init__('BOUND_PREPARATION_REQUIRED: isolated poles have different reservoir Fermi occupations; no unique biased density')


def audit(system,c,mu,vd):
    ref=c.potential_reference_v
    ons,hop,tx,leads=device_blocks(system,c,(ref,ref+vd));off=np.r_[0,np.cumsum([len(h) for h in ons])]
    H=np.zeros((off[-1],off[-1]))
    for i,h in enumerate(ons):H[off[i]:off[i+1],off[i]:off[i+1]]=h
    for i,t in enumerate(hop):
        H[off[i]:off[i+1],off[i+1]:off[i+2]]=t;H[off[i+1]:off[i+2],off[i]:off[i+1]]=t.T
    edge=min(float(e.min()-2*tx) for e,w in leads)
    upper=edge-max(1e-10,abs(edge)*1e-10)
    def effective(e):
        sl,sr=lead_matrices(np.array([complex(e)]),tx,leads);h=H.copy()
        if max(abs(sl.imag).max(),abs(sr.imag).max())>1e-10:raise ValueError('Expected evanescent self-energy')
        h[:len(ons[0]),:len(ons[0])]+=sl[0].real
        h[-len(ons[-1]):,-len(ons[-1]):]+=sr[0].real
        return h
    eigen=eigh(effective(upper),eigvals_only=True,check_finite=False)
    count=int(sum(eigen<upper));found=[]
    lower=float(min(np.diag(H)-np.sum(abs(H),axis=1)+abs(np.diag(H)))-4*tx-1)
    for j in range(count):
        def equation(e):return e-float(eigh(effective(e),eigvals_only=True,subset_by_index=(j,j),check_finite=False)[0])
        root=brentq(equation,lower,upper,xtol=1e-12)
        fl=float(fermi(root,mu,c.temperature));fr=float(fermi(root,mu-vd,c.temperature))
        found.append({'energy_ev':root,'fermi_left':fl,'fermi_right':fr,
            'occupation_difference':abs(fl-fr),'needs_preparation_model':abs(fl-fr)>c.ms_occupation_uncertainty})
    return {'valley':system.valley.name,'lowest_propagating_edge_ev':edge,
        'bound_poles_below_both_leads':found,'bias_dependent_bound_occupation_unresolved':any(x['needs_preparation_model'] for x in found),
        'embedded_dark_states_excluded':False,'method':'nonlinear E-lambda_j[H+Sigma(E)] roots below lowest lead band'}
