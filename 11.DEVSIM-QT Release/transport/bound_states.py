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


from transport.errors import TransportDomainError,TransportNumericalError

class BoundOccupationError(TransportDomainError):
    def __init__(self,diagnostics):
        self.bound_diagnostics=diagnostics
        super().__init__('BOUND_PREPARATION_REQUIRED: isolated poles have different reservoir Fermi occupations; no unique biased density',diagnostics)


def block_inertia(ons,hop,left,right,energy):
    """Negative inertia of H+Sigma(E)-E using a block LDL Sturm chain.

    Valid here only for real energies outside both lead continua. This changes
    the linear algebra used to locate poles, not their occupation or lifetime.
    """
    previous=None;negative=0
    for i,h in enumerate(ons):
        d=h.copy()-energy*np.eye(len(h))
        if i==0:d+=left
        if i==len(ons)-1:d+=right
        if previous is not None:
            scale=max(1.,np.linalg.norm(previous,2))
            if np.min(abs(np.linalg.eigvalsh(previous)))<1e-13*scale:
                # A singular block pivot is a numerical condition, not a pole
                # classification. Caller must bracket it from either side.
                raise TransportNumericalError('SINGULAR_LDL_PIVOT: move the audit energy off this pivot')
            d-=hop[i-1].conj().T@np.linalg.solve(previous,hop[i-1])
        d=(d+d.conj().T)/2
        negative+=int(np.count_nonzero(np.linalg.eigvalsh(d)<0.))
        previous=d
    return negative


def audit(system,c,mu,vd):
    ref=c.potential_reference_v
    ons,hop,tx,leads=device_blocks(system,c,(ref,ref+vd));off=np.r_[0,np.cumsum([len(h) for h in ons])]
    edge=min(float(e.min()-2*tx) for e,w in leads)
    upper=edge-max(1e-10,abs(edge)*1e-10)
    def count(e):
        sl,sr=lead_matrices(np.array([complex(e)]),tx,leads)
        if max(abs(sl.imag).max(),abs(sr.imag).max())>1e-10:raise ValueError('Expected evanescent self-energy')
        try:return block_inertia(ons,hop,sl[0].real,sr[0].real,e)
        except TransportNumericalError:
            delta=1e-11*max(1.,abs(e));counts=[]
            for shifted in (e-delta,e+delta):
                ll,rr=lead_matrices(np.array([complex(shifted)]),tx,leads)
                counts.append(block_inertia(ons,hop,ll[0].real,rr[0].real,shifted))
            if counts[0]!=counts[1]:
                # At a genuine pole use the lower-side count; enclosing root
                # uncertainty is bounded by delta and retained in diagnostics.
                return counts[0]
            return counts[0]
    number=count(upper);found=[]
    minima=[]
    for i,h in enumerate(ons):
        radius=np.sum(abs(h),axis=1)-abs(np.diag(h))
        if i:radius+=np.sum(abs(hop[i-1]),axis=0)
        if i<len(hop):radius+=np.sum(abs(hop[i]),axis=1)
        minima.append(float(min(np.diag(h)-radius)))
    lower=min(minima)-4*tx-1
    if count(lower)!=0:raise RuntimeError('Bound-pole lower bracket is not below the spectrum')
    for j in range(number):
        a,b=lower,upper
        for iteration in range(80):
            middle=(a+b)/2
            if count(middle)>j:b=middle
            else:a=middle
            if b-a<1e-12:break
        root=(a+b)/2
        fl=float(fermi(root,mu,c.temperature));fr=float(fermi(root,mu-vd,c.temperature))
        found.append({'energy_ev':root,'enclosure_ev':[a,b],'fermi_left':fl,'fermi_right':fr,
            'occupation_difference':abs(fl-fr),'needs_preparation_model':abs(fl-fr)>c.ms_occupation_uncertainty})
    return {'valley':system.valley.name,'lowest_propagating_edge_ev':edge,
        'bound_poles_below_both_leads':found,'bias_dependent_bound_occupation_unresolved':any(x['needs_preparation_model'] for x in found),
        'embedded_dark_states_excluded':False,'above_band_bound_states_excluded':False,
        'near_threshold_exclusion_ev':float(edge-upper),
        'classification_scope':'isolated poles of the retained open Hamiltonian; full-basis and contact convergence required to establish physical device states',
        'method':'block LDL inertia and monotone bisection of H+Sigma(E)-E below lowest lead band; no eta'}
