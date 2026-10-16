"""Published spatial-origin rule with rejection; no state rotations or mean μ."""
import numpy as np
from scipy.linalg import eigh
from scipy.special import expit
from config import KB_EV
from multispace.projection import weights,degenerate_groups,origin_spectrum


def fermi(e,mu,temperature): return expit(-(np.asarray(e)-mu)/(KB_EV*temperature))


class AmbiguousOccupation(RuntimeError):
    def __init__(self,state):
        self.state=state
        ids=np.flatnonzero(state['ambiguous'])
        super().__init__(f'MS_OCCUPATION_UNRESOLVED: {len(ids)} active states; indices {ids[:12].tolist()}; no guessed filling')


def classify(e,u,regions,c,alternative_regions=()):
    w=weights(u,regions);reach=w[:,0]+w[:,2]
    contrast=(w[:,0]-w[:,2])/np.maximum(reach,1e-300)
    origin=np.where(contrast>0,0,2)
    uncertain=(reach<c.ms_reach_min)|(abs(contrast)<c.ms_contrast_min)
    boundary_unstable=np.zeros(len(e),bool)
    for labels in alternative_regions:
        ww=weights(u,labels);delta=ww[:,0]-ww[:,2]
        boundary_unstable|=(np.sign(delta)!=np.sign(contrast))|(abs(delta)<c.ms_tie_tol)
    uncertain|=boundary_unstable
    cluster_unstable=np.zeros(len(e),bool);clusters=[]
    for sl in degenerate_groups(e,c.ms_degenerate_ev):
        if sl.stop-sl.start<2:continue
        spectrum=origin_spectrum(u[:,sl],regions)
        supported=(spectrum.min()>c.ms_tie_tol) or (spectrum.max()<-c.ms_tie_tol)
        if not supported:cluster_unstable[sl]=True
        clusters.append({'first':sl.start,'last_exclusive':sl.stop,'projector_min':float(spectrum.min()),
                         'projector_max':float(spectrum.max()),'origin_invariant':bool(supported)})
    uncertain|=cluster_unstable
    origin[uncertain]=-1
    return {'e':e,'u':u,'weights':w,'largest_region':np.argmax(w,axis=1),'origin':origin,
        'uncertain_origin':uncertain,'boundary_unstable':boundary_unstable,
        'cluster_unstable':cluster_unstable,'reach':reach,'contrast':contrast,
        'bound':reach<c.ms_reach_min,'clusters':clusters}


def occupations(state,mu_l,mu_r,temperature,tolerance=1e-7,strict=True):
    fl=fermi(state['e'],mu_l,temperature);fr=fermi(state['e'],mu_r,temperature)
    ambiguous=state['uncertain_origin']&(abs(fl-fr)>tolerance)
    # For equal/tail occupations the difference is explicitly bounded. Never mean μ.
    f=np.where(state['origin']==2,fr,fl)
    out={**state,'occupation':f,'ambiguous':ambiguous,'active':abs(fl-fr)>tolerance,
         'occupation_uncertainty_per_state':np.where(state['uncertain_origin'],abs(fl-fr),0.),
         'mu_l_ev':mu_l,'mu_r_ev':mu_r}
    if strict and ambiguous.any():raise AmbiguousOccupation(out)
    return out


def solve_states(H,regions,c,alternative_regions=()):
    if np.max(abs(H-H.conj().T))>1e-10:raise ValueError('Hamiltonian is not Hermitian')
    e,u=eigh(H,check_finite=False)
    return classify(e,u,regions,c,alternative_regions)


def partition_eigenstates(H,regions,mu_l,mu_r,temperature,tie_tol=1e-7):
    # Compatibility API for independent inherited finite-contact benchmarks.
    from config import Config
    c=Config(ms_tie_tol=tie_tol)
    return occupations(solve_states(H,regions,c),mu_l,mu_r,temperature,strict=False)


def closed_ms_density(H,regions,mu_l,mu_r,temperature,tie_tol=1e-7):
    s=partition_eigenstates(H,regions,mu_l,mu_r,temperature,tie_tol)
    if s['ambiguous'].any():raise AmbiguousOccupation(s)
    return (s['u']*s['occupation'])@s['u'].conj().T,s
