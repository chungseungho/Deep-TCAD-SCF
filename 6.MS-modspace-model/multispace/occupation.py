"""MS-inspired spectral diagnostics. This is NOT a constrained-search DFT.

The supplied PNAS paper assigns channel states using L/R spatial weights.
An ambiguous finite standing wave has no unique origin: report it, do not
silently label it L, split it into artificial travelling states, or average mu.
"""
import numpy as np
from scipy.linalg import eigh
from scipy.special import expit
from config import KB_EV


def fermi(e,mu,temperature): return expit(-(np.asarray(e)-mu)/(KB_EV*temperature))


def partition_eigenstates(H,regions,mu_l,mu_r,temperature,tie_tol=1e-7):
    e,u=eigh(H,check_finite=False)
    # Only machine-degenerate subspaces may be rotated without changing H.
    # Localise them by diagonalising P_L-P_R, invariant to input eigenbasis.
    a=0
    label=(regions==0).astype(float)-(regions==2).astype(float)
    while a<len(e):
        b=a+1
        while b<len(e) and abs(e[b]-e[a])<1e-11: b+=1
        if b-a>1:
            sub=u[:,a:b];_,rot=eigh(sub.T@(label[:,None]*sub));u[:,a:b]=sub@rot
        a=b
    w=np.stack([(abs(u[regions==r])**2).sum(axis=0) for r in range(3)],axis=1)
    origin=np.where(w[:,0]>w[:,2],0,2)
    # States with largest electrode weight belong to that electrode.
    origin=np.where(np.argmax(w,axis=1)==0,0,origin)
    origin=np.where(np.argmax(w,axis=1)==2,2,origin)
    tie=abs(w[:,0]-w[:,2])<tie_tol
    bound=(w[:,0]+w[:,2])<1e-8
    active=(e>min(mu_l,mu_r)-8*KB_EV*temperature)&(e<max(mu_l,mu_r)+8*KB_EV*temperature)
    ambiguous=(tie|bound)&active&(abs(mu_l-mu_r)>1e-12)
    occupations=np.where(origin==0,fermi(e,mu_l,temperature),fermi(e,mu_r,temperature))
    if abs(mu_l-mu_r)<1e-12: occupations=fermi(e,mu_l,temperature)
    return {'e':e,'u':u,'weights':w,'origin':origin,'occupation':occupations,
            'ambiguous':ambiguous,'bound':bound,'active':active}


def closed_ms_density(H,regions,mu_l,mu_r,temperature,tie_tol=1e-7):
    state=partition_eigenstates(H,regions,mu_l,mu_r,temperature,tie_tol)
    if np.any(state['ambiguous']):
        raise RuntimeError(f"MS origin unresolved for {int(state['ambiguous'].sum())} active finite eigenstates; "
                           'use finite driven reservoirs or increase electrode size; no guessed occupation')
    rho=(state['u']*state['occupation'])@state['u'].T
    return rho,state
