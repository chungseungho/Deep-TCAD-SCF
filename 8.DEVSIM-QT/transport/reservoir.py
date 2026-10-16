"""Finite extended-reservoir Lindblad transport, not the MS-DFT algorithm.

A=-iH-Gamma/2; A C + C A^dagger = -Pi, in eV time units.
Each isolated L/R eigenmode is driven at gamma*f(e_k-mu_alpha). This is a
specified Markovian external battery, NOT a closed finite system. No retarded
or lesser Green function, energy integration, or semi-infinite lead is used.

Refs: Elenewski/Gruss/Zwolak JCP 147,151101 (2017),
Gruss/Smolyanitsky/Zwolak JCP 147,141102 (2017).
Reservoir spacing << gamma << kBT and length/gamma convergence remain essential.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh, solve_continuous_lyapunov
from scipy.optimize import brentq
from multispace.occupation import fermi
from config import Q, HBAR_EV_S, KB_EV


@dataclass
class Reservoirs:
    gamma: np.ndarray
    injection_l: np.ndarray
    injection_r: np.ndarray
    levels_l: np.ndarray
    levels_r: np.ndarray


def reservoirs(H,regions,gamma,mu_l,mu_r,temperature):
    G=np.zeros_like(H);ps=[];levels=[]
    for r,mu in ((0,mu_l),(2,mu_r)):
        ind=np.where(regions==r)[0]
        e,u=eigh(H[np.ix_(ind,ind)],check_finite=False)
        p=np.zeros_like(H);p[np.ix_(ind,ind)]=gamma*(u*fermi(e,mu,temperature))@u.T
        G[ind,ind]=gamma;ps.append(p);levels.append(e)
    return Reservoirs(G,*ps,*levels)


def solve_steady(H,baths,split=False):
    A=-1j*H-.5*baths.gamma
    Pi=baths.injection_l+baths.injection_r
    rho=solve_continuous_lyapunov(A,-Pi)
    rho=(rho+rho.conj().T)/2
    residual=np.linalg.norm(A@rho+rho@A.conj().T+Pi)/max(np.linalg.norm(Pi),1e-30)
    eigen=np.linalg.eigvalsh(rho)
    if residual>1e-7 or eigen.min()<-1e-7 or eigen.max()>1+1e-7:
        raise RuntimeError(f'Invalid Lindblad solution: residual={residual:g}, Pauli=({eigen.min():g},{eigen.max():g})')
    components=None
    if split:
        l=solve_continuous_lyapunov(A,-baths.injection_l)
        components=(.5*(l+l.conj().T),rho-.5*(l+l.conj().T))
    return rho,{'lyapunov_residual':float(residual),'pauli_min':float(eigen.min()),
                'pauli_max':float(eigen.max())},components


def bond_currents(sys,rho,degeneracy):
    # Positive means electron flux from left to right, customary positive Id.
    return np.array([-degeneracy*2*Q/HBAR_EV_S*np.imag(np.trace(
        h@rho[sys.sl(i+1),sys.sl(i)])) for i,h in enumerate(sys.hop)])


def contact_currents(sys,rho,baths,degeneracy):
    out=[]
    for r,p in ((0,baths.injection_l),(2,baths.injection_r)):
        ind=np.where(sys.regions==r)[0]
        rate=np.trace(p[np.ix_(ind,ind)]-baths.gamma[np.ix_(ind,ind)]@rho[np.ix_(ind,ind)]).real
        out.append(degeneracy*Q/HBAR_EV_S*rate)
    return out


def reservoir_spacing(levels,mu,temperature):
    e=np.sort(levels);e=e[(e>mu-4*KB_EV*temperature)&(e<mu+4*KB_EV*temperature)]
    # Exact degeneracies are separate channels, NOT a dense transport continuum.
    d=np.diff(np.unique(np.round(e,10)))
    return float(np.max(d)) if len(d) else float('inf')


def neutral_mu(builders,c,grid):
    """One common mu from all valleys and isolated finite source spectrum.

    Match nominal donor number in quantum cells, with the SAME modes/masses/
    grid/finite length as transport. Do not adjust doping or mass to match I-V.
    """
    if c.backend=='negf':
        from transport.negf_reference import neutral_lead_mu
        return neutral_lead_mu(builders,c,grid)
    spectra=[]
    for b in builders:
        sys=b.build(np.zeros(len(grid.xyz)));ind=np.where(sys.regions==0)[0]
        e=eigh(sys.H[np.ix_(ind,ind)],eigvals_only=True,check_finite=False)
        spectra.append((e,2*b.v.multiplicity))
    target=c.donor_sd*c.ls*c.sd_width*c.sd_thickness*1e-21
    capacity=sum(g*len(e) for e,g in spectra)
    if target>=capacity*.95: raise ValueError('Insufficient source modes for neutral doping')
    if c.backend=='spectral':
        from transport.finite_spectral import broadened_occupation
        filling=lambda e,mu:broadened_occupation(e,c.gamma,mu,c.temperature)
    else: filling=lambda e,mu:fermi(e,mu,c.temperature)
    fun=lambda mu: sum(g*filling(e,mu).sum() for e,g in spectra)-target
    low=min(e.min() for e,_ in spectra)-2
    high=max(e.max() for e,_ in spectra)+2
    mu=brentq(fun,low,high,xtol=1e-13)
    physical_mu=float(mu-c.potential_reference_v)
    return physical_mu,{'mu_source_ev':physical_mu,'zero_reference_mu_ev':float(mu),'potential_reference_v':c.potential_reference_v,'donor_electrons_source':float(target),
                     'isolated_source_electrons':float(fun(mu)+target),
                     'neutrality_relative_error':float(abs(fun(mu))/target),
                     'calibration':'isolated finite source, all retained valleys and configured bath filling; coupled contact neutrality checked separately'}
