"""Finite-electrode correlation from complex eigenstates and analytic Fermi integrals.

H_eff = H - i(Gamma_L+Gamma_R)/2 for an EXPLICIT wide-band external battery
on the finite L/R regions. This is a bath embedding self-energy in substance;
it is not a closed finite system and is not MS constrained-search DFT.

The integral int f_alpha(E) G Gamma_alpha G^dagger dE/(2pi) is evaluated
analytically using digamma functions in the eigenbasis of H_eff. SCF evaluates
neither a Green matrix at each energy nor lesser G or RGF. It remains
mathematically equivalent to this finite-bath correlation, not 'no embedding'.

Equilibrium is evaluated as a matrix function. A bias correction avoids an
ill-conditioned population solve for states weakly connected to reservoirs.
All states in the RETAINED mode space are included, not an occupied subset.
"""
import numpy as np
from scipy.linalg import eig,solve
from scipy.special import digamma
from config import KB_EV


def broadened_occupation(e,gamma,mu,temperature):
    a=.5+(gamma/2+1j*(np.asarray(e)-mu))/(2*np.pi*KB_EV*temperature)
    return .5-np.imag(digamma(a))/np.pi


def solve_finite_spectral(H,regions,gamma,mu_l,mu_r,temperature,split=False,diagnostic_detail=True):
    gl=np.diag((regions==0).astype(float)*gamma)
    gr=np.diag((regions==2).astype(float)*gamma);G=gl+gr
    z,V=eig(H-.5j*G,check_finite=False)
    inv=solve(V,np.eye(len(H)),assume_a='gen',check_finite=False)
    # Dissipativity identity provides nonnegative widths even for <1e-16 eV
    # resonances, whose imaginary eigenvalues cannot be obtained by subtraction.
    widths=np.real(np.sum(V.conj()*(G@V),axis=0)/np.sum(abs(V)**2,axis=0))
    z=z.real-.5j*widths
    scale=2*np.pi*KB_EV*temperature
    psi_r=digamma(.5+1j*(z-mu_r)/scale)
    psi_l=digamma(.5+1j*(z-mu_l)/scale)
    M=(V*psi_r)@inv
    equilibrium=.5*np.eye(len(H))-(M-M.conj().T)/(2j*np.pi)
    # For the real, time-reversal-symmetric EMA H and real diagonal Gamma,
    # equilibrium C is exactly real. Enforce this identity before adding the
    # bias correction, so roundoff in f(H_eff) cannot look like off-current.
    equilibrium=(equilibrium.real+equilibrium.real.T)/2
    # Positive factor form avoids cancellation in <left eigenvector|Gamma|left>.
    Fl=inv[:,regions==0]*np.sqrt(gamma)
    Fr=inv[:,regions==2]*np.sqrt(gamma)
    L=Fl@Fl.conj().T;R=Fr@Fr.conj().T
    denominator=2*np.pi*(z[:,None]-z.conj()[None,:])
    dp=psi_l-psi_r
    numerator=dp[:,None]-dp.conj()[None,:]
    kernel=np.divide(numerator,denominator,out=np.zeros_like(numerator),where=abs(denominator)>1e-290)
    # The diagonal has a stable scalar form and is sensitive to very narrow states.
    delta=broadened_occupation(z.real,widths,mu_l,temperature)-broadened_occupation(z.real,widths,mu_r,temperature)
    np.fill_diagonal(kernel,np.divide(delta,widths,out=np.zeros_like(delta),where=widths>1e-290))
    dark=widths<1e-28
    if np.any(dark&(abs(delta)>1e-7)):
        raise RuntimeError('Undetermined occupation of a disconnected bound state in the bias window; increase modes/electrode coupling')
    rho=equilibrium+V@(L*kernel)@V.conj().T
    rho=(rho+rho.conj().T)/2
    eig_rho=np.linalg.eigvalsh(rho)
    if eig_rho.min()<-2e-7 or eig_rho.max()>1+2e-7:
        raise RuntimeError(f'Non-positive analytic spectral density: {eig_rho.min():g}, {eig_rho.max():g}')
    components=None
    if split:
        num=psi_l[:,None]-psi_l.conj()[None,:]-1j*np.pi
        K=np.divide(num,denominator,out=np.zeros_like(num),where=abs(denominator)>1e-290)
        scalar=broadened_occupation(z.real,widths,mu_l,temperature)
        np.fill_diagonal(K,np.divide(scalar,widths,out=np.zeros_like(scalar),where=widths>1e-290))
        left=V@(L*K)@V.conj().T;left=(left+left.conj().T)/2
        components=(left,rho-left)
    diag={'pauli_min':float(eig_rho.min()),'pauli_max':float(eig_rho.max()),
          'minimum_effective_width_ev':float(widths.min()),'eigenvector_condition':float(np.linalg.cond(V)) if diagnostic_detail else None,
          'eigen_residual':float(np.linalg.norm((H-.5j*G)@V-V*z)/max(np.linalg.norm(H),1e-30)) if diagnostic_detail else None,
          'undetermined_dark_states':int(np.sum(dark&(abs(delta)>1e-7))),
          'density_algorithm':'analytic finite-bath Fermi correlation; no energy-resolved GF in SCF'}
    if not diagnostic_detail: diag['_correlation_spectrum']=eig_rho
    return rho,diag,components
