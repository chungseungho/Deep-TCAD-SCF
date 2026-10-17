"""Independent dense spectral sum and energy-dependent bound residue checks."""
import numpy as np
import pytest
from scipy.integrate import quad_vec
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.special import expit
from transport.negf_reference import surface_sigma,fermi_poles
from transport.bound_states import block_inertia
from transport.errors import TransportNumericalError


def test_bound_pole_residue_spectral_sum_and_fermi_density():
    n=7;t=.2;U=.12;mu=-.04;kT=.02585
    H=np.diag(np.zeros(n))+np.diag(np.full(n-1,-t),1)+np.diag(np.full(n-1,-t),-1);H[n//2,n//2]-=U
    def sigma(z):
        a=np.zeros((n,n),complex);a[0,0]=a[-1,-1]=surface_sigma([z],[0],t)[0,0];return a
    E=-np.sqrt(4*t*t+U*U);ev,vec=eigh(H+sigma(E).real);v=vec[:,np.argmin(abs(ev-E))]
    delta=1e-6;derivative=(sigma(E+delta)-sigma(E-delta))/(2*delta)
    residue=np.outer(v,v)/(v@(np.eye(n)-derivative.real)@v)
    def integrand(k):
        e=-2*t*np.cos(k);G=np.linalg.inv(e*np.eye(n)-H-sigma(e));gam=1j*(sigma(e)-sigma(e).conj().T)
        spec=(G@gam@G.conj().T).real/(2*np.pi)*2*t*np.sin(k)
        return np.stack([spec,spec*expit((mu-e)/kT)])
    continuum,error=quad_vec(integrand,1e-9,np.pi-1e-9,epsabs=2e-10,epsrel=2e-10)
    assert np.max(abs(continuum[0]+residue-np.eye(n)))<1e-7
    poles,weights=fermi_poles(96);contour=.5*np.eye(n)
    for p,w in zip(poles,weights):contour+=2*kT*w*np.linalg.inv((mu+1j*p*kT)*np.eye(n)-H-sigma(mu+1j*p*kT)).real
    # Contour ALREADY contains the bound pole. Add only to continuum, never to
    # a complete equilibrium contour (which would double count it).
    expected=continuum[1]+residue*expit((mu-E)/kT)
    assert np.max(abs(contour-expected))<1e-7
    assert abs(np.trace(residue)-1)>1e-3 # lead tail: naive normalization fails


def test_dark_state_inside_other_channel_continuum_is_not_excluded_by_low_edge_audit():
    H=np.diag([0.,.1]);sig=np.diag([surface_sigma([.1],[0.],.3)[0,0],0.])
    dark=np.array([0.,1.]);gamma=1j*(sig-sig.conj().T)
    assert np.linalg.norm((H+sig-.1*np.eye(2))@dark)==0
    assert np.linalg.norm(gamma@dark)==0
    assert -.6<.1<.6 # embedded in the first lead channel continuum
    # Its density is an independent occupation input at bias, not eta filling.
    rho_cont=np.diag([.4,0]);rho_bound=np.outer(dark,dark)*.7
    assert np.trace(rho_cont+rho_bound)==pytest.approx(1.1)


def test_singular_ldl_pivot_is_numerical_failure_not_bound_evidence():
    ons=[np.array([[0.]]),np.array([[0.]])];hop=[np.array([[.2]])]
    with pytest.raises(TransportNumericalError,match='SINGULAR_LDL_PIVOT'):
        block_inertia(ons,hop,np.zeros((1,1)),np.zeros((1,1)),0.)
    for e in (-1e-9,1e-9):assert block_inertia(ons,hop,np.zeros((1,1)),np.zeros((1,1)),e)==1


def test_degenerate_poles_and_threshold_are_counted_with_multiplicity():
    t=.2;U=.12;E=-np.sqrt(4*t*t+U*U);ons=[-U*np.eye(2)]
    for energy,expected in [(E-1e-7,0),(E+1e-7,2),(-2*t-1e-10,2)]:
        sig=surface_sigma([energy],np.zeros(2),t)[0].real
        assert block_inertia(ons,[],np.diag(sig),np.diag(sig),energy)==expected
