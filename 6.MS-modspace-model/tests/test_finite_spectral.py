import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad_vec,quad
from transport.finite_spectral import solve_finite_spectral,broadened_occupation
from multispace.occupation import fermi
from config import Q,HBAR_EV_S


def example():
    H=np.diag([.2,.3,.25])+np.diag([-.07,-.1],1)+np.diag([-.07,-.1],-1)
    regions=np.array([0,1,2]);gamma=.03
    return H,regions,gamma


def test_digamma_density_matches_independent_gf_quadrature():
    H,regions,gamma=example();GL=np.diag([gamma,0.,0]);GR=np.diag([0.,0,gamma])
    for ml,mr in [(.3,.2),(.3,.3)]:
        rho,diag,parts=solve_finite_spectral(H,regions,gamma,ml,mr,300,True)
        def integrand(E):
            # Independent NumPy solve; no reuse of eigendecomposition or kernel.
            G=np.linalg.solve(E*np.eye(3)-H+.5j*(GL+GR),np.eye(3))
            return G@(GL*fermi(E,ml,300)+GR*fermi(E,mr,300))@G.conj().T/(2*np.pi)
        direct,error=quad_vec(integrand,-np.inf,np.inf,epsabs=1e-10,epsrel=1e-10)
        assert_allclose(rho,direct,atol=2e-10)
        assert_allclose(parts[0]+parts[1],rho,atol=1e-13)
        assert diag['pauli_min']>=0 and diag['pauli_max']<=1


def test_bond_current_matches_independent_landauer_integral():
    H,regions,gamma=example();rho,_,_=solve_finite_spectral(H,regions,gamma,.3,.2,300)
    bond=-2*Q/HBAR_EV_S*np.imag(H[0,1]*rho[1,0])
    second=-2*Q/HBAR_EV_S*np.imag(H[1,2]*rho[2,1])
    def integrand(E):
        G=np.linalg.solve(E*np.eye(3)-H+.5j*np.diag([gamma,0,gamma]),np.eye(3))
        return gamma**2*abs(G[0,2])**2*(fermi(E,.3,300)-fermi(E,.2,300))*Q/HBAR_EV_S/(2*np.pi)
    current,_=quad(integrand,-2,2,epsabs=1e-17,epsrel=1e-10)
    assert_allclose(bond,current,rtol=1e-9,atol=1e-18)
    assert_allclose(bond,second,rtol=1e-12)


def test_zero_bias_reversal_and_energy_gauge():
    H,regions,gamma=example()
    rho,_,_=solve_finite_spectral(H,regions,gamma,.3,.2,300)
    zero,_,_=solve_finite_spectral(H,regions,gamma,.3,.3,300)
    reverse,_,_=solve_finite_spectral(H,regions,gamma,.2,.3,300)
    shifted,_,_=solve_finite_spectral(H+.8*np.eye(3),regions,gamma,1.1,1.,300)
    assert_allclose(zero.imag,0,atol=1e-15)
    assert_allclose(rho.imag,-reverse.imag,atol=1e-14)
    assert_allclose(rho,shifted,atol=1e-13)


def test_nearly_bound_barrier_states_do_not_break_pauli():
    H=np.diag([0.,.8,.81,.82,0.])+np.diag([-1e-4,-.01,-.01,-1e-4],1)+np.diag([-1e-4,-.01,-.01,-1e-4],-1)
    rho,diag,_=solve_finite_spectral(H,np.array([0,1,1,1,2]),.01,.03,-.07,300)
    assert diag['pauli_min']>-1e-8 and diag['pauli_max']<1+1e-8
    assert np.max(np.diag(rho)[1:4].real)<1e-6


def test_scalar_lorentzian_filling():
    for e in [-.2,.0,.2]:
        rho,_,_=solve_finite_spectral(np.array([[e]]),np.array([0]),.03,.01,.01,300)
        assert_allclose(rho[0,0].real,broadened_occupation(e,.03,.01,300),atol=1e-14)


def test_post_scf_pole_transmission_matches_direct_inverse():
    from transport.spectral import finite_transmission
    H,regions,gamma=example();E=np.linspace(-.1,.8,15)
    T=finite_transmission(H,regions,gamma,E)
    for e,t in zip(E,T):
        G=np.linalg.solve(e*np.eye(3)-H+.5j*np.diag([gamma,0,gamma]),np.eye(3))
        assert_allclose(t,gamma**2*abs(G[0,2])**2,rtol=1e-12,atol=1e-13)
