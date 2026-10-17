import numpy as np
from scipy.integrate import quad_vec
from transport.resonances import resonance_breaks
from transport.negf_reference import rgf,lead_matrices
from transport.quadrature import integrate_vector


def test_banded_guidance_finds_narrow_resonance_against_dense_reference():
    t=.2;weak=.0007;ons=[np.array([[.4]]) for _ in range(3)];hop=[np.array([[-weak]])]*2;leads=[(np.array([.4]),np.eye(1))]*2
    breaks,diag=resonance_breaks(ons,hop,t,leads,.2,.6,.02)
    def production(E):
        sl,sr=lead_matrices(E,t,leads);_,A=rgf(E,ons,hop,sl,sr,False,0)
        return (abs(A[1][:,0,0])**2/(2*np.pi))[:,None]
    integrated,error=integrate_vector(production,breaks,8,2e-5,1e-9,120000)
    H=np.diag([.4]*3)+np.diag([-weak]*2,1)+np.diag([-weak]*2,-1)
    def reference(E):
        sl,sr=lead_matrices(np.array([E]),t,leads);sigma=np.diag([sl[0,0,0],0,sr[0,0,0]])
        G=np.linalg.inv(E*np.eye(3)-H-sigma);gamma=-2*sl[0,0,0].imag
        return abs(G[1,0])**2*gamma/(2*np.pi)
    width=2*weak**2/t
    points=np.unique(np.r_[.2,.6,.4,.4+width*np.array([-1000,-100,-10,-1,1,10,100,1000])])
    exact=sum(quad_vec(reference,a,b,epsabs=1e-10,epsrel=1e-9)[0] for a,b in zip(points[:-1],points[1:]))
    assert abs(integrated[0]-exact)<2e-6
    assert abs(exact-.5)<2e-4
    assert 'band' in diag['storage']
