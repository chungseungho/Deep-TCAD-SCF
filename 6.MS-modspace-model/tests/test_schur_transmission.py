import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from types import SimpleNamespace
from numpy.testing import assert_allclose
from validation.post_scf_current import finite_rgf_transmission


def test_distributed_finite_bath_schur_transmission_matches_direct_full_inverse():
    rng=np.random.default_rng(3);n=11;t=.2
    H=np.diag(.4+rng.random(n)*.1)+np.diag(np.full(n-1,-t),1)+np.diag(np.full(n-1,-t),-1)
    regions=np.r_[np.zeros(3,int),np.ones(5,int),np.full(3,2)]
    system=SimpleNamespace(H=H,grid=SimpleNamespace(x=np.arange(n)),hop=[np.array([[-t]])]*(n-1),
                           offsets=np.arange(n+1),regions=regions,sl=lambda i:slice(i,i+1))
    energy=np.linspace(-.05,.95,111);gamma=.013
    T=finite_rgf_transmission(system,gamma,energy)
    for e,value in zip(energy,T):
        G=np.linalg.inv(e*np.eye(n)-H+.5j*np.diag((regions!=1)*gamma))
        direct=gamma**2*np.sum(abs(G[np.ix_(regions==0,regions==2)])**2)
        assert_allclose(value,direct,rtol=5e-13,atol=1e-14)
