import numpy as np
from transport.bound_states import block_inertia,audit
from validation.reference_bound_audit import audit as dense_audit
from validation.analytical import chain


def test_block_sturm_inertia_matches_independent_dense_eigenvalues():
    rng=np.random.default_rng(31);sizes=[2,3,1,2,3]
    ons=[];hop=[]
    for m in sizes:
        a=rng.normal(size=(m,m));ons.append((a+a.T)/2)
    for a,b in zip(sizes[:-1],sizes[1:]):hop.append(rng.normal(size=(a,b))*.2)
    off=np.r_[0,np.cumsum(sizes)];H=np.zeros((off[-1],off[-1]))
    for i,h in enumerate(ons):H[off[i]:off[i+1],off[i]:off[i+1]]=h
    for i,t in enumerate(hop):
        H[off[i]:off[i+1],off[i+1]:off[i+2]]=t;H[off[i+1]:off[i+2],off[i]:off[i+1]]=t.T
    for e in np.linspace(-4.,4.,31):
        expected=int(sum(np.linalg.eigvalsh(H)<e))
        actual=block_inertia(ons,hop,np.zeros_like(ons[0]),np.zeros_like(ons[-1]),e)
        assert actual==expected


def test_impurity_bound_pole_matches_analytic_infinite_chain_and_dense_audit():
    c,s,edge=chain(n=15);strength=.12;s.H[7,7]-=strength
    result=audit(s,c,edge,.0);reference=dense_audit(s,c,edge,.0)
    roots=result['bound_poles_below_both_leads'];expected=edge+.4-np.sqrt(.4**2+strength**2)
    assert len(roots)==1
    assert abs(roots[0]['energy_ev']-expected)<2e-11
    assert abs(roots[0]['energy_ev']-reference['bound_poles_below_both_leads'][0]['energy_ev'])<2e-11


def test_positive_barrier_has_no_below_band_poles():
    c,s,edge=chain(n=15,barrier=.12)
    assert audit(s,c,edge,0.)['bound_poles_below_both_leads']==[]
