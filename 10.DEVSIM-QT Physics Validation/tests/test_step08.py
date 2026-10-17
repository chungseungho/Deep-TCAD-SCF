import numpy as np
import pytest
from config import Config
from multispace.constraints import solve_number,number_bounds
from multispace.occupation import classify,AmbiguousOccupation
from multispace.tracking import match_states
from coupling.domain_newton import solve,DomainSearchFailure


def test_number_interval_no_prior_assignment():
    c=Config();s=classify(np.array([0.]),np.array([[1.],[1.]])/np.sqrt(2),np.array([0,2]),c)
    with pytest.raises(AmbiguousOccupation) as e:solve_number([s],[2],1.,.1,c,0.)
    d=e.value.number_diagnostics
    assert d['chemical_potential_interval_ev']==pytest.approx([0.,.1],abs=1e-12)
    assert d['electron_number_interval_at_audit_mu'][0]<1<d['electron_number_interval_at_audit_mu'][1]
    assert not d['unresolved_assignment_in_number_equation']


def test_number_zero_bias_and_reservoir_swap():
    c=Config();s=classify(np.array([0.]),np.array([[1.],[1.]])/np.sqrt(2),np.array([0,2]),c)
    mu,filled,d=solve_number([s],[2],1.,0.,c,0.)
    assert mu==pytest.approx(0,abs=1e-12);assert d['electron_number_uncertainty_bound']==0
    assert number_bounds([s],[2],.03,.1,300)==pytest.approx(number_bounds([s],[2],-.07,-.1,300))


def test_tracking_crossing_and_complex_phase():
    a={'e':np.array([0.,1.,2.]),'u':np.eye(3,dtype=complex)}
    b={'e':np.array([-.1,.2,2.]),'u':np.eye(3,dtype=complex)[:,[1,0,2]]*1j}
    m=match_states(a,b);assert m['new_indices'].tolist()==[1,0,2]
    assert np.allclose(m['overlap_squared'],1)
    assert not m['defines_injection_origin']


def test_tracking_degenerate_rotation_invariant():
    a={'e':np.array([0.,0.,2.]),'u':np.eye(3)}
    q=np.array([[1,1,0],[-1,1,0],[0,0,np.sqrt(2)]])/np.sqrt(2)
    m=match_states(a,{'e':a['e'],'u':q})
    assert m['degenerate_subspaces'][0]['normalized_projector_overlap']==pytest.approx(1)
    assert np.allclose(m['degenerate_subspaces'][0]['singular_values'],1)


class InvalidTrial(Exception):pass


def test_newton_invalid_trial_recovers_physical_root():
    def f(x):
        if x[0]<=0:raise InvalidTrial('positive log domain')
        return np.log(x)
    x,info=solve(f,[3.],invalid=(InvalidTrial,),trust_radius=10,tol=1e-9)
    assert x[0]==pytest.approx(1,abs=1e-8)
    assert any(not e['valid'] for e in info['events'])


def test_newton_no_fake_root_on_domain_boundary():
    def f(x):
        if x[0]>=1:raise InvalidTrial('outside physical domain')
        return x-2
    with pytest.raises(DomainSearchFailure):solve(f,[0.],invalid=(InvalidTrial,),trust_radius=3,maxiter=80,max_backtracks=10)


def test_programming_errors_not_backtracked():
    def bad(x):raise KeyError('implementation fault')
    with pytest.raises(KeyError):solve(bad,[0.],invalid=(InvalidTrial,))


def test_contour_and_transmission_batching_preserve_values():
    from dataclasses import replace
    from validation.analytical import chain
    from transport.negf_reference import solve_negf
    from transport.ms_postprocess import transmission
    c,s,edge=chain(n=13,barrier=.07);mu=edge+.13
    a=solve_negf(s,replace(c,contour_batch_size=1),mu,mu,0)[0]
    b=solve_negf(s,replace(c,contour_batch_size=96),mu,mu,0)[0]
    assert np.allclose(a,b,atol=2e-14,rtol=0)
    e=np.linspace(edge+.001,edge+.5,101)
    assert np.allclose(transmission(s,replace(c,transmission_batch_size=1),e,0),
                       transmission(s,replace(c,transmission_batch_size=200),e,0),atol=2e-14,rtol=0)
