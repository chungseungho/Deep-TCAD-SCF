import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pytest
from numpy.testing import assert_allclose
from config import Config
from multispace.projection import weights
from multispace.occupation import solve_states,classify,occupations,AmbiguousOccupation
from multispace.constraints import solve_number
from device.regions import five_regions
from device.definition import Frame


def test_projectors_positive_partition_and_number():
    rng=np.random.default_rng(4);u=np.linalg.qr(rng.normal(size=(9,9)))[0]
    w=weights(u,np.repeat([0,1,2],3))
    assert np.min(w)>=0;assert_allclose(w.sum(axis=1),1,atol=1e-14)
    f=np.arange(9)/9;rho=(u*f)@u.T
    assert_allclose(np.trace(rho),f.sum(),atol=1e-14)
    assert np.linalg.eigvalsh(rho).min()>-1e-14


def test_electrode_origins_no_mean_chemical_potential_and_gauge():
    c=Config();H=np.diag([0.,.9,.01]);labels=np.array([0,1,2])
    a=occupations(solve_states(H,labels,c),.1,-.1,300)
    b=occupations(solve_states(H+.73*np.eye(3),labels,c),.83,.63,300)
    assert_allclose(a['occupation'],b['occupation'],atol=1e-14)
    assert set(a['origin'])=={-1,0,2}
    assert a['occupation'][0]>.95 and a['occupation'][1]<.02


def test_bound_bias_window_is_rejected():
    c=Config();s=solve_states(np.diag([-.4,0.,.5]),np.array([0,1,2]),c)
    with pytest.raises(AmbiguousOccupation):occupations(s,.1,-.1,300)


def test_symmetric_delocalized_standing_state_is_rejected():
    c=Config();H=np.diag([.2,.2,.2])+np.diag([-.08,-.08],1)+np.diag([-.08,-.08],-1)
    s=solve_states(H,np.array([0,1,2]),c)
    with pytest.raises(AmbiguousOccupation):occupations(s,.25,.15,300)
    eq=occupations(s,.2,.2,300)
    rho=(s['u']*eq['occupation'])@s['u'].T
    assert np.max(abs(rho.imag))==0


def test_degenerate_subspace_rotation_never_manufactures_origin():
    c=Config();e=np.zeros(2);labels=np.array([0,2])
    for angle in [0,.13,.7,np.pi/4]:
        u=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        s=classify(e,u,labels,c)
        assert s['cluster_unstable'].all()
        with pytest.raises(AmbiguousOccupation):occupations(s,.1,-.1,300)
        eq=occupations(s,0.,0.,300)
        assert_allclose((u*eq['occupation'])@u.T,.5*np.eye(2),atol=1e-14)


def test_near_degenerate_opposite_origins_are_flagged():
    c=Config(ms_degenerate_ev=1e-6)
    s=classify(np.array([0.,1e-8]),np.eye(2),np.array([0,2]),c)
    assert s['cluster_unstable'].all()


def test_degenerate_all_left_subspace_is_safe_and_basis_invariant():
    c=Config();s=classify(np.zeros(2),np.eye(2),np.array([0,0]),c)
    eq=occupations(s,.1,-.1,300)
    assert (eq['origin']==0).all();assert not eq['ambiguous'].any()


def test_boundary_origin_flip_is_rejected():
    c=Config();u=np.eye(3);e=np.arange(3)*.1
    s=classify(e,u,np.array([0,1,2]),c,[np.array([2,1,0])])
    assert s['boundary_unstable'].all()


def test_finite_number_constraint_and_bias_split():
    c=Config();s=solve_states(np.diag([-.1,.1,.2]),np.array([0,0,2]),c)
    mu,filled,diag=solve_number([s],[4],5.,.05,c,0.)
    assert abs(4*filled[0]['occupation'].sum()-5)<1e-10
    assert abs(diag['bias_constraint_error_ev'])<1e-15


def test_five_physical_regions_and_axis_permutations():
    c=Config();x=np.array([1,13,20,42,53])
    assert_allclose(five_regions(x,c),[0,1,2,3,4])
    xyz=np.arange(15).reshape(5,3)
    for axis in (0,1,2):
        frame=Frame(axis);assert_allclose(frame.to_canonical(frame.to_world(xyz)),xyz)
        assert_allclose(frame.to_world(xyz)[:,axis],xyz[:,0])
