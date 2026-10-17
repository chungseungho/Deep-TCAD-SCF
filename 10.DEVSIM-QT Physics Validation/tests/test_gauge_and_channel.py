import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import replace
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.linalg import null_space
from config import Config
from multispace.occupation import solve_states,occupations,AmbiguousOccupation,fermi
from validation.analytical import chain
from transport.negf_reference import solve_negf
from transport.ms_postprocess import transmission


def test_c_dominant_supported_left_state_uses_left_fermi():
    vector=np.sqrt([.18,.80,.02]);u=np.column_stack([vector,null_space(vector[None,:])])
    H=u@np.diag([0.,.7,1.])@u.T
    s=occupations(solve_states(H,np.array([0,1,2]),Config()),.05,-.05,300)
    assert s['largest_region'][0]==1 and s['origin'][0]==0
    assert_allclose(s['occupation'][0],fermi(0.,.05,300),atol=1e-14)


def test_weak_nonzero_reservoir_reach_is_not_an_origin():
    vector=np.sqrt([8e-10,1-1e-9,2e-10]);u=np.column_stack([vector,null_space(vector[None,:])])
    H=u@np.diag([0.,.7,1.])@u.T
    s=solve_states(H,np.array([0,1,2]),Config())
    assert s['reach'][0]>0
    with pytest.raises(AmbiguousOccupation):occupations(s,.05,-.05,300)


def test_open_density_and_transmission_common_potential_gauge():
    c,system,edge=chain(n=7);mu=edge+.2;vd=.03
    blocks,flux,_,_=solve_negf(system,c,mu,mu-vd,vd)
    offset=.4;shifted=replace(system,H=system.H-offset*np.eye(system.dimension),
        subbands=[e-offset for e in system.subbands])
    cc=replace(c,potential_reference_v=offset)
    blocks2,flux2,_,_=solve_negf(shifted,cc,mu-offset,mu-vd-offset,vd)
    assert_allclose(np.asarray(blocks),np.asarray(blocks2),atol=2e-10)
    assert_allclose(flux,flux2,rtol=1e-8,atol=1e-17)
    E=np.array([edge+.1,edge+.2,edge+.4])
    assert_allclose(transmission(system,c,E,vd),transmission(shifted,cc,E-offset,vd),atol=1e-12)
