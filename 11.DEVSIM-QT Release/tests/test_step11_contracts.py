import numpy as np
import pytest
from dataclasses import replace
from validation.analytical import chain
from transport.backends.negf import NEGFBackend
from config import Config,VALLEYS
from modespace.hamiltonian import Grid,ModeBuilder


def test_uncomputed_negf_contact_density_is_none_and_not_zero_export():
    c,s,e=chain(n=5);c=replace(c,ls=1,ld=1,lg=1,spacer=1,reservoir_left_nm=1,reservoir_right_nm=1)
    class Builder:
        def build(self,potential):return s
    backend=NEGFBackend(c);r=backend.solve(np.zeros(5),0.,{'builders':[Builder()],'grid':s.grid,'fixed_mu':e+.1,'bias':0.,'detail':True})
    assert r['density_l'] is None and r['density_r'] is None and r['density_unresolved'] is None
    exported=backend.export_fields(r)
    assert str(exported['source_drain_density_status'])=='UNSUPPORTED'
    assert 'density_from_l_cm3' not in exported


def test_mutated_config_cannot_reuse_stale_basis_cache():
    c=Config();b=ModeBuilder(c,Grid.build(c),VALLEYS[0]);c.modes+=1
    with pytest.raises(ValueError,match='Configuration changed'):b.build(np.zeros(len(b.g.xyz)))


def test_scf_delegates_model_scope_assessment_to_backend():
    from pathlib import Path
    text=(Path(__file__).resolve().parents[1]/'solver.py').read_text()
    assert 'backend.physical_validity(' in text
    assert 'validity=evaluate(' not in text
