"""Regression checks for code identity and restored electrostatic exports."""
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
import pytest
import provenance


def test_source_identity_ignores_environment_but_detects_simulator_changes(tmp_path,monkeypatch):
    monkeypatch.setattr(provenance,'__file__',str(tmp_path/'provenance.py'))
    core=tmp_path/'solver.py';core.write_text('physical_parameter = 1\n')
    baseline=provenance.source_digest()
    for name in ['.venv/lib/site.py','venv/lib/site.py','results/generated.py','tests/example.py','validation/example.py']:
        file=tmp_path/name;file.parent.mkdir(parents=True,exist_ok=True);file.write_text('dependency = 2\n')
    assert provenance.source_digest()==baseline
    core.write_text('physical_parameter = 3\n')
    assert provenance.source_digest()!=baseline


def test_restored_export_recovers_all_material_fields_and_reports_units(tmp_path):
    try:import gmsh,devsim
    except (ImportError,OSError,RuntimeError) as exc:pytest.skip(str(exc))
    from config import Config
    from geometry.build import build_mesh
    from poisson.device import Poisson
    from export_fields import restore_last_poisson
    c=Config(ls=2,ld=2,lg=4,spacer=1,width=2,thickness=2,sd_width=6,sd_thickness=6,tox=1,
             dx=.5,dy=1,dz=1,mesh_h=1.8,channel_h=.8,junction_h=.7)
    p=Poisson(c,build_mesh(c,tmp_path/'mesh'))
    try:
        density=p.donors*.45;p.set_charge(density);saved=p.solve(.9,.1);all_saved=p._get('Potential').copy()
        states=tmp_path/'states';states.mkdir()
        np.savez_compressed(states/'vg0.900000_vd0.100000.npz',converged=True,
            node_xyz_nm=p.xyz['silicon'],poisson_density_cm3=density,node_phi_v=saved,
            vg_v=.9,vd_v=.1,source_sha256='synthetic-electrostatic-test',configuration_fingerprint='fixture')
        p.set_charge(p.donors);p.solve(c.flatband,0.)
        assert np.max(abs(p._get('Potential')-all_saved))>.1
        report=restore_last_poisson(p,tmp_path)
        assert_allclose(p._get('Potential'),all_saved,atol=2e-13)
        assert report['maximum_saved_node_phi_error_v']<2e-13
        metadata=json.loads((tmp_path/'field_export.json').read_text())
        assert metadata['native_devsim_coordinates']=='cm'
        assert metadata['quantum_vtu_coordinates']=='nm'
    finally:p.close()
