import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pytest
try:import gmsh,devsim
except (ImportError,OSError,RuntimeError) as exc:pytest.skip(str(exc),allow_module_level=True)
from validation.poisson_benchmarks import slab


def test_shared_dof_layered_dielectric_exact_capacitor(tmp_path):
    r=slab(tmp_path,.55)
    assert r['maximum_potential_error_v']<5e-13
    assert r['gauss_relative_error']<5e-13


def test_charged_dielectric_gauss_and_mesh_refinement(tmp_path):
    a=slab(tmp_path/'coarse',.7,11.7,11.7,1e18)
    b=slab(tmp_path/'fine',.3,11.7,11.7,1e18)
    assert a['gauss_relative_error']<5e-13 and b['gauss_relative_error']<5e-13
    assert b['weighted_l2_error_v']<a['weighted_l2_error_v']*.7
