import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
import pytest
try:
    import gmsh,devsim
except (ImportError,OSError,RuntimeError) as exc:
    pytest.skip(str(exc),allow_module_level=True)
from config import Config
from geometry.build import build_mesh
from poisson.device import Poisson


def test_actual_three_material_gaa_shared_potential_and_gauss(tmp_path):
    c=Config(ls=2,ld=2,lg=4,spacer=1,width=2,thickness=2,sd_width=6,sd_thickness=6,tox=1,
             dx=.5,dy=1,dz=1,mesh_h=1.8,channel_h=.8,junction_h=.7)
    mesh=build_mesh(c,tmp_path);p=Poisson(c,mesh)
    try:
        p.set_charge(p.donors);zero=p.solve(c.flatband,0.)
        assert_allclose(zero,0,atol=1e-13)
        p.set_charge(p.donors*.25);p.solve(.9,.1)
        balance=p.gauss_balance()
        assert balance['gauss_relative_error']<5e-13
        assert balance['contact_sum_relative_error']<5e-13
        assert balance['maximum_interior_equation_residual_c']<1e-30
        for r0,r1 in [('silicon','oxide'),('silicon','spacer'),('oxide','spacer')]:
            common,a,b=np.intersect1d(p.indices[r0],p.indices[r1],return_indices=True)
            assert len(common)>0
            assert_allclose(p.values('Potential',r0)[a],p.values('Potential',r1)[b],atol=0)
    finally: p.close()
