import json
from pathlib import Path
import numpy as np
import pytest
from config import Config
from provenance import validate_restart,source_digest


def test_restart_rejects_source_config_mesh_and_negative_charge(tmp_path):
    c=Config();xyz=np.zeros((2,3));p=tmp_path/'restart.npz'
    base=dict(config_json=json.dumps(c.to_dict()),source_sha256=source_digest(),node_xyz_nm=xyz,poisson_density_cm3=np.ones(2))
    np.savez(p,**base);assert validate_restart(p,c,xyz)['poisson_density_cm3'].shape==(2,)
    for key,val,match in [('source_sha256','wrong','source'),('config_json',json.dumps({**c.to_dict(),'flatband':.7}),'config'),('node_xyz_nm',xyz+1,'coordinate'),('poisson_density_cm3',np.array([-1,2]),'density')]:
        np.savez(p,**{**base,key:val})
        with pytest.raises(ValueError,match=match):validate_restart(p,c,xyz)


def test_common_scf_has_only_neutral_transport_domain_exception():
    root=Path(__file__).resolve().parents[1]
    for path in ('solver.py','coupling/domain_newton.py'):
        text=(root/path).read_text()
        assert 'BoundOccupationError' not in text and 'AmbiguousOccupation' not in text
        assert 'bound_diagnostics' not in text


def test_native_boundary_profile_changes_data_not_stiffness(tmp_path):
    pytest.importorskip('gmsh');pytest.importorskip('devsim')
    from geometry.build import build_mesh
    from poisson.device import Poisson
    c=Config(ls=2,ld=2,lg=2,spacer=1,width=2,thickness=2,sd_width=4,sd_thickness=4,tox=.5,mesh_h=2,channel_h=1,junction_h=1,dx=1,dy=1,dz=1)
    p=Poisson(c,build_mesh(c,tmp_path))
    try:
        K=p.K.copy();p.set_charge(p.donors);p.solve(.5,0)
        profile=lambda xyz:.003*xyz[:,1]+.005*xyz[:,2]
        p.set_boundary_profile(profile);p.solve(.5,.03)
        assert (p.K-K).nnz==0
        for x,bias in [(0,0),(c.length,.03)]:
            m=np.isclose(p.allxyz[:,0],x,atol=1e-8)
            np.testing.assert_allclose(p._get('Potential')[m],profile(p.allxyz[m])+bias,atol=1e-12)
    finally:p.close()
