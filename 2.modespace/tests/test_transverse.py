"""Physical numerical checks for transverse mode-space basis (no Gmsh/DEVSIM)."""
import unittest
import tempfile
import sys
from pathlib import Path
from dataclasses import replace
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import run as s
Config, CONST, build_mode_slices, H = s.load_reference(ROOT.parent/'1.example_build')


class TestStep02(unittest.TestCase):
    def setUp(self):
        self.cfg=Config()

    def test_exact_dirichlet_fd_eigenvalues(self):
        cfg=replace(self.cfg,ny=11,nz=7,modes=5)
        for v in cfg.valleys:
            es,us=s.lowest_eigenpairs(H(cfg,v,np.zeros((cfg.ny,cfg.nz))),cfg.modes)
            exact=s.sorted_analytic_levels(cfg,v,CONST)[:cfg.modes]
            self.assertLess(np.max(np.abs(es-exact)),1e-11)

    def test_orthogonality_and_hermitian(self):
        cfg=replace(self.cfg,ny=11,nz=7,modes=4)
        v=cfg.valleys[0]
        yy,zz=np.meshgrid(*s.axes(cfg)[1:],indexing='ij')
        h=H(cfg,v,.03*np.exp(-yy**2/20-zz**2))
        e,u=s.lowest_eigenpairs(h,cfg.modes)
        self.assertLess(np.max(np.abs(h-h.T.conj())),1e-13)
        self.assertLess(np.max(np.abs(u.T@u-np.eye(cfg.modes))),1e-12)

    def test_potential_energy_shift(self):
        cfg=replace(self.cfg,ny=9,nz=7,modes=4)
        v=cfg.valleys[0]
        e0,_=s.lowest_eigenpairs(H(cfg,v,np.zeros((cfg.ny,cfg.nz))),cfg.modes)
        e1,_=s.lowest_eigenpairs(H(cfg,v,.15*np.ones((cfg.ny,cfg.nz))),cfg.modes)
        np.testing.assert_allclose(e0-e1,.15,atol=1e-12)

    def test_thinner_sheet_increases_confinement(self):
        c1=replace(self.cfg,ny=9,nz=9,modes=2,thickness_nm=5.)
        c2=replace(self.cfg,ny=9,nz=9,modes=2,thickness_nm=3.)
        v=c1.valleys[0]
        e1,_=s.lowest_eigenpairs(H(c1,v,np.zeros((9,9))),2)
        e2,_=s.lowest_eigenpairs(H(c2,v,np.zeros((9,9))),2)
        self.assertGreater(e2[0],e1[0])

    def test_neighbour_overlap_flat_potential(self):
        cfg=replace(self.cfg,nx=7,ny=7,nz=5,modes=3)
        v=cfg.valleys[0]
        e,u,t,onsite,hop=build_mode_slices(cfg,v,np.zeros((cfg.nx,cfg.ny,cfg.nz)))
        o=u[0].T@u[1]
        self.assertLess(np.max(np.abs(o-np.eye(cfg.modes))),1e-12)
        self.assertLess(np.max(np.abs(hop[0]+t*np.eye(cfg.modes))),1e-12)
        self.assertLess(np.max(np.abs(onsite[0]-onsite[0].T.conj())),1e-14)

    def test_import_step01_poisson_npz(self):
        cfg=replace(self.cfg,nx=7,ny=5,nz=3)
        x,y,z=s.axes(cfg)
        values=x[:,None,None]*.01 + y[None,:,None]*.001 + z[None,None,:]*.003
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'poisson_only.npz'
            np.savez(p,x_nm=x,y_nm=y,z_nm=z,phi_v=values)
            imported,tag=s.load_phi(p,cfg)
            np.testing.assert_allclose(imported,values,atol=1e-13)
            self.assertIn('Poisson',tag)

    def test_density_normalization(self):
        cfg=replace(self.cfg,ny=9,nz=7,modes=2)
        v=cfg.valleys[0]
        _,u=s.lowest_eigenpairs(H(cfg,v,np.zeros((9,7))),2)
        dy=cfg.width_nm/(cfg.ny+1);dz=cfg.thickness_nm/(cfg.nz+1)
        probability=np.abs(u[:,0])**2/(dy*dz)
        self.assertAlmostEqual(float(probability.sum()*dy*dz),1.0,places=12)


if __name__=='__main__': unittest.main()
