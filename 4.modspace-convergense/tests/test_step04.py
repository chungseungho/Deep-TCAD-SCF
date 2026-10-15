"""CPU-only physics/numerics tests, independent of gmsh and DEVSIM."""
import sys
import unittest
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
from diagnostics import manufactured_cube, direct_reference
from quantum import spectral_batched, lead_matrix_self_energy

class TestStep04(unittest.TestCase):
    def test_affine_exact_on_original_tetra(self):
        t=manufactured_cube(); self.assertLess(t.diagnostics()['linear_phi_error_V'],1e-12)
    def test_charge_mass_projection(self):
        t=manufactured_cube()
        rng=np.random.default_rng(99)
        n=rng.random(t.shape)*2e18
        stat=t.totals(n,t.to_nodes(n))
        self.assertLess(stat['relative_charge_mismatch'],1e-12)
    def test_constant_and_linear_interpolation(self):
        t=manufactured_cube();ones=t.to_regular(np.ones(8))
        self.assertLess(abs(ones-1).max(),1e-13)
    def test_exact_surface_matrix_quadratic(self):
        t=.38;eta=1e-4;H=np.array([[1.2,.11],[.11,1.46]],complex)
        E=np.linspace(.2,2.3,21)
        S,_=lead_matrix_self_energy(E,H,t,eta)
        residual=[]
        for e,ss in zip(E,S):
            rr=ss@((e+1j*eta)*np.eye(2)-H-ss)-t*t*np.eye(2)
            residual.append(np.max(abs(rr)))
        self.assertLess(max(residual),1e-11)

    def test_fullmatrix_lead_retarded_and_passive(self):
        H=np.array([[1.2,.2],[.2,1.4]])
        E=np.array([.4,.8,1.2,1.8]);s,g=lead_matrix_self_energy(E,H,.45,1e-6)
        self.assertLess(np.max(s.imag.diagonal(axis1=1,axis2=2)),1e-10)
        self.assertGreater(np.min(np.linalg.eigvalsh(g)),-1e-10)
    def test_RGF_dense_transmission_and_injectivity(self):
        rng=np.random.default_rng(31);n=7;m=3;t=.3;eta=1e-4
        base=np.diag(np.array([.2,.45,.6])+2*t).astype(complex)
        base[0,1]=base[1,0]=.08
        h=np.tile(base,(n,1,1))
        for j in range(1,n-1):
            a=rng.normal(size=(m,m))*.01;a=(a+a.T)/2
            h[j]+=a
        E=np.linspace(.2,1.1,12)
        T,AL,AR=spectral_batched(h,t,base,base,E,eta,True)
        TD,ADL,ADR=direct_reference(h,t,base,base,E,eta)
        self.assertLess(np.max(abs(T-TD)),1e-8)
        self.assertLess(np.max(abs(AL-ADL)),1e-8)
        self.assertLess(np.max(abs(AR-ADR)),1e-8)
    def test_fixed_energy_nodes_anchored(self):
        from quantum import energy_grid
        class C:
            efermi_source_ev=.3;temperature_k=300
        class B:
            hopping_ev=.6
        E1,_=energy_grid(C(),.1,B(),np.array([.0]),np.array([-.1]),.005)
        E2,_=energy_grid(C(),.1,B(),np.array([.00001]),np.array([-.10001]),.005)
        shared=np.intersect1d(E1,E2)
        self.assertGreater(len(shared),len(E1)-4)
        self.assertTrue(np.allclose(np.diff(E1),.005,atol=1e-10))

if __name__=='__main__':unittest.main()

class TestGmshPhysicalConnectivity(unittest.TestCase):
    def test_msh2_parser_selects_si_tetra(self):
        import tempfile
        from transfer import load_si_tetra_msh2
        txt='''$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
2
3 4 "Silicon"
3 8 "GateOxide"
$EndPhysicalNames
$Nodes
4
11 0 0 0
12 1 0 0
13 0 1 0
14 0 0 1
$EndNodes
$Elements
1
1 4 2 4 7 11 12 13 14
$EndElements
'''
        with tempfile.TemporaryDirectory() as tmp:
            f=Path(tmp)/'mesh.msh';f.write_text(txt)
            x,t=load_si_tetra_msh2(f)
        self.assertEqual(x.shape,(4,3))
        self.assertTrue(np.array_equal(t,np.array([[0,1,2,3]])))
