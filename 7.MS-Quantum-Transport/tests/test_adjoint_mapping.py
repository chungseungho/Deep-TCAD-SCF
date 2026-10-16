import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
from scipy.sparse import diags
from coupling.mapping import affine_matrix,Mapping,balanced_matrix


def test_affine_adjoint_is_positive_charge_and_work_conservative():
    nodes=np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]])
    sites=np.array([[.1,.2,.3],[.3,.15,.2],[.2,.2,.2]])
    p=affine_matrix(nodes,np.array([[0,1,2,3]]),nodes,sites)
    v=np.array([.04,.06,.04,1/6-.14]);w=np.full(3,1/18)
    m=Mapping(p,p.T@diags(w),v,w,nodes,sites,{'method':'adjoint'})
    n=np.array([.2,1.7,.7]);phi=nodes@np.array([.3,-.7,.2])+.4
    assert np.min(m.charge(n))>=0
    assert_allclose(m.charge(n)@v,n@w,atol=2e-16)
    assert_allclose((m.potential(phi)*n)@w,(phi*m.charge(n))@v,atol=2e-16)
    assert_allclose(m.potential(phi),sites@np.array([.3,-.7,.2])+.4,atol=2e-16)


def test_balanced_adjoint_preserves_both_constants_and_work():
    from scipy.spatial import cKDTree
    rng=np.random.default_rng(4);nodes=rng.random((13,3));sites=rng.random((17,3))
    nv=np.ones(13)/13;sv=np.ones(17)/17
    c,report=balanced_matrix(nodes,sites,nv,sv)
    p=diags(1/sv)@c.T;m=Mapping(p,c,nv,sv,nodes,sites,report)
    assert_allclose(m.charge(np.ones(17)),1,atol=3e-9)
    assert_allclose(m.potential(np.ones(13)),1,atol=3e-14)
    n=rng.random(17);phi=rng.random(13)
    assert_allclose((m.potential(phi)*n)@sv,(phi*m.charge(n))@nv,atol=2e-16)


def test_exact_tetra_cell_integral_has_both_marginals_and_linear_precision():
    from types import SimpleNamespace
    from itertools import product
    from coupling.cell_integrals import integrated_matrix
    from poisson.discretization import assemble
    nodes=np.array(list(product([0.,1.],repeat=3)))
    tet=np.array([[0,4,6,7],[0,6,2,7],[0,2,3,7],[0,3,1,7],[0,1,5,7],[0,5,4,7]])
    _,nv=assemble(nodes,tet)
    x=(np.arange(3)+.5)/3;y=(np.arange(4)+.5)/4;z=(np.arange(2)+.5)/2
    X,Y,Z=np.meshgrid(x,y,z,indexing='ij');sites=np.column_stack([X.ravel(),Y.ravel(),Z.ravel()])
    grid=SimpleNamespace(x=x,y=y,z=z,mask=np.ones((3,4,2),bool),xyz=sites,cell_volume=1/24)
    C,report=integrated_matrix(nodes,tet,nodes,nv,grid);sv=np.full(24,1/24)
    P=diags(1/sv)@C.T
    assert np.min(C.data)>0
    assert_allclose(C@np.ones(24),nv,atol=3e-13)
    assert_allclose(C.T@np.ones(8),sv,atol=3e-13)
    assert_allclose(P@nodes,sites,atol=3e-12)
    phi=nodes@np.array([.3,-.7,.2])+.8;n=np.arange(24)/24
    assert_allclose((P@phi*n)@sv,phi@(C@n),atol=2e-14)
