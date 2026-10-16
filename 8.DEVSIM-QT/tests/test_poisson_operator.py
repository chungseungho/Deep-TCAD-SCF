import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
from poisson.discretization import assemble


def test_tetra_stiffness_energy_affine_volume():
    # Deliberately obtuse tetrahedron: signed edge conductances are necessary.
    nodes=np.array([[0.,0,0],[2,0,0],[.2,.1,0],[.3,.02,.15]])
    K,V=assemble(nodes,np.array([[0,1,2,3]]))
    vol=abs(np.linalg.det(np.column_stack([nodes[1]-nodes[0],nodes[2]-nodes[0],nodes[3]-nodes[0]])))/6
    assert_allclose(V.sum(),vol,rtol=1e-14)
    assert_allclose(K.toarray(),K.toarray().T,atol=1e-14)
    assert_allclose(np.asarray(K.sum(axis=1)),0,atol=1e-14)
    grad=np.array([.3,-.1,.7]);phi=nodes@grad+2
    assert_allclose(phi@K@phi,vol*np.dot(grad,grad),rtol=1e-12)
    assert np.linalg.eigvalsh(K.toarray()).min()>-1e-12


def test_internal_shared_face_flux_cancels():
    nodes=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
    K,V=assemble(nodes,np.array([[0,1,2,3],[0,2,1,4]]))
    phi=np.array([.2,.7,-.2,.3,-.5]);charge=K@phi
    assert abs(charge.sum())<1e-14
    assert_allclose(V.sum(),1/3,rtol=1e-14)
