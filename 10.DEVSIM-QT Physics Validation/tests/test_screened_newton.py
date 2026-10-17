import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.sparse import csr_matrix
from numpy.testing import assert_allclose
from coupling.mixing import screened_inverse


def test_screened_inverse_matches_analytical_local_response_with_oxide_dofs():
    K=np.array([[3.,-1.,0.,-.2],[-1.,4.,-1.,-.3],[0.,-1.,3.,-.1],[-.2,-.3,-.1,2.]])
    indices=np.array([0,2]);C=np.array([.7,0.,1.3,0.]);boundary=np.zeros(4,bool)
    J=-np.eye(2)-np.linalg.inv(K)[np.ix_(indices,indices)]@np.diag(C[indices])
    M=screened_inverse(csr_matrix(K),C,indices,boundary)
    for v in [np.array([.3,-.2]),np.array([1.,1.])]:
        assert_allclose(J@(M@v),v,atol=1e-13)
