import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
from transport.quadrature import integrate_vector


def test_integrable_band_edge_singularity_is_resolved():
    def f(e): return np.column_stack([1/np.sqrt(e),np.sqrt(e),np.exp(e)])
    value,diag=integrate_vector(f,[0.,1.],8,1e-9,1e-12)
    assert_allclose(value,[2.,2/3,np.e-1],rtol=1e-10)


def test_narrow_resonance_is_adaptively_resolved():
    width=.00013;centre=.231
    def f(e): return (width/((e-centre)**2+width**2))[:,None]
    value,diag=integrate_vector(f,np.linspace(0,1,51),8,2e-7,1e-10)
    truth=np.arctan((1-centre)/width)-np.arctan(-centre/width)
    assert_allclose(value[0],truth,rtol=2e-8)
