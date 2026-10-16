import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
from transport.quadrature import integrate_vector


def test_ten_picovolt_width_resonance_resolved_without_eta():
    centre=.123456789;width=1e-11
    # Independent analytic Lorentzian area, including finite endpoints.
    fun=lambda e:(width/((e-centre)**2+width**2))[:,None]
    breaks=np.unique(np.r_[-1,1,centre,centre+width*np.array([-1e5,-1e3,-10,-1,1,10,1e3,1e5])])
    value,report=integrate_vector(fun,breaks,8,2e-5,1e-8,120000)
    exact=np.arctan((1-centre)/width)-np.arctan((-1-centre)/width)
    assert_allclose(value[0],exact,rtol=2e-5)
    assert report['maximum_normalized_error_estimate']<=1
