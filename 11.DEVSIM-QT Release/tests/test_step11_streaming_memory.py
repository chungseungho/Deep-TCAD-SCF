import tracemalloc
import numpy as np
from transport.quadrature import integrate_vector


def test_many_component_integral_does_not_store_all_energy_or_panel_values():
    # The previous vector algorithm retained >100 MiB for this modest fixture.
    # This is a measured resource regression, not a mock of the implementation.
    tracemalloc.start()
    try:
        result,diag=integrate_vector(lambda e:np.repeat((e*e)[:,None],1024,axis=1),np.linspace(0,1,401),8,1e-7,1e-10,20000)
        _,peak=tracemalloc.get_traced_memory()
    finally:tracemalloc.stop()
    np.testing.assert_allclose(result,1/3,rtol=1e-9)
    assert peak<8*1024*1024
