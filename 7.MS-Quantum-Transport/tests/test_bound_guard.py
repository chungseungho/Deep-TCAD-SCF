import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pytest
from validation.analytical import chain
from transport.bound_states import audit
from transport.negf_reference import solve_negf


def test_bias_window_bound_pole_requires_preparation_model():
    c,system,edge=chain(n=9,barrier=-.5)
    result=audit(system,c,edge,0.)
    assert len(result['bound_poles_below_both_leads'])>0
    e=result['bound_poles_below_both_leads'][0]['energy_ev'];mu=e+.05
    result=audit(system,c,mu,.1)
    assert result['bias_dependent_bound_occupation_unresolved']
    with pytest.raises(RuntimeError,match='BOUND_PREPARATION_REQUIRED'):
        solve_negf(system,c,mu,mu-.1,.1)
