"""Fast structural checks for the modular transport boundary."""
from pathlib import Path
import pytest

from transport.registry import create_transport_backend, registered_backends, register_backend


def test_builtin_registry_and_backend_protocol():
    from config import Config
    from transport.api import QuantumTransportBackend
    assert {'ms','negf','spectral'}.issubset(registered_backends())
    for name in ('ms','negf','spectral'):
        backend=create_transport_backend(name,Config())
        assert isinstance(backend,QuantumTransportBackend)
        for method in ('prepare','calibrate_contacts','solve','convergence_metrics','validate'):
            assert callable(getattr(backend,method))


def test_poisson_has_no_transport_backend_dependency():
    root=Path(__file__).parents[1]
    for path in (root/'poisson').glob('*.py'):
        assert 'transport.backends' not in path.read_text()


def test_registry_accepts_future_backend_without_scf_edit():
    import transport.registry as registry
    class FutureBackend: pass
    register_backend('test_future',lambda config:FutureBackend())
    assert isinstance(create_transport_backend('test_future',None),FutureBackend)
    registry._FACTORIES.pop('test_future')


def test_mock_quantum_runs_through_common_scf_controller():
    import numpy as np
    from config import Config
    from coupling.scf_controller import evaluate_transport_step
    from transport.backends.mock_quantum import register_mock_backend
    import transport.registry as registry
    from transport.registry import create_transport_backend

    class Poisson:
        def set_charge(self,value):self.charge=np.asarray(value)
        def solve(self,vg,vd):return np.full_like(self.charge,float(vg))
    class Mapping:
        def potential(self,value):return np.asarray(value)
        def charge(self,value):return np.asarray(value)
    register_mock_backend()
    backend=create_transport_backend('mock_quantum',Config(backend='mock_quantum').validate())
    poisson=Poisson();mapping=Mapping();grid=object()
    backend.prepare({'grid':grid,'mapping':mapping})
    mu,calibration=backend.calibrate_contacts({'grid':grid})
    assert mu==pytest.approx(.2) and calibration['status']=='mock_calibrated'
    step=evaluate_transport_step(poisson,mapping,backend,[],grid,np.ones(4),.2,.1,.05)
    assert backend.solve_calls==1
    assert np.allclose(step.transport_result.electron_density_cm3,1.0)
    assert np.allclose(step.mapped_density,step.transport_result['n'])
    assert step.poisson_potential.shape==(4,) and step.poisson_image.shape==(4,)
    assert step.transport_result.mu_right_ev==pytest.approx(.15)
    metrics=backend.convergence_metrics(None,step.transport_result)
    assert backend.converged(1,metrics,step.potential_residual,0.,0.)
    assert backend.finalize(step.transport_result,.05,grid,True)['transport_status']=='CALCULATED'
    assert backend.export_fields(step.transport_result)=={}
    registry._FACTORIES.pop('mock_quantum')


def test_mock_quantum_completes_solver_scf_without_poisson_runtime(tmp_path,monkeypatch):
    import numpy as np
    from config import Config
    from transport.backends.mock_quantum import register_mock_backend
    from transport.registry import create_transport_backend
    import solver
    import analysis.validity
    import export_fields
    import device.regions
    import transport.registry as registry

    class Poisson:
        donors=np.ones(4);volume=np.ones(4);mesh_path=tmp_path/'mesh.msh'
        xyz={'silicon':np.zeros((4,3))}
        def set_charge(self,value):self.charge=np.asarray(value)
        def solve(self,vg,vd):return np.full(4,.1)
        def gauss_balance(self):return {'gauss_relative_error':0.}
        def export(self,path):pass
    class Mapping:
        def potential(self,value):return np.asarray(value)
        def charge(self,value):return np.asarray(value)
    class Grid:
        xyz=np.zeros((4,3));x=np.arange(4,dtype=float);y=np.arange(1,dtype=float);z=np.arange(1,dtype=float)
        mask=np.ones((4,1,1),dtype=bool);cell_volume=1.
    register_mock_backend()
    c=Config(backend='mock_quantum',max_iterations=3).validate()
    backend=create_transport_backend('mock_quantum',c)
    monkeypatch.setattr(solver,'diagnose',lambda *a,**k:{'reservoir_screening_diagnostic_pass':True})
    monkeypatch.setattr(analysis.validity,'element_fields',lambda *a,**k:{})
    monkeypatch.setattr(analysis.validity,'evaluate',lambda *a,**k:{'physical_status':'NOT_EVALUATED'})
    monkeypatch.setattr(export_fields,'export_state',lambda *a,**k:None)
    poisson=Poisson()
    result,_,_=solver.scf(c,poisson,Mapping(),[],Grid(),.2,.1,.05,tmp_path,'mock',backend=backend)
    assert result['numerical_status']=='SCF_CONVERGED'
    assert backend.solve_calls==2  # one iteration plus the mandated final reevaluation
    registry._FACTORIES.pop('mock_quantum')


def test_common_scf_controller_has_no_backend_name_switches():
    root=Path(__file__).parents[1]
    for rel in ('solver.py','coupling/scf_controller.py'):
        source=(root/rel).read_text()
        assert "backend ==" not in source and "backend=='" not in source
        assert 'c.backend' not in source

