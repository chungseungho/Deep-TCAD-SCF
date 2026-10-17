"""Fast structural checks for the modular transport boundary."""
from pathlib import Path

from transport.registry import create_transport_backend, registered_backends, register_backend


def test_builtin_registry_and_backend_protocol():
    from config import Config
    assert registered_backends() == ('ms', 'negf', 'spectral')
    for name in registered_backends():
        backend=create_transport_backend(name,Config())
        for method in ('prepare','calibrate_contacts','solve','convergence_metrics','validate'):
            assert callable(getattr(backend,method))


def test_poisson_has_no_transport_backend_dependency():
    root=Path(__file__).parents[1]
    for path in (root/'poisson').glob('*.py'):
        assert 'transport.backends' not in path.read_text()


def test_registry_accepts_future_backend_without_scf_edit():
    class FutureBackend: pass
    register_backend('test_future',lambda config:FutureBackend())
    assert isinstance(create_transport_backend('test_future',None),FutureBackend)

