"""Extensible transport backend registry; no SCF edits required to register."""
from collections.abc import Callable
from typing import Any

_FACTORIES: dict[str, Callable[[Any], Any]] = {}


def register_backend(name: str, factory: Callable[[Any], Any], *, replace: bool = False) -> None:
    key = name.strip().lower()
    if not key:
        raise ValueError('Backend name cannot be empty')
    if key in _FACTORIES and not replace:
        raise ValueError(f'Backend already registered: {key}')
    _FACTORIES[key] = factory


def create_transport_backend(name: str, config: Any):
    # Importing this package registers the built-in adapters once.
    import transport.backends  # noqa: F401
    key = name.strip().lower()
    try:
        return _FACTORIES[key](config)
    except KeyError as exc:
        raise ValueError(f'Unsupported DEVSIM-QT backend: {key}') from exc


def registered_backends() -> tuple[str, ...]:
    import transport.backends  # noqa: F401
    return tuple(sorted(_FACTORIES))

