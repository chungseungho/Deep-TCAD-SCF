"""Stable boundary between quantum transport and the SCF framework.

Backends receive a device-independent context and return common observables.
Optional backend diagnostics live in ``TransportResult.optional``.
"""
from typing import Any, Protocol, runtime_checkable
from transport.results import TransportResult


@runtime_checkable
class QuantumTransportBackend(Protocol):
    name: str

    def prepare(self, context: Any) -> None: ...
    def calibrate_contacts(self, context: Any) -> dict[str, Any]: ...
    def solve(self, potential: Any, bias: float, context: Any) -> TransportResult: ...
    def convergence_metrics(self, previous: TransportResult | None,
                            current: TransportResult) -> dict[str, float]: ...
    def validate(self, result: TransportResult) -> None: ...

