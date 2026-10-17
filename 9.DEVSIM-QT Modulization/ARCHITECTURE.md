# DEVSIM-QT modular architecture

## Dependency direction

The runtime path is `Gmsh geometry → DEVSIM Poisson → coupling.Mapping → coupling.scf_controller → QuantumTransportBackend → Mapping.charge → Poisson`. The Poisson equations and the established potential and charge mappings are retained. `poisson/` contains no transport imports; backend modules receive quantum potential, bias, grid/builders and configuration through a context, and never receive a DEVSIM Poisson object.

`solver.scf` owns run bookkeeping, charge mixing, physical Poisson residual checks, checkpointing and final file coordination. It holds one backend instance and calls its protocol hooks for contact calibration, transport solve, validation, convergence metrics, finalization and backend-specific exported fields. It contains no MS/NEGF/Spectral name tests or equations. `coupling.scf_controller.evaluate_transport_step` implements one backend-agnostic Poisson–mapping–transport–mapping–Poisson step; the mock test runs this production function.

## Backend boundary

`transport.api.QuantumTransportBackend` specifies the adapter contract. `transport.results.TransportResult` defines density, drain current, bond current, source/drain chemical potentials, diagnostics and an optional payload for internal state required by convergence or serialization. The built-in implementations live in `transport/backends/ms.py`, `negf.py` and `spectral.py`. Their occupation tracking, backend convergence rules, contact calibration, MS open-current postprocessing and result fields live with the backend, aided by `transport/backend_hooks.py`. Physics formulas are carried over from STEP08.

## Add a backend

1. Add `transport/backends/<name>.py` with a class implementing `QuantumTransportBackend`; return a `TransportResult` from `solve` and implement all convergence, finalization, contact-calibration and export hooks.
2. Register it in `transport/backends/__init__.py` with `register_backend('<name>', <Class>)`.
3. Select it in a configuration using `backend: '<name>'` or `--backend <name>`. Config validation and CLI choices are populated from the registry.

No edits are required in `solver.py`, `coupling/scf_controller.py`, `poisson/device.py`, `coupling/mapping.py`, or existing backend implementations to add a backend. Add backend-specific configuration fields only when that physical method needs them; retain defaults so existing configurations remain valid.

## Preserved physical semantics

The Hamiltonian, valley and spin multiplicities, reservoir filling, density/current normalization, Poisson assembly, boundary conditions, mapping operators and SCF mixing/Newton equations remain as in STEP08. MS density remains closed-system stationary occupation; its separately postprocessed open scattering current is explicitly owned by `MSBackend`. Existing unresolved bound-state occupation and incomplete biased GAA quantitative validation are unchanged.
