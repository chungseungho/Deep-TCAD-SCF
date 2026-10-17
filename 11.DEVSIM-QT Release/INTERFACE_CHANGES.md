# STEP11 interface changes

- `transport/contracts.py`: ContactModel, HamiltonianProvider and DensitySolver protocols; RGF density solver is selected inside NEGFBackend. Future predictor/kinetic-correction contracts are data-only placeholders.
- `transport/contact_model.py`: LeadReference owns profile, density, neutral μ, full eigenchannels and overlap embedding. The legacy contact model remains selectable for regression.
- `poisson/device.py`: generic `set_boundary_profile(provider)` changes end-face data only. It has no backend imports. Mapping and Poisson operators are unchanged.
- `modespace/hamiltonian.py`: System stores diagonal blocks/hops, optional explicit dense materialization. NEGF builders set block_only; finite MS/Spectral and small tests keep dense access.
- `QuantumTransportBackend.physical_validity`: model-scope interpretation is now backend-owned; common SCF only supplies Poisson field measurements.
- `transport/errors.py`: common TransportDomainError (MODEL_UNDERSPECIFIED) and TransportNumericalError. Common Newton/SCF do not import physical backend exception classes.
- `provenance.validate_restart`: source/config/node/density checks. Bias is intentionally independent of the config to support continuation.
- `validation/release_cli.py`: doctor, run, sweep, resume, select-modes, refine, compare, export, certify.
- `validation/release_certificate.py`: mandatory items use PASS/FAIL/NOT_RUN/INCONCLUSIVE/OUT_OF_SCOPE; RELEASE_READY is computed, never manually asserted.

To add a transport backend, implement every method in `transport/api.py`, return `TransportResult`, own calibration/domain checks/convergence/export metadata in the backend, register a factory with `register_backend(name, factory)` at application startup, and select that name in Config. A boundary provider can be installed from the backend's `prepare` hook when the environment supplies Poisson. The solver only catches common transport-domain exceptions. Test replacement through `mock_quantum` before device validation.

The generic future-learning contract requires physical raw residual verification of predictions; no learned predictor or DFT correction is present in this release candidate.
