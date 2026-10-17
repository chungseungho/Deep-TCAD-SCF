# DEVSIM-QT modular architecture

## Dependency direction

`run.py` and `config.py` retain the STEP08 CLI and configuration contract. Gmsh geometry and the DEVSIM Poisson implementation remain in `geometry/` and `poisson/`; quantum-grid interpolation and conservative charge transfer remain in `coupling/mapping.py` and `coupling/cell_integrals.py`. Transport algorithms remain under `transport/` and do not import DEVSIM objects.

The new `transport.api` protocol and `transport.results.TransportResult` define the intended replaceable boundary. `transport.registry` registers MS, NEGF, and Spectral adapters and supports future factories without changing registry callers. Adapter diagnostics and STEP08 internal observables are retained in the optional result payload so existing output/export code continues to work.

## Current migration seam

`solver.electronic` now resolves a backend from the registry. The adapter invokes the preserved STEP08 evaluator, which still contains the original numerical branches. This compatibility seam preserves physical formulas and saved output semantics while making the backend entry point explicit. It is an intermediate migration boundary: backend-specific SCF stability, MS open-current postprocessing, and export fields remain in `solver.py` and must move behind backend hooks before claiming complete SCF independence.

## Numerical invariants

The Poisson FEM/DEVSIM equations, boundary conditions, `phi_Q=P phi_P` potential mapping, integrated charge mapping, Hamiltonian construction, valley/spin multiplicities, reservoir occupations, energy integration, and MS postprocessing are carried over from STEP08 unchanged. The existing distinction between closed stationary MS density and separately postprocessed open scattering current remains in the output metadata.
