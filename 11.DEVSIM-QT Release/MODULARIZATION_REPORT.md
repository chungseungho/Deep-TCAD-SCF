> Historical baseline document retained for research provenance. Its results and paths are not STEP11 validation evidence. See RESULTS.md and RELEASE_CERTIFICATE.json for this candidate.

# STEP09 modularization report

## Refactor completed

The existing STEP09 working tree was continued in place. The previous implementation left model branches and a legacy evaluator in `solver.py`; those were removed. MS, NEGF and Spectral density calculations now run from their respective backend classes. The common SCF path calls the transport protocol for solve, validation, state/convergence metrics, finalization/current checks, metadata and export fields. MS occupation matching, MS current postprocessing, backend diagnostics, backend contact calibration and backend-specific state export are owned by the backend layer. `coupling.scf_controller.evaluate_transport_step` is the production shared-map step and has no backend-specific logic.

Poisson and mapping implementations were carried forward unchanged. `Config` validation and CLI backend choices are sourced from the registry. The mock backend is registered only in its test and is exercised through the production SCF controller function.

## Tests actually run

- STEP08 baseline checkout: `python -m pytest -v` — **60 passed, 3 skipped**.
- Updated STEP09 tree: `python -m pytest -v` — **67 passed, 3 skipped**.
- The three skips in each tree are `test_devsim_integration.py`, `test_shared_poisson_capacitor.py`, and `test_triple_junction.py`. All skip because the Gmsh Python module cannot load native `libXft.so.2` in this runtime.
- New tests cover built-in registry/API compliance, a future registry entry, structural dependency direction, `mock_quantum` through both the shared SCF step and the complete `solver.scf` iteration/final reevaluation/export path, and STEP08 numerical goldens for all three transport backends.
- Numerical regression uses the same reduced 3D GAA quantum grid, zero mapped potential, fixed source chemical potential −0.1 eV and zero drain bias. For MS, NEGF and Spectral, the maximum absolute difference was **0** for electron density at every quantum cell, bond current, source chemical potential, tracking signature, final current (including the independent MS postprocess), postprocessed bond currents, and input mapped potential. Golden data and its fixture are in `validation/modularization_reference.npz` and `validation/modularization_reference.json`.
- Existing mapping conservation/work tests passed in both suites. Full STEP08-vs-STEP09 Gmsh/DEVSIM Poisson potential and end-to-end SCF comparisons were **NOT RUN** because Gmsh cannot import. Thus the quantum backend regression is exact for this fixture, while end-to-end electrostatic/SCF equivalence remains unverified.

## Environment

Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and pytest 9.1.1 were used. The `gmsh==4.15.2` and `devsim==2.11.0` wheels installed. DEVSIM imports when `DEVSIM_MATH_LIBS='libblas.so.3:liblapack.so.3'` is set. Gmsh import remains blocked by missing system `libXft.so.2`; installing it with apt was blocked by the execution container's `setgroups`/privilege restrictions. Native Poisson tests therefore remain skipped, not passed.

## STEP10 replacement implications

To replace the NEGF algorithm, add or replace a backend implementation and registry entry, then select it in config/CLI. The common SCF controller, Poisson device, and mapping code require no modification. Current validation does not claim that a replacement algorithm has been implemented, nor that STEP08's unresolved bound-state issue or biased GAA quantitative-validation status has changed. `quantitative_validation_pass` remains false.
