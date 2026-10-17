# Modularization report

## Source and changes

This package was copied from the locally available `restored_step08/8.DEVSIM-QT` tree. STEP01–08 repository files were not edited. The copy excludes bulky historical `results/` runs and cache directories; the executable source, configs, tests, and existing validation scripts are included. Added the common transport protocol/result schema, factory registry, three backend adapters, a controller-facing registry seam, and structural tests.

The STEP08 calculation evaluator and Poisson/mapping equations were preserved verbatim. Existing CLI options remain in place. Existing physical-status fields, including `quantitative_validation_pass: false`, are retained. This refactor does not resolve the documented bound-state occupation issue or claim that biased GAA fully self-consistent transport is quantitatively validated.

## Verification

- Python 3.12.14: available.
- `python -m compileall -q .`: PASS.
- `python run.py --help`: PASS; existing backend choices and commands are exposed.
- Registry smoke check: PASS for MS, NEGF, Spectral factories.
- Existing pytest suite and new structural pytest: NOT RUN; `pytest` is not installed in the available Python environment.
- Gmsh/DEVSIM runtime regression and STEP08-vs-STEP09 numerical comparison: NOT RUN. The environment lacks pytest and the project-local virtual environment; no numerical equivalence is asserted from static checks.

## Remaining work

The registry seam deliberately preserves numerical behavior but does not yet complete the requested full migration: backend-specific branches still reside in the preserved evaluator, SCF acceptance/postprocessing/export still inspect backend names, and `TransportResult` is wrapped for compatibility rather than being the sole object consumed by SCF. Completing those moves requires running the STEP08 numerical regression suite in its supported Gmsh/DEVSIM environment, especially for domain-Newton recovery and final physical-residual reevaluation. This package is an executable modularization foundation, not a verified claim that all Definition-of-Done gates passed.
