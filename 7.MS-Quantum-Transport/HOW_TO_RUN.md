# How to run STEP07

Run commands from `7.MS-Quantum-Transport/`. The package is self-contained;
STEP01–06 are not imported at runtime. Small prior-state inputs actually used
for regression are copied under `results/inputs/`, with hashes. A seed is not
an accepted solution for a new Hamiltonian or ensemble.

## Private CPU environment

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
.venv/bin/python validation/run_tests.py
```

Linux Gmsh requires its normal OpenGL/X shared-library dependencies; DEVSIM
requires a compatible BLAS/LAPACK. If those are already installed outside the
default loader path, set LD_LIBRARY_PATH to that installation and optionally
DEVSIM_MATH_LIBS=libopenblas.so.0. Do not modify a shared Python/native setup.
The Work run used private native dependencies; they are not vendored in this ZIP.
The test runner isolates native DEVSIM/Gmsh checks from the quantum tests and
requires JUnit completion files. A mixed-process run on this runtime exited
before its completion summary, so its zero exit status was not counted as a
test pass. Both isolated groups must complete; logs are retained.

## Single point and continuation

```bash
.venv/bin/python run.py point --config configs/gaa_ms.json --vg .8 --vd 0 --out results/my_ms_eq
.venv/bin/python run.py point --config configs/gaa_negf.json --vg .8 --vd .1 --out results/my_negf
.venv/bin/python run.py point --config configs/gaa_spectral.json --vg .8 --vd .1 --out results/my_spectral
.venv/bin/python run.py idvgs --config configs/gaa_ms.json --vd .1 --gate-range 0:1.2:.05 --out results/my_transfer
.venv/bin/python run.py idvds --config configs/gaa_ms.json --gates .4,.8,1.2 --drain-range 0:1.2:.05 --out results/my_output
```

Range syntax is `start:stop:positive_step` or comma-separated values. Endpoints
are included when aligned with the step. Both gate and drain are supported
from 0 to 6 V. `--dry-run` reports configuration and requested point count
without generating a mesh or solving physics.

Accepted points supply density seeds for the next bias. Numerical SCF failures
can trigger bisection from the last accepted bias down to
`--minimum-bias-step`. Origin ambiguity is missing model evidence and is not
bypassed by changing the occupation rule. `--resume` reuses accepted saved
states only when output settings and source identity match. An explicit
`--initial-state path.npz` must use identical Poisson nodes. Use a fresh output
folder when changing configuration or source; existing results are protected.

The delivered MS transfer/output attempts contain rejected origins. Repeating
these commands is expected to reproduce failure statuses until the finite
reservoir/origin model is improved. Do not disable rejection to obtain a curve.

`potential_reference_v` shifts all electrostatic contact potentials together.
Automatically calibrated μ shifts by its negative. An explicitly supplied μ
must be shifted consistently by the caller. Frozen open-density/transmission
tests and an actual GAA SCF regression verify this gauge convention.

## Frozen and independent checks

```bash
.venv/bin/python validation/analytical.py results/my_analytical
.venv/bin/python run.py frozen --state results/gaa_ms_lowbias/states/vg0.800000_vd0.000000.npz --out results/my_frozen
.venv/bin/python validation/frozen_step06_regression.py results/inputs/device_lowbias.npz results/my_repaired_negf
.venv/bin/python validation/audit_saved_bound_states.py results/gaa_spectral_regression/states/vg0.800000_vd0.100000.npz results/my_bound_audit.json
.venv/bin/python validation/frozen_convergence.py results/gaa_ms_lowbias/states/vg0.800000_vd0.000000.npz results/my_subspace
.venv/bin/python validation/reduced_realspace.py results/my_reduced_reference
.venv/bin/python validation/numerical_campaign.py results/my_numerical_plan
```

The last command creates a plan/configs only. Add `--execute` to run the
listed fresh 3D calculations locally. It includes three levels for dx,
transverse spacing, mode counts, reservoir boundaries, contact lengths and
Poisson meshes, plus fixed/local and plain/nested comparisons. No plateau is
certified from an unexecuted plan. A frozen-field subspace study is not a
self-consistent convergence campaign.

## Stress tests and exports

The finite-number default and `configs/gaa_ms_fixed_mu.json` are different
ensembles. The latter matches the independently calibrated lead mu for a
same-contact comparison. A finite-N mu shift must not be presented as a
same-contact NEGF equivalence test.

```bash
.venv/bin/python run.py stress --config configs/gaa_ms.json --axis gate --gate-range 0:6:.25 --vd .1 --out results/my_gate_stress
.venv/bin/python run.py stress --config configs/gaa_ms.json --axis drain --drain-range 0:6:.25 --vg .8 --out results/my_drain_stress
```

Expect physical-scope warnings and possible electronic-model failures. These
are stress tests, not validated 6 V silicon/oxide characteristics. Each accepted
state exports quantum hexahedra in nm (`*_nm.vtu`) and native DEVSIM fields in
cm (`*_electrostatics_cm.vtk` and companion files). Potential is V and density
is cm^-3 in both. NPZ stores the canonical device coordinates and all numerical
settings. For world-axis geometry interpretation use `device/definition.py`.

## User-submitted batch jobs

`job.sh` is a template. It contains no cluster account, partition or automatic
submission. For example, the user can supply their own resource options to
`sbatch ... job.sh stage-a-negf`. Supported tasks are `test`, `stage-a-ms`,
`stage-a-negf`, `stress-gate`, `stress-drain`. STEP07_PYTHON can select an
existing private interpreter. No school HPC connection or job submission was
performed for this deliverable.

Exit code 0 means requested numerical solves completed, not that Level4
physical validation passed. Exit code 2 records failed/rejected points. Read
`numerical_status`, `transport_status`, `physical_status` and
`quantitative_validation_pass` separately.
