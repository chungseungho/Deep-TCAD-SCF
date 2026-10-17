# Reproduce STEP10

## Environment

Tested: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, pytest 9.1.1, Gmsh 4.15.2, DEVSIM 2.11.0.

```bash
python -m pip install -r requirements.txt
```

Gmsh requires native Xft/GLU libraries. DEVSIM requires BLAS/LAPACK. On normal Ubuntu these can be installed with the system package manager (`libxft2`, `libglu1-mesa`, `libopenblas0-pthread` and dependencies). Existing locally extracted libraries were used in this Work runtime.

Set `QT_NATIVE_LIB_DIR` to the directory containing `libXft.so.2`, with an `openblas-pthread/` child when using extracted libraries. The execution wrapper discovers this session's existing `../native_libs/root/usr/lib/x86_64-linux-gnu` if present. Otherwise use the normal system `LD_LIBRARY_PATH` and `DEVSIM_MATH_LIBS`. Tested settings were `DEVSIM_MATH_LIBS=libopenblas.so.0` and, separately, `libblas.so.3:liblapack.so.3`.

## Tests and independent benchmarks

```bash
python validation/execute.py --log-dir results/reproduce_tests -- python -m pytest -v
python validation/execute.py --log-dir results/reproduce_independent/run -- python validation/step10_independent.py results/reproduce_independent
```

The wrapper records source snapshots, hashes, full command/log, runtime, exit code and source stability. Exit 124 is a timeout, never a converged state. Use a fresh output directory for each changed source/config.

## Actual 3D calculations

```bash
python validation/execute.py --log-dir results/reproduce_raised/run --timeout 180 -- python run.py point --config configs/step10_raised.json --out results/reproduce_raised --vg 0.4 --vd 0
python validation/execute.py --log-dir results/reproduce_uniform/run --timeout 180 -- python run.py point --config configs/step10_uniform.json --out results/reproduce_uniform --vg 0.4 --vd 0
python validation/execute.py --log-dir results/reproduce_uniform_001/run --timeout 180 -- python run.py point --config configs/step10_uniform.json --out results/reproduce_uniform_001 --vg 0.4 --vd 0.01 --initial-state results/reproduce_uniform/states/vg0.400000_vd0.000000.npz --max-iterations 24
```

Exact commands for mesh refinement, longer contacts, restart and bounded failed attempts are in `results/final_*/run/execution.json`. To reproduce a particular code version, extract that run's `source_snapshot.zip` into a separate checkout and use its recorded config/command.

## Numerical diagnostics and plots

```bash
python validation/execute.py --log-dir results/reproduce_modes/run --timeout 300 -- python validation/step10_frozen_audit.py results/final_raised_equilibrium/states/vg0.400000_vd0.000000.npz results/reproduce_modes
python validation/execute.py --log-dir results/reproduce_leads/run -- python validation/step08_uniform.py results/reproduce_leads
python validation/execute.py --log-dir results/reproduce_energy/run --timeout 180 -- python results/energy_convergence/reproduce.py
python results/assessment/reproduce_analysis.py
python results/assessment/reproduce_charge_current.py
python results/assessment/reproduce_attempt_flags.py
```

New checkpoints can be supplied directly to `--initial-state results/<run>/checkpoints/<bias-key>.npz`; their coordinate metadata is checked by the existing loader. The legacy conversion for the earlier 0.05 V run is retained separately with an exact mesh hash check.

Frozen-field diagnostics do not create new SCF solutions. Density/current need joint mode/grid/contact/energy convergence at a biased representative point before quantitative certification. No state in this package is Teacher Data. No missing I–V point is fabricated.
