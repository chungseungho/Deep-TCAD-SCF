# Install and reproduce

Python 3.12 was tested. Create a fresh environment and install pinned direct dependencies:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python validation/release_cli.py doctor --out environment.json
python -m pytest -v
```

On Linux, Gmsh may require system GLU/Xft libraries and DEVSIM requires a usable BLAS/LAPACK. These native libraries are not redistributed. If using an extracted native library tree, set `QT_NATIVE_LIB_DIR` to its `usr/lib/x86_64-linux-gnu` directory. The CLI/recorded executor adds it and its `openblas-pthread` subdirectory to the subprocess library path. Use `DEVSIM_MATH_LIBS=libopenblas.so.0` for that tree. `doctor` reports actual imports, not just wheel installation.

Run from this directory. Use fresh output directories when source/config changes:

```bash
python validation/release_cli.py run -- point --config configs/release_raised.json --out output/eq --vg .4 --vd 0
python validation/release_cli.py resume -- point --config configs/release_raised.json --out output/bias --vg .4 --vd .01 --initial-state output/eq/checkpoints/vg0.400000_vd0.000000.npz
```

Checkpoint reuse verifies source hash, all configuration fields, mesh coordinates and finite nonnegative density. `--resume` reuses converged points only within a matching output directory. A saved numerical iterate is a starting guess; final raw Poisson/charge residuals are reevaluated. Bias is separate from model config, so continuation can change bias.

Full requested sweep commands are available, but **the packaged configs have not passed the required device accuracy gates**. Do not label these commands' outputs certified merely because they converge:

```bash
python validation/release_cli.py sweep -- idvgs --config configs/release_raised.json --out output/idvg005 --vd .05 --gate-range 0:1.2:.05
python validation/release_cli.py sweep -- idvgs --config configs/release_raised.json --out output/idvg01 --vd .1 --gate-range 0:1.2:.05
python validation/release_cli.py sweep -- idvds --config configs/release_raised.json --out output/idvd --gates .4,.8,1.2 --drain-range 0:1.2:.05
```

No interpolation replaces a failed point. Requested and adaptive continuation points are distinguished. MODEL_UNDERSPECIFIED is not silently converted into a numerical continuation success.

```bash
python validation/release_cli.py select-modes --config configs/release_raised.json --out proposed.json
python validation/release_cli.py refine --config configs/release_raised.json --out output/modes --pairs 8/32,12/48,16/64 --vg .4 --vd .01
python validation/release_cli.py compare output/modes/modes_8_32/states/vg0.400000_vd0.010000.npz output/modes/modes_16_64/states/vg0.400000_vd0.010000.npz --out comparison.json
python validation/release_cli.py export --state output/eq/states/vg0.400000_vd0.000000.npz --out output/export --transmission
python validation/release_cli.py certify
```

Exports: compressed NPZ, CSV, and VTU/VTK. Transmission export is an energy sample of the stored state; it is not a new SCF or an integration-convergence certificate. Regions are compared separately; physical geometry changes are rejected by the numerical comparison command.

For recorded execution with source snapshot, command, runtime, peak child RSS and complete log:

```bash
python validation/execute.py --log-dir results/my_run -- python run.py point --config configs/release_raised.json --out output/my_run --vg .4 --vd 0
```

The executor has **no default walltime cutoff**. An explicitly supplied `--timeout` is recorded as TIMEOUT and is not convergence. Operating-system resource stops are documented separately. Avoid concurrent high-mode runs until the local RAM/CPU budget is measured.

`validation/revalidate_equilibrium_seed.py` is a restricted, explicit migration utility used during this work: only equilibrium states, explicitly recorded old/new source hash pairs, and a successful complete equilibrium density reevaluation are accepted. It produces an uncertified initial guess and retains original provenance; it does not bypass normal restart validation or modify the original state.

Reproduce standalone validations with `validation/step11_lead_convergence.py`, `validation/step11_frozen_modes.py`, `validation/step11_bound_localization.py` and the commands in each `results/*_run/execution.json`. STEP10-derived reference tests are retained; no new STEP11 result is replaced by a STEP10 device state.
