# DEVSIM-QT STEP09 modularization

This is a modularized working copy of STEP08. It retains the existing Poisson, mapping, MS, NEGF, Spectral, SCF, CLI, and validation implementations and adds a backend protocol, common result schema, registry, and backend adapters. Physics validation status is unchanged: biased GAA fully self-consistent I–V remains quantitatively unvalidated.

## Run

Install the dependencies listed in `requirements.txt`, plus the native Gmsh and DEVSIM runtime dependencies required by the device solver. Existing invocations remain supported, for example:

```bash
python run.py point --backend negf --out results/example --vg 0.4 --vd 0.05
python run.py point --backend ms --out results/example_ms --vg 0.4 --vd 0.05
python run.py point --backend spectral --out results/example_spectral --vg 0.4 --vd 0.05
```

Run tests with `pytest -q`. See `ARCHITECTURE.md` for the new extension seam and `MODULARIZATION_REPORT.md` for migration status and checks actually run. Do not treat a successful software run as quantitative validation of device physics.
