# DEVSIM-QT STEP09 modularization

This is the modularized STEP09 working copy of STEP08. The SCF path calls a replaceable transport interface; Poisson and mapping stay independent of transport algorithms. The existing MS, NEGF, Spectral and physics-validation implementations are preserved. Physics validation status is unchanged: biased GAA fully self-consistent I–V remains quantitatively unvalidated.

## Run

Install the dependencies listed in `requirements.txt`, plus native libraries required by Gmsh and DEVSIM. Existing invocations remain supported, for example:

```bash
python run.py point --backend negf --out results/example --vg 0.4 --vd 0.05
python run.py point --backend ms --out results/example_ms --vg 0.4 --vd 0.05
python run.py point --backend spectral --out results/example_spectral --vg 0.4 --vd 0.05
```

Run tests with `pytest -q`. See `ARCHITECTURE.md` for the new extension seam and `MODULARIZATION_REPORT.md` for migration status and checks actually run. Do not treat a successful software run as quantitative validation of device physics.

To add a transport model, implement the interface and register its factory in `transport/backends/__init__.py`; see `ARCHITECTURE.md`. The registry supplies CLI choices and configuration validation.
