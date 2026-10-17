# DEVSIM-QT 1.0 Release Candidate

Contact-consistent ballistic silicon GAA development based on STEP10 commit `28e5a26c16ceedeb88602a0216a478e979fea247`. Existing STEP01–10 source is preserved.

**This is a Release Candidate. Device accuracy and the complete raised-S/D gate/drain sweeps are not certified. Teacher-ready remains false.** Read `RELEASE_CERTIFICATE.json` for the automated mandatory-gate result, `MODEL_SCOPE.md` for the physical contract and limits, and `RESULTS.md` for measured outcomes.

Implemented here: a full-basis transverse Schrödinger–Poisson electrode, common electrode boundary/self-energy data, rectangular basis-overlap coupling, block production Hamiltonians, banded resonance guidance, constant-memory adaptive integration, neutral failure contracts, strict checkpoint validation, and an environment/refinement/export/certificate CLI. Independent references explicitly test bound spectral residues, cut invariance and full real-space/contact density.

```bash
python -m pip install -r requirements.txt
python validation/release_cli.py doctor
python -m pytest -v
python validation/release_cli.py run -- point --config configs/release_raised.json --out output/eq --vg .4 --vd 0
python validation/release_cli.py certify
```

See `RUN_INSTRUCTIONS.md` for native libraries, sweep commands, resume and numerical refinement. A nonzero certificate exit code means the mandatory release criteria are not all PASS; it is expected for this candidate.

Primary measured data and logs are under `results/`. Failed attempts, resource stops, original source snapshots and scope limitations are retained. `results/assessment/*actual_iv.csv/png` contains only actually computed SCF points; a plot with equilibrium alone is not an I–V validation. No failed point is interpolated or fitted.

The main implementation boundary is `ContactModel / HamiltonianProvider / DensitySolver -> QuantumTransportBackend -> TransportResult -> common SCF`. New backends register in `transport/registry.py` through a backend factory and implement `transport/api.py`; they do not require Poisson, mapping or SCF edits. `mock_quantum` replacement is exercised by tests. MS/Spectral remain their original research models.

The current `configs/release_*.json` files are diagnostic refinement levels. The word release identifies this package, not a passed accuracy preset. `select-modes` proposes a starting basis; `refine` and representative-bias checks are still required.

No remote push or external HPC submission was performed. Dependency/source notices are in `DEPENDENCY_LICENSES.md`. Historical STEP10 documents are separated under `docs/baseline_history/` and are not STEP11 evidence.
