# DEVSIM-QT STEP10 — physical validation

**Partial success. `DEVICE_SIMULATION_VERIFIED=false`.**

This package includes the complete STEP09 modular source, STEP10 independent transport benchmarks, a verified block-inertia bound-pole audit, real Gmsh/DEVSIM–NEGF SCF runs, and failed-convergence evidence. No validated Raised S/D I–V or training-ready data is claimed.

- Native tests: 74 passed.
- Independent wave-matching, dense-GF, full-3D and multimode cases: 13 passed.
- Real Raised S/D and uniform GAA equilibrium SCF: converged.
- Uniform GAA Vg=0.4 V, Vd=0.01/0.05/0.1 V: converged, Id≈36.15/102.98/121.60 nA; diagnostic samples only.
- Device mode, mesh and contact convergence gates: failed/not established.
- Isolated-state preparation in biased Raised S/D GAA: unresolved.

Read `PHYSICS_VALIDATION.md`, `RUN_INSTRUCTIONS.md`, and `STATUS.json`. `results/final_*` contains exact commands, source snapshots, logs, checkpoints and density/potential/current data. `results/assessment` contains plots and CSV from actual states only.

The MS, NEGF and Spectral registry is retained. Poisson, mapping and the SCF numerical algorithm remain unchanged from STEP09. Two checkpoint metadata fields were added so the existing coordinate-checked restart loader can read new checkpoints directly; see `INTERFACE_CHANGES.md`. The explicit `geometry_profile='uniform'` option permits a diagnostic wire with equal silicon source/channel/drain cross sections; the default remains the original Raised S/D model.
