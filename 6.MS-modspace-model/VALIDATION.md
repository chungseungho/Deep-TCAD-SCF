# Physics and numerical validation

All numbers below come from saved tests or executed solver outputs. `quantitative_validation_pass` is false throughout. The three-level convergence tables are one-bias sensitivities; they do not certify the full geometry/bias domain.

## Scope and independently executed checks

| Check | Evidence |
|---|---|
| Geometry / 4-sided GAA / volumes / shared faces | Per-run geometry_checks.json and mesh/geometry.json |
| P1 FEM, not native FVM | Signed material-weighted tetra stiffness in one DEVSIM potential space |
| Homogeneous and layered capacitor | results/poisson_benchmarks/summary.json; three meshes, analytic potential/normal displacement |
| Uniform charged slab / Gauss / refinement | Same JSON; original source/control volume in actual DEVSIM equations |
| Shared Si/Ox/spacer potential | test_triple_junction.py; exactly one shared unknown, charged real geometry |
| Potential and positive charge exchange | Exact clipped cell/tetra integrals; marginals, affine precision and work tests |
| Full mode basis vs real 3D FD | Independently assembled variable-section real-space matrix in test_physics.py |
| Analytic subbands / six valleys and spin | Unchanged masses and explicit degeneracies; subband_convergence.csv |
| Finite spectral vs direct GF density | results/benchmarks/summary.json |
| Bond vs transmission current | Reduced direct inversion and actual-GAA post_scf_current/summary.json |
| Semi-infinite RGF vs dense inversion | Unequal block-size and complete small open-system density/current tests |
| Pauli / equilibrium / gauge / reversal | Finite matrix tests; no biased standing-state origin guessing |
| Numerical failure retention | Per-bias ERROR/NOT_CONVERGED and complete tracebacks; no interpolated current |
| Overlapping transfer/output biases | results/cross_sweep_consistency.json; independently continued saved states at the same three biases |

## STEP05 audit

Read-only reruns use its original saved density and unmodified source. See step05_audit.json, step05_interface_audit.json and associated log. Its reported grid/mode/length/basis sensitivities match the CSV audit. The material triple-line potential discontinuity was an additional code defect, separate from bath/model limitations. STEP06's shared-DOF solution fixes that defect.

## Fully coupled one-bias comparisons

Vgs=.8 V, Vds=.1 V. Potential residual is unmixed volts; charge residual is volume-weighted relative L2. Mode counts are per representative valley pair, with +/- multiplicity and spin applied separately once. Channel charge includes the suspended bridge and both spacer extensions; the reported barrier window is the gated 24 nm. Neutrality is the outer-half S/D number difference relative to donors, not net device neutrality. A gate-controlled transistor need not have zero total net charge.

| Case | Setting | Id (uA) | Status | Phi residual (V) | Charge residual | S/D neutrality |
|---|---|---:|---|---:|---:|---:|
| dx_0 | 0.5; dx (nm) | 1.32094 | SCF_CONVERGED_PROTOTYPE | 3.97e-05 | 0.000173 | -0.149 |
| dx_1 | 0.25; dx (nm) | 0.476494 | SCF_CONVERGED_PROTOTYPE | 4.34e-05 | 0.000146 | -0.149 |
| gamma_0 | 5; Gamma (meV) | 1.10961 | SCF_CONVERGED_PROTOTYPE | 9.8e-05 | 0.000158 | -0.163 |
| gamma_1 | 10; Gamma (meV) | 1.32094 | SCF_CONVERGED_PROTOTYPE | 3.97e-05 | 0.000173 | -0.149 |
| gamma_2 | 20; Gamma (meV) | 1.24505 | SCF_CONVERGED_PROTOTYPE | 4.25e-05 | 6.53e-05 | -0.121 |
| gamma_3 | 40; Gamma (meV) | 0.748575 | SCF_CONVERGED_PROTOTYPE | 5.7e-05 | 6.5e-05 | -0.0697 |
| basis_0 | fixed 3/6; Basis and channel/S-D modes | 1.32094 | SCF_CONVERGED_PROTOTYPE | 3.97e-05 | 0.000173 | -0.149 |
| basis_1 | local 3/6; Basis and channel/S-D modes | 0.806594 | SCF_CONVERGED_PROTOTYPE | 7.39e-05 | 0.000124 | -0.146 |
| basis_2 | fixed 9/18; Basis and channel/S-D modes | 2.50827 | SCF_CONVERGED_PROTOTYPE | 2.34e-05 | 4.67e-05 | -0.0991 |
| basis_3 | local 9/18; Basis and channel/S-D modes | 2.11956 | SCF_CONVERGED_PROTOTYPE | 6.87e-05 | 0.000291 | -0.0976 |
| length_0 | 12/12; Ls/Ld (nm) | 1.32094 | SCF_CONVERGED_PROTOTYPE | 3.97e-05 | 0.000173 | -0.149 |
| length_1 | 24/24; Ls/Ld (nm) | 1.42604 | SCF_CONVERGED_PROTOTYPE | 8.71e-05 | 5.22e-05 | -0.0981 |
| length_2 | 36/36; Ls/Ld (nm) | 1.49274 | SCF_CONVERGED_PROTOTYPE | 4.99e-05 | 3.08e-05 | -0.0698 |
| mesh_0 | 3.5/1.8/1.3; h/channel/junction (nm) | 1.32094 | SCF_CONVERGED_PROTOTYPE | 3.97e-05 | 0.000173 | -0.149 |
| mesh_1 | 2.5/1.2/0.9; h/channel/junction (nm) | 1.27458 | SCF_CONVERGED_PROTOTYPE | 2.42e-05 | 7.44e-05 | -0.171 |
| mesh_2 | 1.8/0.8/0.6; h/channel/junction (nm) | 1.23515 | SCF_CONVERGED_PROTOTYPE | 5.64e-05 | 0.000202 | -0.158 |
| modes_0 | 3/6; Channel / S-D modes | 1.32094 | SCF_CONVERGED_PROTOTYPE | 3.97e-05 | 0.000173 | -0.149 |
| modes_1 | 6/12; Channel / S-D modes | 2.18576 | SCF_CONVERGED_PROTOTYPE | 6.96e-05 | 0.000206 | -0.122 |
| modes_2 | 9/18; Channel / S-D modes | 2.50827 | SCF_CONVERGED_PROTOTYPE | 2.34e-05 | 4.67e-05 | -0.0991 |
| scf_0 | 0.0001/0.002/0.002; Phi (V) / charge / current limits | 1.32094 | SCF_CONVERGED_PROTOTYPE | 3.97e-05 | 0.000173 | -0.149 |
| scf_1 | 5e-05/0.001/0.001; Phi (V) / charge / current limits | 1.32094 | SCF_CONVERGED_PROTOTYPE | 3.97e-05 | 0.000173 | -0.149 |
| scf_2 | 2.5e-05/0.0005/0.0005; Phi (V) / charge / current limits | 1.32111 | SCF_CONVERGED_PROTOTYPE | 1.98e-05 | 6.85e-05 | -0.149 |
| transverse_0 | 1/1; dy/dz (nm) | 1.32094 | SCF_CONVERGED_PROTOTYPE | 3.97e-05 | 0.000173 | -0.149 |
| transverse_1 | 0.5/0.5; dy/dz (nm) | 1.2421 | SCF_CONVERGED_PROTOTYPE | 8.51e-05 | 0.000425 | -0.15 |
| transverse_2 | 0.25/0.25; dy/dz (nm) | 1.22414 | SCF_CONVERGED_PROTOTYPE | 9.66e-05 | 0.000268 | -0.15 |

Plateau decisions and adjacent current/charge/gated-barrier differences are in each campaign's summary.json. Missing or failed levels fail the plateau test. Fixed/local comparisons are basis comparisons, not monotonic refinement levels. The executed transverse-grid and SCF-tolerance cases appear in the table; incomplete levels appear below. These comparisons remain restricted to one bias, rather than a validated off/threshold/on and drain-bias domain.

## Incomplete campaign cases

These cases have no final accepted point table. A checkpoint or an SCF iteration is not a completed bias. Interrupted computations can be resumed with the same source/configuration; planned but unstarted levels have no physical result. This table records the state when the report was generated.

| Case | Execution status | Last saved iteration | Last phi residual (V) |
|---|---|---:|---:|
| convergence/dx_2 | NO_ACCEPTED_FINAL_STATE | 18 | 0.30816272923994015 |

## Separately executed junction-compatible basis

Low wide-section modes are enriched by zero-extended narrow-section modes and the original full operator is reprojected. The independent extension fingerprint is retained in each state. These results do not replace the baseline sweep or use fitted current/hopping factors.

| Case | Id (uA) | Status | Phi residual (V) | Charge residual |
|---|---:|---|---:|---:|
| nested_spectral_dx0p5 | 4.4914915 | SCF_CONVERGED_PROTOTYPE | 6.582317869134168e-05 | 0.00016079306607103152 |
| nested_spectral_dx0p25 | 4.5128997 | SCF_CONVERGED_PROTOTYPE | 4.5065169298794117e-05 | 0.0001583675321939606 |

The two-grid current difference is 0.474379% relative to the finer grid. This is a diagnostic improvement, with two levels and finite-bath/contact/transverse limitations still present; it is not a certified plateau. The semi-infinite biased attempts are separate: the first attempts have explicit quadrature ERROR records and `nested_negf_refined` retains its unconverged SCF history/checkpoint. None supplies an accepted biased reference.

## Actual-device frozen-field contact comparison

The actual finite-bath SCF field at Vgs=.8 V, Vds=.1 V, dx=.5 nm, modes=3/6 was held fixed. Semi-infinite continuation replaces the terminal x hard walls and distributed finite baths; the interior projected operator and potential are retained. No semi-infinite density or Poisson–quantum fixed point is claimed here.

| Contact / chemical potential | Current (uA) |
|---|---:|
| Finite bath, its calibrated mu=0.0043366198 eV | 1.32093954 |
| Semi-infinite leads, identical mu | 1.22953764 |
| Semi-infinite leads, separately neutral mu=0.0029844975 eV | 1.17282132 |

Three independently initialized adaptive energy panels (.02/.01/.005 eV) are saved in `results/frozen_contact_comparison/summary.json`. This separates a contact substitution from chemical-potential recalibration. It demonstrates a current difference on a fixed field, not agreement of two fully coupled contact models. Sampled transmission plots are distinct from the adaptive current integrals.


## Provisional transfer metrics

| Metric | Actual saved value | Definition / limit |
|---|---:|---|
| Vth | 0.371554 V | Constant current 1e-7 A/um of gated perimeter; 0.034 um |
| Minimum sampled SS | 48.999 mV/dec | Adjacent samples below threshold; not a validated subthermal device claim |
| Maximum gm | 2.13407 uA/V | Numerical derivative of provisional transfer curve |
| Ion / Ioff | 2.83099e+08 | Vgs=1.2 / 0 V, Vds=.1 V |

Definitions and absolute currents are in `results/iv/transfer_metrics.json`. These are properties of the executed numerical configuration, not calibrated experimental predictions. All sweep points retain their individual contact neutrality and dispersion flags.


## Reference limitations and next steps

1. Resolve biased semi-infinite density integration, near-dark states and convergence against independently refined real-energy panels; do not insert artificial eta to obtain a successful plot.
2. Establish the contact model before selecting a teacher: finite constant-Gamma tails and connected neutrality remain important. Increasing transverse basis size changes high-energy tail population and is not automatically a physical contact cure.
3. Abrupt section changes also require a shared retained trace space: results/junction_overlap.json quantifies the missing narrow-to-wide norm and its tx amplification. validation/junction_compatible.py implements and separately fingerprints a nested-subspace comparison; it does not replace the baseline I–V. Complete transverse mesh, mode/basis, Poisson, dx, contact length and energy convergence at off/threshold/on, equilibrium/.1/.5/1 V. Current panels above the actual FD dispersion-valid window remain provisional.
4. Only after a restricted numerical/reference domain is validated should dataset labels become teacher-eligible. No neural model was trained in STEP06.

## Reproduction and cost

All run folders include settings, source/config identity, resource estimate, mesh and raw SCF histories. Runtime per accepted/failed point is in its CSV. Requirements and environment are recorded. Finite dense spectral cost is cubic; local results do not establish speedup relative to RGF. School partition/resource settings are chosen by the user; no remote job was submitted.

The final delivery also includes the exact earlier simulator core in `results/source_snapshots/`, indexed by its original source digest. Final changes to export restoration and exclusion of private virtual environments from the code identity do not rewrite the source identity of existing physical states. The archive is independently re-hashed by `validation/package_audit.py`. Use a separate copy with that archived core to resume the original guarded calculations, or use a new output folder with the delivered code. Field re-exports are checked against the exact saved silicon potential; native DEVSIM coordinates are cm and quantum VTU coordinates are nm.
