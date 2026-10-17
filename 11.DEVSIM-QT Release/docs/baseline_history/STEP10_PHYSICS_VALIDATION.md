# STEP10 physical validation: partial success, device certification failed

## Scope and verdict

The open ballistic solver passes the independent benchmark set in `results/independent_final`. Real Gmsh → DEVSIM → conservative mapping → NEGF → Poisson calculations were executed. Raised-source/drain and uniform-cross-section 3D GAA equilibrium solutions converged. A uniform 3D GAA biased point also converged. These are numerical SCF results, not a quantitatively validated GAA model. `DEVICE_SIMULATION_VERIFIED`, `QUANTITATIVE_VALIDATED`, and `TRAINING_READY` remain false.

The original device uses Lg=24 nm, width=12 nm, thickness=5 nm, tox=2 nm, S/D doping=1e19 cm^-3, T=300 K, and the inherited anisotropic silicon masses and valley/spin weights. Flatband=0.5 V is an inherited work-function assumption, not an experimental fit. No mass, doping, broadening or threshold was adjusted to shape an I–V curve.

## Physical audit

- The cell-centered EMA stencil uses the original kinetic constant and anisotropic masses. Missing transverse neighbors implement hard walls; the extra longitudinal half-cell hard-wall term is removed once at each open endpoint. Full real-space stencil projection was independently checked.
- Electron potential energy is `-phi`. Contacts use phiL=reference, phiR=reference+Vd and muR=muL−Vd. The existing common gauge and source/drain reversal tests passed.
- Surface self-energies solve the retarded quadratic for a semi-infinite FD wire. Causality, the quadratic residual, propagation thresholds and wave matching were checked without artificial broadening.
- Density uses the common-mu contour contribution plus real-axis reservoir injectivity. Current and density use the same open Hamiltonian and self-energies. The backend applies `2 * valley.multiplicity` exactly once; volume reconstruction converts nm^-3 to cm^-3.
- The Poisson discretization is the existing conformal P1 FEM stiffness represented through DEVSIM edge/node equations. The integrated potential/charge transfer preserves the discrete exchange work and total charge. Poisson, mapping, SCF numerical iterations and the three existing backend interfaces were preserved. Checkpoints gained coordinate/density aliases needed by the existing restart loader; no iteration equation changed.
- A contact-cut limitation remains: production leads assume a flat reservoir transverse potential while the finite S/D develops a transverse/longitudinal potential profile. No asymptotic contact plateau has been established. Moving the contact cut changes the solution significantly.

## Independent quantum benchmarks actually run

`python validation/step10_independent.py results/independent_final` passed 13 cases.

| Check | Observed maximum error |
|---|---:|
| Uniform/barrier/step transmission vs independent wave matching | 1.33e-14 |
| Dense lesser-GF integrated density, Vd=0…0.1 V | 1.51e-9 per orbital |
| Dense-GF vs RGF bond current | 1.10e-13 A per degeneracy |
| Surface quadratic residual | 2.34e-17 eV² |
| Full 3D stencil vs projected Hamiltonian | 8.89e-16 eV |
| Full 3D open dense-GF vs reconstructed mode density spectrum | 1.61e-14 eV^-1 |
| Uniform multichannel transmission vs open-channel count | 8.89e-15 |
| Transmission outside the lead bands | 0 |

The independent wave-matching solver solves amplitudes in the lattice difference equation and contains no Green function or self-energy calculation. The small 3D reference uses an independently assembled real-space stencil and a dense inverse. Dense inversion is not used as the device production algorithm.

## Bound states and numerical improvement

The old bound audit repeatedly diagonalized a global dense Hamiltonian while finding nonlinear poles. STEP10 replaces that linear algebra with block-LDL inertia and monotone real-axis bisection of `H + Sigma(E) − E`. The Hamiltonian, pole equation, occupation rejection criterion and absence of artificial eta are preserved. New tests compare block inertia with dense eigenvalue counts, and the impurity bound energy with the analytic infinite-chain result and the frozen dense audit.

At Vg=0.4 V, Vd=0.01 V in the Raised S/D device, Poisson images contain isolated poles below both propagating leads with differing reservoir occupations. The solver rejected those images; the diagnostic run reached its 180 s budget without a converged physical final state. Trial logs and checkpoints are retained. This is a computationally bounded failed attempt, not a proof that every possible physical solution is absent.

The audit proves poles only in the retained projected Hamiltonian. It does not establish whether each pole survives full real-space/mode/contact convergence. Embedded dark states and above-band isolated states are explicitly **not excluded**. Near-threshold poles inside the stated exclusion interval are also not certified. No mean chemical potential, added lifetime, omitted energy interval or invented population was used to bypass this issue.

An isolated-state population cannot generally be deduced from reservoir continuum injection alone. Preparation/history dependence is discussed in G. Stefanucci, Physical Review B 75, 195115 (2007), DOI [10.1103/PhysRevB.75.195115](https://doi.org/10.1103/PhysRevB.75.195115). This motivates preserving the failure status, not a claim that the present device poles have already been proven physical.

## Real DEVSIM–NEGF results

The measured points and exact commands are in `results/final_*/run/`. Every final run stores a source snapshot, hash manifest, complete log, exit code and a check that source files were unchanged during execution. Raw NPZ states contain both grids, density, node potential, quantum potential, current, config and units/provenance. VTK/VTU fields are retained for converged states.

Earlier `gaa_*`, `uniform_*` exploratory folders and the first benchmark attempt are retained as explicitly labeled exploratory evidence. They are not selected for final state assessment; `final_*`, `independent_final` and the final native test record are authoritative.

| Device / bias | Numerical result |
|---|---|
| Raised S/D, Vg=0.4 V, Vd=0 | SCF converged; phi residual 1.744e-5 V; density residual 2.173e-5; Gauss error 3.35e-17 |
| Uniform, Vg=0.4 V, Vd=0 | SCF converged; phi residual 7.047e-7 V; density residual 1.254e-6 |
| Uniform, Vg=0.4 V, Vd=0.01 V | SCF converged; Id=3.6153e-8 A; phi residual 8.498e-5 V; density residual 2.984e-4; bond-current defect 5.50e-16 relative |
| Uniform, Vg=0.4 V, Vd=0.05 V | SCF converged after checkpoint continuation; Id≈1.0298e-7 A; exact residuals and current conservation are retained in its state JSON |
| Uniform, Vg=0.4 V, Vd=0.1 V | SCF converged after direct checkpoint restart; Id≈1.2160e-7 A; exact final-state diagnostics are retained in its state JSON |
| Raised S/D, Vg=0.4 V, Vd=0.01 V | Not converged in the bounded run; bound-occupation rejection in physical images |

Additional mesh, contact-length, restart and attempted higher-bias runs are recorded separately. The final `STATUS.json` and `results/assessment` are authoritative for all completed states and failed attempts. An exit-code-124 run is a timeout and never counted as converged.

## Numerical convergence gates

Criteria were recorded in `validation/step10_criteria.json` before the final runs. Density and current require 1% changes across numerical refinements; chemical-potential changes require 0.1 meV. Independent kernel tolerances and SCF tolerances are stated separately.

- **Mode convergence failed.** At the saved Raised S/D potential, increasing (channel, S/D) modes from (8,16) to (12,24) changes spatial density by 18.54% in L2. The source chemical potential also changes by 1.33 meV. Counts and energies of retained-model bound poles change with the basis. This is a frozen-potential diagnostic, not a self-consistent mode plateau.
- **Poisson mesh convergence not established.** One refinement gives a density change of 2.268% and maximum mapped-potential change 17.77 mV.
- **Contact-length convergence failed.** Increasing each S/D from 12 to 24 nm changes density in the matched central domain by 11.93% and potential by 124.13 mV. Contact screening diagnostics remain false.
- **Restart consistency passed within the SCF accuracy.** The recalculated Raised S/D equilibrium state differs by 9.79e-6 in density L2 and 9.89e-7 V in potential.
- A complete separable uniform-lead mode/grid/length diagnostic is in `results/uniform_lead_convergence`. At the full-basis chemical potential, four pure transverse modes miss 45.2% of the 20×15 nm wire's line density, whereas 48 miss about 0.296%. This diagnostic is distinct from the nested enriched GAA basis. Changing the full-basis transverse spacing from 0.5 to 0.25 nm still changes the neutral chemical potential by about 0.29 meV. None of these frozen lead results certifies a device grid/contact plateau.
- Frozen biased-uniform-field energy refinements passed. Tightening panel/order/tolerances/pole count in three stages changed density and current by at most 5.94e-9 and 1.01e-9 relative in the tested final Vd=0.01 V field. These results are in `results/energy_convergence`; they do not replace re-solved device convergence in all other dimensions.

The requested Raised S/D Id–Vgs and Id–Vds sweeps are not certified or fabricated after these failed gates. The uniform diagnostic output CSV/PNG contains only the actual Vd=0, 0.01, 0.05 and 0.1 V converged points at Vg=0.4 V. Missing points, Vth, SS, gm and a physically certified saturation trend are not invented from insufficient data.

## Tests and limitations of the baseline regression

The current native suite reports **74 passed** with no skipped native integration modules. The unmodified STEP09 native baseline first reported 70 passed and one failure in the historical fixed MS NPZ golden (maximum relative density discrepancy about 6.54e-4 under the current BLAS/thread environment). That failed baseline log is retained. STEP10 compares the unchanged STEP08 evaluator and backend in the same runtime with the original strict tolerances; the old golden is retained as historical data. The likely sensitivity of a truncated degenerate MS eigenspace is not reclassified as a physical pass.

The SIESTA tar attachment was not available in the execution workspace and was not read or used. Existing STEP01–09 source was not edited or pushed. No user HPC job was submitted. No AI training was performed, and no result is training-ready.
