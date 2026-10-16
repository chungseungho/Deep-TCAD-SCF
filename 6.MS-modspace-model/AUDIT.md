# STEP01–05 audit and STEP06 priorities

The existing stages were read without changing their source or saved results. The pre-work inventory is `results/step01_05_sha256.json`; the final `results/package_audit.json` checks that inventory again. The new simulator is standalone. Only the optional historical audit scripts require the old stage directories.

## Development sequence

| Stage | Existing implementation | Consequence for STEP06 |
|---|---|---|
| 01 | Gmsh/DEVSIM device construction, EMA constants and transverse solver | Preserve the original anisotropic masses and degeneracies; independently include the needed code |
| 02 | Transverse eigenmode, subband and overlap checks, importing STEP01 | A transverse benchmark alone does not verify a driven transistor |
| 03 | Fully coupled GAA NEGF prototype, fixed mode projection and high-bias stress sweeps | Distinguish implemented equations from actually executed end-to-end validation and discrete bandwidth |
| 04 | Original tetra connectivity, weighted adjoint transfer, improved contact self-energy and SCF mixing | Preserve the conservative intent; distinguish P1 stiffness/lumped charge from native circumcentric FVM |
| 05 | Raised S/D, finite-electrode spectral correlation, explicit external reservoirs and saved I–V | Separate physical contact approximations from numerical defects and from closed MS-inspired diagnostics |

STEP03's documented two-cycle and STEP04's proposed remedies are historical observations, not evidence that one unique mechanism explains every later SCF failure.

## STEP05 raw-result audit

`validation/audit.py` reads the original CSV files. Baseline current is **1.6332459617805795e-6 A**. The values below are relative to that saved baseline and reproduce the previously reported sensitivities.

| Change | Actual relative current change |
|---|---:|
| Finer longitudinal grid | −36.8146% |
| More transverse modes | +29.8694% |
| S/D length 24 nm | +24.8270% |
| Local rather than fixed basis | −24.3499% |
| Gamma 5 meV | −21.8366% |
| Gamma 20 meV | +3.6631% |
| Medium Poisson mesh | +1.5289% |
| Finer transverse grid | −4.1424% |

Transfer acceptance was 25/25. Output acceptance was 71/75; the four failed `(Vgs,Vds)` pairs were `(0.8,0.55)`, `(0.8,1.2)`, `(1.2,0.3)`, `(1.2,0.4)` V. Those are historical STEP05 outcomes, not STEP06 results. One small mesh change or one gamma pair does not establish a plateau.

## Additional reproducible code defect

`validation/audit_step05_interfaces.py` loads the original mesh and saved density at Vgs=.4, Vds=.1 V, then solves **the original STEP05 Poisson code**. It does not claim a fresh full STEP05 quantum–Poisson sweep.

| Shared interface | Shared nodes | Maximum potential jump |
|---|---:|---:|
| Silicon / oxide | 419 | 0.0263086663 V |
| Silicon / spacer | 409 | 0.0135351567 V |
| Oxide / spacer | 216 | 0.124627276 V |

The overlapping pairwise interface equations at a three-material junction did not enforce one continuous potential. The previous two-material capacitor tests could miss this topology. STEP06 uses one shared conformal P1 potential unknown, original element permittivities and a silicon-volume source fraction. It enforces potential continuity exactly and displacement continuity in the assembled weak equation. Charged triple-junction tests and independent homogeneous/layered capacitor checks exercise the new formulation.

## Defect versus model limitation

| Issue | Classification | STEP06 action and remaining evidence |
|---|---|---|
| Triple-line potential discontinuity | Code/discretization defect | Shared potential DOF, per-material tetra stiffness, actual Gauss check |
| Independently chosen potential/density maps | Discrete exchange inconsistency | Exact tetra/cell shape-function integrals; positive charge, both volume marginals, affine potential and work adjoint |
| Residual of a mixed update or stale saved field | Acceptance/provenance defect | Unmixed fixed-point checks and a quantum re-evaluation on the exact saved field |
| Coarse FD band and dispersion | Discretization validity | Per-valley local populated energy-window and velocity-error diagnostics; not just comparison of Vds with 4tx |
| Abrupt raised/narrow section at fixed low rank | Joint truncation limitation | Exact overlaps retained; missing trace norm quantified; nested wide/narrow subspace extension separately implemented and tested |
| Finite wide-band bath tails, length and gamma | Contact-model limitation | Same-bath direct GF and independent actual-device transmission checks; actual coupled gamma/length sensitivities; semi-infinite alternative |
| Isolated source neutrality versus connected electrodes | Contact/electrostatic limitation | Both diagnostics saved; no donor/mass/current fitting |
| Closed spatial-weight occupations and nonzero DC | Physical interpretation | Closed MS diagnostics remain distinct from a driven reservoir steady state |

## Remaining priorities

1. Establish an independently converged biased semi-infinite contact reference, including narrow resonances and active weakly connected states. Equilibrium success does not settle this problem.
2. Refine the joint transverse basis, trace compatibility and longitudinal/transverse grids. The two-grid nested-basis improvement is a separate experiment, not a replacement for the baseline sweep or a three-level plateau.
3. Complete coupled mesh/contact/length/basis checks at equilibrium, low/moderate/high drain bias and off/threshold/on gate bias, within the actual discrete energy window. Record missing and failed cases explicitly.
4. Select a restricted validated configuration before using any state as an AI teacher. STEP06 includes no training, surrogate prediction or fitted current correction.

Current acceptance counts, numerical tables and limits are generated from the actual STEP06 outputs in [STATUS.md](STATUS.md) and [VALIDATION.md](VALIDATION.md).
