# Executed validation

## Test completion

56 tests pass: 52 quantum/numerical and four actual native
DEVSIM/Gmsh tests in isolated processes, each with a completed JUnit XML.
`results/test_report/summary.json` is the completion authority. An earlier mixed
run exited before its summary; its zero exit code was not counted as a pass.
Native logs retain a NumPy long-double warning; model arrays use float64.

Coverage includes positive projector sums, electron count, C-dominant origin,
weak/bound-state rejection, degenerate rotations, equilibrium/gauge limits,
Pauli bounds, masses, full 3D projection, nested trace, independent density/GF,
RGF/dense agreement, Poisson/Gauss/interface/capacitor checks and the screened
Newton inverse. The native/quantum test isolation is also the recommended
portable test command, not a suppressed failure.

## Independent small models

Chain, barrier and well checks pass (`results/analytical/`). Uniform T error is
3.553e-15; equilibrium population error
2.220e-16; current relative error
4.926e-08. The discrete well spectrum error is
3.331e-16 eV.

`results/reduced_realspace/` executes three grids (1,.5,.25 nm), with full,
fixed-low, nested-low and local-nested-low bases, against an independent 3D
stencil. Full-basis density and projection identities pass. This is not a
continuum or self-consistent GAA convergence plateau.

## Actual GAA and conservation

| Case | Recorded result | Scope |
|---|---|---|
| MS finite-N .8/0 V | phiR 6.218e-08 V; nR 8.826e-08; Gauss 4.399e-17 | Equilibrium only |
| MS transfer | 25/25 origins rejected | No current curve |
| MS output | 3 equilibrium converged, 72 biased rejected | No biased current curve |
| Finite spectral .8/.1 V | I 4.4922834576e-06 A; conservation 2.168e-14 | Separate finite-bath model |
| Independent finite-bath T integral | max bond mismatch 8.175e-09 | Same frozen H, three panel settings |
| Repaired semi-infinite integral | density difference 1.011e-13; I difference 3.100e-15 | Two panels, no eta; not coupled SCF |

All accepted-state Gauss checks and raw unmixed residuals are recorded per
state. One potential DOF spans material interfaces; the operator is the
configured tetrahedral P1 discretization, not native circumcentric FVM.

## Density and ensemble comparisons

At the finite-N MS field and mu: N_MS=72.6,
N_NEGF=74.97370028; relative volume density L2 error
6.06523%. That mu differs from the independently
neutral infinite-lead mu. These are frozen-field results, not same-contact
coupled equivalence.

At the **same neutral-lead mu**, separately coupled MS and NEGF fixed points
differ by 24.06567% in density and
9.52178 mV in maximum nodal potential.
The frozen open-density error on this MS field is
24.78782%. Regional populations and all-valley
lowest-subband channel barriers are in `matched_mu_coupled_comparison/`.
Equilibrium I=0 does not establish density equivalence.

Finite-N MS reservoir deficits at Vg=.8 are
-47.546% (L) and
-48.081% (R).
Screening targets fail. Boundary placement is a sensitivity test, not an IV fit.

## Numerical limits

Actual frozen GAA subspaces 2/4, 4/8 and 6/12, plus plain/local comparisons,
record dimensions, orthogonality, junction missing norm, density and population
above mu+8 kBT. Large added-mode density is not called a negligible high-energy
tail. No mode plateau is certified.

The donor field has different discrete measures: quantum cells contain
72.6 donor electrons, while Poisson nodal lumping gives
72.8409188782, a 0.33184% quadrature discrepancy.
The finite-N constraint uses the declared quantum-cell count; Poisson Gauss
uses its actual nodal load. This error is reported, not hidden by donor fitting
or a claim of exact global neutrality.

`results/numerical_plan/` contains 22 **unexecuted** configurations for three
levels of longitudinal/transverse/mode/Poisson/reservoir/contact studies and
fixed/local/plain/nested comparisons. A plan is not a result. A joint coupled
GAA current/density plateau and a bias-domain contact limit remain unfinished.
STEP06's two-grid 0.47% difference is not treated as a plateau.

## Stress, gauge and restart

Six gate points through 6 V are all origin-rejected. Six drain points through
6 V give one equilibrium solve and five rejections. Numerical status and
available-iterate field/scope flags are separate. Fields during a rejected
Newton trial can be the last available Poisson field, not the unevaluated
trial Hamiltonian; they are not SCF predictions.

A .4 V common potential shift is tested in actual GAA SCF and in frozen open
density/current/transmission. The finite-tolerance SCF differences are in
`gauge_scf_comparison.json`; exact frozen gauge tests are separate. Resume
preserves accepted NPZ bytes; regenerated CSV physical values match that state.

## Not executed or not established

- No accepted biased MS GAA fixed point or reproduced MS transport IV.
- No certified biased GAA NEGF fixed point under the bound-preparation contract.
- No constrained-search functional or variational energy proof.
- No full GAA joint mesh/mode/contact convergence certificate.
- No same-geometry DD solve: supplied U-channel topology, dimensions, channel
  acceptor doping and oxide differ (see DD_COMPARISON).
- No SIESTA binary compilation/run or atomistic published-device reproduction.
- No validated vertical/planar device family or missing high-field physics.

Every quantitative_validation_pass flag remains false. Failures, source/config
IDs, source snapshots and state residuals remain available for review.
