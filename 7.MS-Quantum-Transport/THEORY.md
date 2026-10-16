# STEP07 theory and numerical contract

## Model and units

This program implements an electron-only silicon effective-mass/Hartree model.
It is **not variational MS-DFT**. The original MS constrained-search functional,
exchange-correlation energy, atomistic orbital basis and total-energy
minimization are absent. See `MS_IMPLEMENTATION_STUDY.md` for the actual
published occupation algorithm and the limits of transferring it to EMA.

Geometry is in nm, H and μ in eV, Poisson coordinates in cm, potential in V,
density in cm^-3 and current in A. A positive electrostatic potential lowers
the electron band energy: Ec=−φ relative to the chosen band zero. Positive VDS
sets μL−μR=VDS in eV. Spin is a factor of 2; each listed silicon valley has a
separate ± multiplicity of 2. These factors appear once in both density and
current. The canonical longitudinal/transverse frame is mapped to a world-axis
permutation; the block transport algorithms do not depend on a Cartesian axis.

## Electrostatics and exchange

Gmsh constructs actual raised source/drain, rectangular silicon sheet, oxide
wrap and spacers with conformal shared faces. DEVSIM solves one shared
potential DOF per physical mesh vertex. Material-specific P1 tetrahedral
stiffness is assembled before expressing the operator as signed edge balance.
This is **lumped-source P1 FEM through DEVSIM**, not a claim of native
circumcentric FVM. Negative individual edge couplings do not imply negative
element energy; a discrete maximum principle is not asserted.

The default `integrated` exchange matrix is

\[
C_{ij}=\int_{\text{Si}\cap\text{cell }j}N_i(r)\,dr,\quad
P=\mathrm{diag}(V_j^{-1})C^T,\quad
n_i=\sum_j C_{ij}n_j/V_i^{\rm Si}.
\]

Exact tetrahedron/rectangular-cell intersections preserve constants, total
charge and potential-density exchange work to floating precision. Other
legacy mapping options are retained as explicit comparison options and are
not the validated default. Dirichlet node volume loads are excluded from the
interior Gauss balance, with boundary reactions included independently.
Physical element electric field is −grad φ; displacement additionally includes
that element's material permittivity.

## EMA and reduced basis

The cell-centered anisotropic FD kinetic operator uses silicon transverse hard
walls and half-cell Dirichlet distances. Three pairs of masses are inherited
unchanged: (.19,.19,.916), (.916,.19,.19), (.19,.916,.19) m0. All slices share a
common physical transverse grid. Missing silicon is zero extension. Onsite
potential, shoulder hard walls and interslice hopping are all projected from
the same real-space operator; full offdiagonal mode coupling is retained.

The nested option enriches wide source/drain subspaces with zero-extended
narrow-sheet modes. It preserves the narrow trace at the abrupt junction
without rescaling hopping. It does not prove completeness. Increasing retained
modes can change the low-energy DOS and eigenstates, not just add a negligible
unoccupied tail. Fixed μ and fixed finite electron number are different
convergence experiments; report the μ change when the latter is used.

For each valley, the longitudinal lattice bandwidth is 4t, t=KIN/(mx dx²).
At kinetic energy Ekin above the local retained subband edge, the lattice to
parabolic velocity ratio is sqrt(1−Ekin/(4t)). `analysis/validity.py` evaluates
that error on the actual populated energy window, rather than checking only
VDS. Having Ekin below 4t is weaker than a small dispersion error.

## A: finite MS-inspired density

The five-region geometry is merged to L/(extensions+C)/R for projectors.
Eigenvectors remain those of the unmodified Hermitian finite H. State weights
are positive, sum to one, and largest-region class is recorded separately from
electrode origin. The published L/R channel rule is applied only with adequate
reservoir reach, contrast and boundary/subspace stability. No rotation creates
an artificial direction. A mixed-sign restricted PL−PR spectrum in a nearly
degenerate cluster rejects an active origin assignment.

The default finite-number mode solves a common μ offset at every SCF step so
Σν gi fi equals the integrated donor electron count. A fixed-contact-μ option
isolates the alternative reservoir ensemble. Neither condition implies local
neutrality of both reservoirs. An unresolved active/bound state stops the
calculation. Inactive/tail uncertainty is explicitly bounded. The density
matrix is U diag(f) U†, with Pauli bounds and electron count checked before
reconstruction on the 3D quantum grid. A real stationary closed density carries
zero bond current. It is not the open current-carrying scattering density.

GF/RGF is called only after A's electronic SCF converges. Terminal half-cell
longitudinal hard walls are removed and homogeneous semi-infinite FD leads
are attached. Their projected transverse Hamiltonian and hopping are identical
to those used by the open reference. T=Tr ΓL Gr ΓR Ga and I=(q/h)Σg∫T(fL−fR)dE.
Equilibrium gives I=0; that alone does not validate biased transport.

## B and C

B computes semi-infinite NEGF density **inside SCF**. A continued-fraction
Fermi contour gives the equilibrium base; the bias correction integrates the
high-μ reservoir injectivity. RGF diagonal blocks reconstruct density, while
injected bond flux and terminal transmission provide independent conservation
checks. The rational Fermi approximation is checked over a conservative H
energy interval and its pole count is increased if necessary.

Real-axis integration uses lead band-edge transformations and estimated
resonance positions/endpoint widths to seed panels. A componentwise q/2q
adaptive error budget accepts or rejects the integral. No artificial eta,
lifetime or broadened SCF tolerance is introduced. Panel-seed estimates are not
rigorous guarantees against every missed pole; independent panel/order/tail
changes remain necessary. A nonlinear bound-pole audit below all lead bands
rejects states whose two reservoir fillings differ materially. The absence of
embedded dark states has **not** been generally proved. Necessary local Pauli
tests are not a full-density-matrix PSD proof.

C retains STEP06's explicit wide-band baths distributed over finite source and
drain blocks. Analytic complex-eigenstate/digamma integration gives its density
and current. Its broadening Γ is a real contact-model parameter; it is neither
the numerical eta of B nor part of A's closed MS density. Γ, reservoir length
and basis convergence must be measured before comparing it with an infinite
lead model.

## Nonlinear solve and acceptance

The raw density residual is a volume-weighted relative L2 norm. The raw
potential residual is max|Poisson[n_quantum]−φ|, evaluated without damping or
screening. Anderson mixing and positive approximate susceptibility are only
preconditioners. On sustained residual, Newton-Krylov solves the original
potential map. Its screened inverse uses the same assembled Poisson matrix
and a local quantum-capacitance approximation for the linear solve; it does
not alter physical source charge or residual tolerances.

Final density is recomputed on the saved potential; the Poisson field is
restored from its saved input charge before export. Source/config SHA values,
all settings, raw residuals, populations, neutrality, conservation and numerical
versus physical statuses accompany every accepted state. Failure CSV/NPZ
entries remain failures. A small numerical residual is not physical validation.

## Scope

This model omits oxide tunneling, breakdown, impact ionization, high-field
scattering, heating and nonparabolic bands. The configured 0.5 V/nm field
threshold is a diagnostic flag, not a breakdown law. Actual tetrahedron fields
and dispersion determine warnings. A failed-iterate field is labeled as such;
it is not a converged high-voltage prediction. Geometric/contact interfaces
support extension work, but other device families and a physical vertical
transistor have not been validated.
