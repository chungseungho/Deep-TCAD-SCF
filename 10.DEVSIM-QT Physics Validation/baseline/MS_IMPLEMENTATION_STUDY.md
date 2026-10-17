# Practical MS implementation study

Research completed before STEP07 solver implementation, 2026-10-16. MS means
**Multi-Space**. Transverse mode space is an independent Galerkin reduction.

## Evidence and access limits

1. Lee, Kim, Kim, *Advanced Science* 7, 2001038 (2020),
   https://doi.org/10.1002/advs.202001038. The user supplied the six-page paper;
   equations 2–4, 8–10, Figure 1 and Experimental Section were inspected.
2. Lee, Yeo, Kim, *PNAS* 117, 10142–10148 (2020),
   https://doi.org/10.1073/pnas.1921273117. User-supplied seven-page paper,
   Results: Electron Occupation Rule, equations 1–4 and Materials and Methods.
3. User-supplied review *A New Direction for First-Principles Device Simulations*
   and 2022 vertical graphene device paper are supporting context, not validation
   of a silicon EMA model.
4. Actual supplied `siesta-5.4.2.tar.gz` was inspected locally. Its inventory,
   archive SHA and selected source SHA values are in `research/`. It contains
   additional MS modules. Their presence does not authenticate their authorship
   or prove identity with the 2020 implementation. No source is redistributed.
5. Publisher SI links were located and attempted through the publisher and PMC.
   Wiley returned 403; PNAS/PMC returned cookie/challenge HTML instead of the
   PDFs. **SI contents were not obtained and are not claimed as inspected.**
   Unpublished implementation details cannot be reconstructed from the main
   papers with certainty. This limit is part of the study, not filled by guesses.

## What the papers actually do

The finite model contains a sufficiently large electrode slab at either end.
Only a bulk-like outer part is a reservoir. Interfacial metal layers that carry
voltage drop belong to the extended channel. In the BDT example there are eight
Au layers in each electrode and the outer three are reservoirs; this is a
molecular example, not a recipe for semiconductor reservoir lengths. The ground
state supplies the initial KS Hamiltonian and density. A finite-bias solve then
updates the Hartree and exchange-correlation contribution to that Hamiltonian.
The excitation picture moves electronic occupation between spatial spaces under
the constraint μL−μR=eVb. The two reservoirs have Fermi distributions at their
respective chemical potentials. Electrostatic potential and electrochemical
potential are distinct quantities.

For a normalized KS state, equation 8 compares its integrated probability in L,
C and R; the largest region defines its spatial class. For a C-dominant state,
equation 10 compares its remaining L and R weights to select the corresponding
reservoir Fermi distribution. PNAS rule I retains these two occupations; its
rule II instead uses a mean chemical potential and fails for conducting
junction examples. That result does not justify a universal theorem that every
finite real eigenstate has an injection direction. The papers explicitly assume
no strictly channel-localized state in the bias window. They do not supply a
general, basis-independent prescription for an arbitrary degenerate standing
wave subspace. A zero or near-zero L/R weight difference is therefore a failure
of origin evidence, not evidence for either reservoir.

The density is the sum of occupied state projectors, with temperature, spin,
k-point weights and electron count treated consistently. The supplied modified
source makes a practical extra detail visible: `msfermid` first brackets an
equilibrium Fermi level, partitions states, then bisects a common μ offset while
keeping the prescribed μ difference so that the sum of occupations equals
`qtot`. `WKE` is then used by the diagonalization routines to reconstruct the
density matrix. This enforces a **total finite-system electron count**, not a
separate equation proving local neutrality of each reservoir. Local neutrality,
screening and electrode-size convergence still require measurement. The source
also has legacy per-orbital averaged coefficient weights and an optional
Mulliken path. Neither is silently equated with an orthogonal real-space
projector. STEP07 uses normalized spatial projectors on its orthonormal grid.

During SCF, the density/DM and effective KS potential change, eigenstates are
recomputed and reclassified, the Fermi offset is solved again, and density or
Hamiltonian mixing stabilizes the iteration. The original constrained-search
functional additionally refers to the ground-state density and orthogonality
to lower excited states. The public equations establish that DFT formulation;
the observed occupation loops alone are not an executable minimizer over all
many-electron constrained wavefunctions. In particular, the provided source
does not turn an arbitrary EMA fixed-point iteration into such a minimizer.

After electronic SCF, the papers attach scattering boundaries using surface
Green functions from separate electrode calculations. Matrix GF transmission
and Landauer current are then postprocessing. Their quantitative comparisons
use the same L/C/R partition and electrode surface GF as the NEGF calculation.
TranSIESTA instead constructs the open-system correlation/DM using energy
contours and reservoir injectivities **inside SCF**. TBtrans postprocesses a
saved Hamiltonian; it does not supply a missing electronic SCF by itself.

## STEP07 numerical contract

The implementation is **Multi-Space-Inspired EMA/Hartree**, not MS-DFT. It has
anisotropic silicon masses, explicit donor charge and a shared-potential 3D
Poisson discretization. It has no XC functional, no constrained excited-state
universal functional, and no assertion of variational energy convergence.

1. Define five physical regions from a device longitudinal coordinate: L
   reservoir, L extension, gated/scattering channel, R extension, R reservoir.
   Only the outer configured lengths enter reservoir projectors. Merge the
   extensions with C for occupation. Measure potential spread/gradient,
   populations and donor neutrality separately; sweep boundary locations.
2. Solve the finite, closed, Hermitian Hamiltonian. Verify orthonormality and
   that wL+wC+wR=1. Record largest-region labels and electrode-origin labels
   separately. Use the paper's channel-origin comparison only where its
   evidence is stable. Require finite reservoir reach, nonzero L/R contrast,
   stable sign when moving the boundary one slice, and a stable near-degenerate
   subspace. The thresholds are rejection/sensitivity controls, not current fits.
3. Never rotate eigenstates to manufacture origins. Within a degenerate cluster,
   examine the spectrum of the restricted PL−PR operator: indefinite or small
   eigenvalues imply basis-dependent origin. A mixed-origin cluster is rejected
   if reservoir occupations differ materially. Equal equilibrium occupation
   gives a basis-invariant projector sum and is safe without assigning a lead.
4. At finite bias solve μL=μ0, μR=μ0−VDS (eV units), optionally solving μ0 for
   the donor-integrated finite electron count. A separate fixed-contact-μ mode
   exists for ensemble comparisons. Do not silently force both neutrality and
   two independent μ values with one unknown. Gauge shifts apply to H, μ and
   contact band edges together. No averaged μ is used.
5. Ambiguous states with |fL−fR| above the declared tolerance stop the electronic
   solve. Outside the occupied/active window, using fL is permitted only with a
   reported worst-case total electron-number uncertainty bounded by the sum
   of |fL−fR| times degeneracy. This is a controlled tail approximation, not a
   new rule for bound states. Biased dark states remain undetermined.
6. Form ρ=U diag(f) U†. Require 0≤f≤1, PSD, electron-number equality, positive
   real-space density, and zero closed stationary bond current. Construct L/R
   density contributions only for states with a supported origin; retain an
   unresolved/equal-occupation part explicitly.
7. Couple this density alone to Poisson. Test the unmixed density and potential
   fixed-point residuals on the exact saved H/field. GF is absent from A's
   electronic SCF. Evaluate transport once it converges, with analytically
   continued uniform semi-infinite FD leads and RGF. Current-carrying open
   density is a separate object from the finite stationary density.
8. Backend B performs independent open NEGF SCF; C preserves the finite
   wide-band spectral contact model. Compare A density with nNEGF[φA] at the
   **same** μ and boundary bands, and also compare independently coupled fixed
   points. If finite-N and infinite-reservoir ensembles differ, report that fact.

## Validity and rejection cases

Large screened reservoirs, a constriction, a resolvable origin and converged
electrode/basis/mesh limits are prerequisites to testing this approximation.
Finite real standing-wave occupations do not automatically equal scattering
state injectivities. A transmitting symmetric wire is an essential
counterexample to naive weight classification. Resonant or weakly coupled
states need particular care; a truly disconnected bias-window state has no
unique steady occupation without a preparation/relaxation model. Neither an
arbitrary half occupation nor artificial broadening fixes that missing physics.

A closed finite model can meet Poisson residuals while giving incorrect open
density. Agreement of postprocessed current alone is insufficient. Density
norms, regional charges, barrier shifts, reservoir diagnostics, independent
fixed points and three-or-more-level mesh/mode/contact studies determine the
actual validation level. STEP06's two nested-grid values are only motivation.

High voltage adds separate EMA and dielectric limitations. Finite-difference
bandwidth/velocity errors and actual element electric fields are evaluated;
the model has no tunneling through the oxide, dielectric breakdown, impact
ionization, high-field scattering, heating or nonparabolic bands. Numerical
completion and physical validity are separate status fields.
