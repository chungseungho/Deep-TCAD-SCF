# Theory boundary and derivation

## Supplied papers, reviewed directly

| Source | Directly supported | Not established by that paper for this code |
|---|---|---|
| Kim & Lee, **A New Direction for First-Principles Device Simulations** (supplied review), sec. 2.2 | finite L/C/R, constrained excited-state density, explicit QFLs; GF transport is a postprocessing step | an EMA variational functional, mode-truncated semiconductor contacts, learned density accuracy |
| Lee, Kim & Kim, **Multi-Space Excitation as an Alternative to the Landauer Picture**, Advanced Science 7,2001038 (2020), DOI [10.1002/advs.202001038](https://doi.org/10.1002/advs.202001038), Eqs. 2-7 | SCF without bulk self-energy replacement; GF transmission and Landauer current AFTER SCF | removal of all GF/embedding from current evaluation; arbitrary finite-state occupation as a unique DC steady state |
| Lee, Yeo & Kim, **Quasi-Fermi level splitting in nanoscale junctions from ab initio**, PNAS 117,10142 (2020), DOI [10.1073/pnas.1921273117](https://doi.org/10.1073/pnas.1921273117), Eqs. 3-4 and Fig. 1 | keep separate L/R-origin occupations (rule I); average chemical potential fails for conducting examples; no strictly channel-localized bias-window state is assumed | a universal origin rule for symmetric standing waves, quasi-degenerate finite EMA spectra or arbitrary bound states |
| Kim et al., **Gate- versus defect-induced voltage drop and negative differential resistance in vertical graphene heterostructures**, npj Computational Materials 8,50 (2022), DOI [10.1038/s41524-022-00731-9](https://doi.org/10.1038/s41524-022-00731-9), Methods | finite graphene electrodes, gate/bias constraints and state-origin tracking in MS-DFT | transferring graphene-specific results, NDR, DFT energetics or screening lengths to this Si GAA model |

MS-DFT total-energy minimization and its variational theorem have NOT been
implemented or proved for this Hartree-EMA model. Naming it an implementation
of MS-DFT would therefore be inaccurate.

## Geometry and discretization

Transport is x. Raised source: `[0,Ls]`, gate begins at `Ls+Lsp`, ends after
`Lg`, raised drain begins after the second spacer. A continuous silicon bridge
of width W and thickness T joins the raised blocks. Its four faces are covered
by oxide under the gate and SiN under the spacers. The bottom gate is present.
Uncovered Si surfaces are hard walls for the electron EMA; holes and oxide
wavefunction penetration/tunnelling are omitted. Electrostatic exterior surfaces
other than contacts have natural zero normal displacement flux.

Quantum space uses cell centres, with `Vcell=dx dy dz`. Si surfaces coincide
with cell faces. Half-cell Dirichlet hard walls are imposed using doubled
boundary diagonal contributions. Transverse eigenvectors have Euclidean cell
normalization; density divides by Vcell. Spin 2 and +/- valley multiplicity 2
are included once each, using the unchanged masses from STEP01.

For a common full transverse cell grid and zero-extended slice modes U_i,

`H_i,i+1 = -tx U_i^dagger U_i+1`, `tx = hbar^2/(2 mx dx^2)`.

The missing-neighbour longitudinal face diagonal is also projected. Thus
the coupling is a Galerkin projection of the actual variable-cross-section
real-space operator. Retaining every transverse degree of freedom recovers
that real-space matrix exactly (tested). Truncation to a small number of modes
is a controlled numerical approximation only after convergence testing.

## Two occupation routes

**Closed MS-inspired route.** Diagonalize the finite Hermitian H, compute
`w_L/C/R = <psi|P_L/C/R|psi>`, and assign distinct Fermi occupations as in
the supplied rule I. Rotate only machine-degenerate subspaces with the
projector `P_L-P_R`. A state with indistinguishable L/R origin in the active
window is reported as unresolved and rejected. No averaging of mu or imposed
travelling-wave conversion fixes this ambiguity. For real nondegenerate
eigenstates, a diagonal occupation matrix commutes with H and cannot by itself
carry a directed DC bond current. This backend is an occupation diagnostic.

**Default finite spectral route (new algorithmic choice here).** The finite
L/R volumes are coupled to external wide-band Fermi environments:

`Gamma_L = gamma P_L`, `Gamma_R = gamma P_R`,
`H_eff = H - i(Gamma_L+Gamma_R)/2`.

This is a phenomenological **bath embedding self-energy** `Sigma_env=-i Gamma/2`.
There is no semi-infinite spatial lead or separate bulk-Hamiltonian replacement,
but an external source of dissipation and electrons is present. A true DC
steady state would not otherwise exist in the closed finite system.

The model fixes occupations through the Fermi functions at muL and muR,
with `muR=muL-Vds` (energies in eV, biases in V). It includes reservoir-induced
coherences and avoids arbitrary finite-state origin labels. It is neither
the MS constrained-search principle nor a proof that MS and EMA are equivalent.

## Analytic density without energy-by-energy Green matrices

The finite-bath correlation is mathematically

`C = sum_alpha int f_alpha(E) G(E) Gamma_alpha G(E)^dagger dE/(2pi)`.

This identity defines the reference, but SCF does not evaluate G(E).
Diagonalize `H_eff = V diag(z_m) V^-1` and write
`B_alpha = V^-1 Gamma_alpha V^-dagger`.

Let `p_m(mu)=digamma(1/2+i(z_m-mu)/(2pi kBT))`. Then

`K_mn(mu) = [p_m(mu)-conj(p_n(mu))-i pi] / [2pi (z_m-conj(z_n))]`,

and `C = V [sum_alpha B_alpha .* K(mu_alpha)] V^dagger`.
The single-level diagonal reduces to

`n = 1/2 - Im digamma[1/2+(Gamma/2+i(epsilon-mu))/(2pi kBT)]/pi`.

The code evaluates the equilibrium part as the matrix function

`Ceq(muR) = I/2 - [F-F^dagger]/(2pi i)`,
`F=V diag(p(muR)) V^-1`,

and adds `V [B_L .* (K(muL)-K(muR))] V^dagger`. This avoids solving an
ill-conditioned population equation for almost isolated high-barrier states.
Tiny widths are obtained from the positive dissipativity identity
`-2 Im(z_m)=v_m^dagger Gamma v_m/(v_m^dagger v_m)` instead of subtracting
nearly real eigenvalues. Disconnected states with unresolved bias-window
occupation are rejected. Pauli bounds and eigensystem residuals are recorded.

All finite retained states contribute, including Lorentzian tails; there is
no omitted occupied-state cutoff. The retained transverse Hilbert space still
requires convergence. Dense complex eigendecomposition is cubic, so this is
not automatically faster than RGF; timing must be measured before AI claims.

## Current and independent check

Using the same C, electron particle flux converted to positive drain-current
convention across slice i is

`I_i = -2 g e/hbar Im Tr(H_i,i+1 C_i+1,i)`.

Here g is spin times valley multiplicity; the leading 2 comes from the
commutator, not another spin factor. Reservoir source terms are absent within
C, so all undriven channel bonds must agree. In equilibrium H and Gamma are
real/time-reversal symmetric and C is real, giving zero current exactly.

Postprocessing may explicitly evaluate finite `G(E)`,
`T(E)=Tr[Gamma_L G Gamma_R G^dagger]`, and the Landauer integral. The
small-system tests compare this independent NumPy inversion/quadrature to the
analytic density and bond current. Sampled full-device T(E) is separately
labelled; narrow resonances may require much finer integration.

The optional Lindblad route instead uses f(e_k) in each isolated reservoir
eigenmode and solves `A C+C A^dagger=-Pi`. It is NOT identical to f(E)
finite-bath filling. Its current integrand is not in general just
`T(E)[fL(E)-fR(E)]`. See Elenewski, Gruss & Zwolak, JCP 147,151101 (2017),
DOI [10.1063/1.5000747](https://doi.org/10.1063/1.5000747), and Gruss et al.,
JCP 147,141102 (2017), DOI [10.1063/1.4997022](https://doi.org/10.1063/1.4997022).

## Remaining physical limits

Finite silicon reservoirs are strongly confined compared with bulk doped Si.
The common source mu is calibrated to their isolated donor electron count
using the configured filling, masses and modes; it is not fitted to measured I.
Connected outer-electrode neutrality is an independent diagnostic and can still
fail. Fixed source/drain Dirichlet voltages and a constant finite bath width are
phenomenological contacts. Lead length, gamma and discrete-level spacing must
be jointly examined; a gamma plateau at one mesh alone can be misleading.

Poisson mesh, quantum spacings, mode count, basis choice, mapping smoothing,
contact length/broadening, and the inherited 0.5 V flatband assumption can all
affect I-V. No phonons, traps, recombination, holes, nonparabolicity, exchange-
correlation, Si/oxide wavefunction leakage or calibrated metal workfunction
are included. There is no proof of a nonequilibrium variational total energy.

An additional discretization limitation is the finite longitudinal FD band:
`W_v=4 KIN/(mx dx^2)`. At dx=1 nm the Delta_x band is only about 0.166 eV
wide. Therefore falling output current at Vds approaching 1.2 V on that coarse
grid cannot be interpreted as a physical Si NDR or validated saturation curve.
The default dx=.25 nm improves this bound to approximately 2.66 eV; multiple
successive longitudinal refinements are still needed.
