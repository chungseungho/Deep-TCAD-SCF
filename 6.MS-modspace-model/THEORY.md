# Theory, discretization and scope

## Electrostatic model and triple junctions

We solve div(epsilon grad phi) = -q(ND - n), with p=0, fully ionized donors and no acceptors. The band reference is Ec=0 at flat, uncharged silicon; electron potential energy is -phi in eV. Source phi=0, drain phi=Vds, gate phi=Vgs-Vfb, Vfb=.5 V. Other exterior surfaces have the natural zero displacement condition. The electron-only approximation is intended for an n-channel ballistic EMA device and is not suitable for ambipolar/inversion-hole, breakdown or recombination physics.

STEP05 used separate DEVSIM material regions with continuous interface equations. A read-only rerun of its original saved device exposed nonzero potential differences at shared material nodes. Pairwise two-material capacitor tests had not detected the three-material junction failure.

STEP06 identifies every conformal geometric vertex with one potential unknown. For each original material tetrahedron,

    K_ij^(t) = epsilon_t V_t grad(N_i) dot grad(N_j)
    V_i^all = sum_t V_t/4
    V_i^Si = sum_(t in Si) V_t/4.

The material-weighted stiffness is expressed as antisymmetric edge flux inside a single DEVSIM numerical region `electrostatic`. Gmsh material volumes and physical surface groups remain separate. Physical interfaces share their nodal potential exactly; normal displacement continuity is enforced weakly by the common assembled equation. There are no overlapping DEVSIM interface equations at triple lines.

This is conforming linear-tetra **P1 FEM with lumped source integration**, equivalently a conservative signed edge balance. It is not DEVSIM's native circumcentric FVM. Obtuse tetrahedra can yield negative edge coefficients; a discrete maximum principle is not claimed. Dirichlet constraints remove the stiffness nullspace. The actual node load is q(ND_i-n_i) V_i^Si, implemented with `SiliconVolumeFraction=V_i^Si/V_i^all`. At material boundaries this factor is essential.

The native model named `DField` stores the algebraic factor epsilon0*(phi_i-phi_j)/edge_length. Its product with `SignedMaterialCouple` is the material-weighted assembled edge contribution. Physical element fields are E_t=-sum_i phi_i grad(N_i), D_t=epsilon_t E_t; displacement-continuity claims refer to the assembled weak equation and analytical capacitor tests. Native edge-factor values alone should be read with this definition.

For Gauss verification, the Dirichlet rows are replaced by the contact equations. We compare the boundary reactions sum(K phi)_B with minus the **interior** free charge, and independently with DEVSIM contact displacement charges. A volume charge at a replaced Dirichlet row belongs to the ideal contact, not an interior equation. The fully assembled residual is also measured on every interior node.

## Exact potential–charge exchange

Let N_i be the P1 potential shape function and Q_j a rectangular quantum cell. We compute

    C_ij = integral_(Si intersect Q_j) N_i(r) dV
    phi_j = sum_i C_ij phi_i / V_j
    n_i   = sum_j C_ij n_j / V_i^Si.

The tetrahedron/cell intersection is a convex polyhedron. Its exact volume and centroid integrate an affine shape function; no fitted smoothing is used. Therefore C>=0, sum_i C_ij=V_j and sum_j C_ij=V_i^Si, up to geometry roundoff. A uniform density, constant potential and affine potential at a rectangular cell's centre are reproduced. Also

    sum_i V_i^Si phi_i n_i = sum_j V_j phi_j n_j.

This establishes the **discrete electrostatic work adjoint** for piecewise constant cell density and cell-averaged P1 potential. It does not prove a nonequilibrium total-energy variational principle. Piecewise constant quantum density, finite Poisson mesh and finite quantum cell sizes still require convergence.

`adjoint`, `balanced_adjoint` and `legacy` mapping options are diagnostic alternatives. The simple point-sampled affine adjoint is positive/conservative but not exactly uniform-density reproducing on an underresolved mesh. Balanced adjoint preserves constants/work but changes affine interpolation. Legacy reproduces STEP05's non-adjoint exchange. Production attempts should use the default `integrated` map.

## Cell-centred EMA and mode space

Geometry and quantum coordinates use nm. KIN=hbar^2/(2m0)=.0380998212 eV nm^2. Each Si valley pair uses the original anisotropic masses; spin 2 and +/- valley multiplicity 2 are applied once each. A cell-centred hard wall is at half a grid interval: a missing neighbour contributes 2t rather than t to the diagonal. This applies to transverse surfaces and raised-S/D shoulders.

Slice modes are zero-extended to the common transverse cell grid. Coupling is

    T_(i,i+1) = -tx U_i^H U_(i+1),  tx=KIN/(mx dx^2).

Different sections are never joined by a guessed identity or by equal mode index alone. The projected longitudinal wall diagonal is retained. In fixed basis, every electrostatic update uses the full U^H diag(-phi) U matrix, including off-diagonal elements. Local basis diagonalizes the current transverse Hamiltonian and updates overlaps. With the full transverse basis, the mode-space matrix is a unitary representation of the actual 3D real-space FD Hamiltonian; a separate test assembles that real-space operator directly.

The discrete longitudinal band is W=4tx, not infinite. For actual local kinetic window K, v_FD/v_cont=sqrt(1-K/W). We record K/W and a 2% velocity-error diagnostic per valley using max(muL,muR)+8kBT and the lowest local subband. This conservative diagnostic is not a proof of continuum convergence. A current decrease near a band cutoff is not evidence of physical NDR.

## Two distinct occupation and contact problems

The closed MS-inspired diagnostic diagonalizes finite H and evaluates L/C/R spatial weights. Exactly degenerate subspaces may be localized by the L-R projector. Ambiguous active-window standing/bound states are refused rather than assigned a guessed origin. Diagonal stationary occupations of a real closed Hamiltonian do not automatically sustain directed DC current. No MS-DFT constrained-search functional or excited-state variational theorem is implemented here.

Default finite spectral transport has explicit external Fermi reservoirs:

    Gamma_L = gamma P_L, Gamma_R = gamma P_R
    H_eff = H - i(Gamma_L+Gamma_R)/2
    C = sum_a integral f_a(E) G^r Gamma_a G^a dE/(2pi).

The SCF algorithm evaluates this correlation through complex eigenstates and analytic digamma Fermi integrals, without energy-resolved G/RGF loops. This still contains a physical bath self-energy. It is not a closed finite system and is not identical to MS-DFT. All retained states and Lorentzian tails contribute. The density matrix is checked for Hermiticity and 0<=C<=I; reservoir-origin parts are positive and additive.

The donor electron number fixes the **isolated finite source** chemical potential using the same grid, modes, masses and broadened Fermi filling. MuR=MuL-Vds. Neither mass nor doping nor current is fitted. Connected electrode neutrality is a separate diagnostic. Uniform wide-band tails, finite electrodes, length, transverse truncation and gamma can affect that neutrality and the current; isolated calibration does not remove those effects. A gamma value producing a plausible I–V is not a contact validation. The strict spacing/gamma/kBT hierarchy flag is a conservative diagnostic inherited for the extended-reservoir study; absence of that hierarchy alone is not a proof of failure for every non-Markovian wide-band model.

The optional semi-infinite NEGF reference replaces the two longitudinal hard walls by analytic nanowire lead continuation. For each lead mode,

    Sigma(z) = (z-a-sqrt((z-a)^2-4t^2))/2

with the retarded analytic branch. Infinite-wire FD k integration calibrates donor line density. Equilibrium density uses continued-fraction Fermi poles; nonequilibrium uses adaptive real-axis lead injectivity. Band-edge sin^2 coordinates remove integrable square-root DOS singularities. This path **uses GF inside SCF**. Its local Pauli checks are necessary and do not prove full-matrix PSD; small-system dense-GF tests provide an independent full-matrix comparison. Unresolved real-energy integration raises an error, rather than adding a fitted eta or changing SCF tolerances. The quadrature estimator does not rigorously exclude missed extremely narrow resonances. Biased dark-state uniqueness and full-device real-axis convergence remain restricted by the actual validation results. The contour-base filling of bound states must not be interpreted as a uniquely established driven occupation for an active disconnected state.

Finite bath and semi-infinite contacts are different physical approximations. Their density/current need not agree at a particular finite gamma/length. Compare identical Hamiltonians/chemical potentials first, then document separately calibrated contacts and fully coupled changes.

## Current, SCF and provenance

For degeneracy g, a positive left-to-right electron drain current is

    I_i = -2 g q/hbar Im Tr(T_(i,i+1) C_(i+1,i)).

The leading 2 is from the commutator and is not another spin degeneracy. Independent transmission current uses g q/h times the energy integral. The finite bath acts in S/D, so current continuity is imposed across the undriven channel and its bounding interfaces, not within a volume containing reservoir injection/extraction.

SCF uses original quantum density and Poisson at every evaluation. Default acceptance needs unmixed phi residual <1e-4 V, volume-weighted charge residual <2e-3 and current change <2e-3. Channel/source–drain conservation additionally requires absolute mismatch <max(1e-15 A, 1e-6*abs(Id)); this absolute floor must be considered for tiny off currents. Gauss relative error must be <1e-7. All actual absolute/relative diagnostics are saved. Screened n/kBT is only a preconditioner. Anderson positivity guards and Newton–Krylov/Armijo do not alter the physical fixed point. The exact saved potential is checked again with its own quantum output density. Missing/failed bias values are never substituted by interpolation or looser tolerances.

Every saved state carries configuration JSON/fingerprint and simulation-source SHA256. The delivered digest excludes generated output, tests, validation drivers and private virtual environments; the final package inventory hashes all delivered files. Earlier core snapshots preserve the exact source identities of existing states when final export/provenance fixes change the delivered code digest. States remain explicitly unvalidated until numerical/reference/contact tests justify a restricted teacher domain.

## Primary sources

- T. Ozaki, PRB 75, 035123 (2007), continued-fraction Fermi contour: https://www.openmx-square.org/tech_notes/tech4-1_1/tech4-1_1.html
- T. Ozaki, K. Nishio, H. Kino, PRB 81, 035116 (2010), contour plus nonequilibrium NEGF: https://arxiv.org/abs/0908.4142
- Lee, Kim & Kim, Advanced Science 7,2001038 (2020), finite multi-space excitation: https://doi.org/10.1002/advs.202001038
- Lee, Yeo & Kim, PNAS 117,10142 (2020), distinct state-origin occupations: https://doi.org/10.1073/pnas.1921273117
- Elenewski, Gruss & Zwolak, JCP 147,151101 (2017), extended reservoirs: https://doi.org/10.1063/1.5000747

These references motivate specific methods. They do not certify this Si GAA implementation or transfer a DFT variational theorem to its Hartree EMA model.

## Abrupt-junction trace-space diagnostic and research alternative

The overlap formula is exact for the chosen retained spaces, but two independently truncated section bases can have almost no common physical trace at a discontinuous shoulder. Increasing tx as dx^-2 then amplifies the mismatch. This is a **nonuniform joint truncation limit**, rather than a license to replace overlap by an identity. `results/junction_overlap.json` records the singular values and missing squared norms. A large `tx*(1-sigma^2)` is a diagnostic scale, not a fitted barrier energy or a proof of one specific I–V error.

`validation/junction_compatible.py` implements a separately fingerprinted alternative. The wide basis is an orthonormal span of its retained low modes and the adjacent narrow basis zero-extended onto the same transverse grid. Every transverse kinetic, full potential, longitudinal shoulder diagonal and overlap matrix is reprojected from the original real-space FD operator. No hopping multiplier, effective mass, donor density or current scale is adjusted. This yields U_narrow^H P_wide U_narrow=I at the step; the full-basis limit is unchanged. Extra S/D ranks and extension SHA/config identity are recorded. Its unit tests establish those identities, while actual SCF outcomes are in the nested_* result folders. Those runs are separate from the baseline sweep and do not certify a complete mode/contact/grid plateau.
