# DEVSIM-QT 1.0 Release Candidate — model scope

This package is a release candidate, not a certified transistor model or teacher dataset. `RELEASE_CERTIFICATE.json` is authoritative. A numerical SCF status is not a physical validity or accuracy certificate.

## Defined model

- Silicon n-channel raised-S/D nanosheet GAA; inherited dimensions, doping, effective masses and flatband are unchanged in the primary configs.
- Anisotropic six-valley EMA, with spin degeneracy; coherent ballistic RGF/NEGF, 300 K.
- Cell-centred hard-wall transverse electron confinement and coupled slice bases; the nested basis retains actual overlap couplings across section changes.
- One full-transverse-basis, translationally invariant electrode Schrödinger–Poisson reference. Exposed S/D side walls have zero normal electrostatic field. Reference potential has zero spatial mean; charge neutrality determines μ.
- The same reference profile supplies Poisson end-face data and lead eigenchannels. The contact uses `V=-tx U_device† U_lead` and `Σ=V g_surface V†`, including rectangular overlaps. No phenomenological broadening or scattering is added.
- Bias gauge: `phi -> phi+V`, electron Hamiltonian `H -> H-V` eV, `mu -> mu-V` eV. Source reference is `potential_reference_v`; drain adds VD; gate retains VG-flatband.
- Geometry in nm, energy in eV, density in cm^-3, Poisson device coordinates in cm internally, exported quantum/field coordinates in nm, current in A.

## Preserved and changed

The Poisson stiffness, material assembly, charge source, and conservative mapping are preserved. Poisson gained a generic boundary profile setter; native regression verifies unchanged stiffness and actual nonconstant Dirichlet data. SCF now catches a backend-neutral domain exception. It does not identify MS/NEGF/Spectral-specific failures or interpret EMA band diagnostics; model validity is a backend hook.

MS and finite Spectral remain research backends. Their physical limitations have not been solved. `legacy_flat` remains available for regression; primary STEP11 configs use `transverse_sp`.

## Bound density policy and remaining gap

Equilibrium rational-contour density already includes bound spectral weight. A separate bound contribution must not be added again. Independent dense tests verify continuum plus the correctly normalized pole residue against that contour. The residue includes `I-dΣ/dE`; a vector normalized only inside the device misses its lead tail.

At nonzero bias, the production guard rejects detected below-all-lead poles when the reservoir occupations differ materially. It does not invent a preparation history, use Newton iteration order as time, or insert a phonon/probe bath. Left/right contour decompositions agree in tested cases with no isolated states.

**The production audit does not globally exclude symmetry-dark states embedded in another channel continuum or all above-band poles.** A test explicitly demonstrates that limitation. Toy sum-rule, Pauli and current-conservation passes do not close this device-level gap. An actual VG=0.4, VD=0.01 initial-trial sampled pole search also found an above-global-edge Delta_x narrow/dark candidate with source localization 99.884% and Gamma expectation 4.65e-34 eV. This is not an exhaustive audit or proof of exact physical isolation. Unknown occupations and unconverged model choices remain non-teacher data. `rho_continuum+rho_bound` is fully checked only in the independent small references, not certified for the full raised device.

## Numerical limits

- Full transverse lead reference is distinct from the device's reduced basis. Contact basis/DOS convergence, frozen-potential projection at fixed full-reference μ, and Fully-SCF device changes are reported separately.
- The current configs are diagnostic refinement levels, not converged production presets. `select-modes` creates a conservative starting proposal covering occupied bands, the gate energy window, a declared evanescent guard and degenerate multiplets. `refine` must still establish observable convergence.
- Actual nested dimensions are recorded in frozen-mode reports and differ from requested counts.
- Banded resonance guidance retains a closed-system spectral search. It avoids global dense H allocation but is not a proof of every extremely narrow resonance. Constant-memory adaptive integration controls its estimated componentwise error; missed-peak and observable-refinement checks remain required.
- No quantitative off-current, SS, Vth, gm, Ioff floor or full sweep accuracy is certified. Internal current conservation is an algebraic consistency test.
- Unsupported source/drain density decompositions are marked UNSUPPORTED and omitted from exports.

## Contact experiments

A: kernel cut-invariance holds for adding ideal lead layers to the *same* open Hamiltonian; test includes mode rotation. This does not certify the finite device's electrostatic boundary.

B: buffer convergence must keep gate/channel/junction fixed in physical coordinates and change only the computational extension into the same reference lead. A complete biased-device B study is not certified here.

C: changing physical S/D dimensions is a different device. The comparison CLI refuses to relabel such a change as numerical convergence. STEP10 length changes are not counted as B evidence.

## Outside this release

AI training, density predictors, DFT kinetic correction, SIESTA/MS source integration, phonon/probe scattering, drift diffusion, oxide tunnelling, breakdown, impact ionization, self heating and fitted effective masses/flatband. No school-server jobs or GitHub push were made. SIESTA originals are excluded from this archive.
