# Supplied SIESTA / TranSIESTA code study

Inspected locally on 2026-10-16. Archive SHA256:
`43f05274cd94a7a966a0164e11ad9b4b5576d987925844e00289b580a5b33fc0`.
The archive contains 20,558 members, including a source checkout and build
directories. No supplied binaries were executed, no server was contacted and
no Fortran code is included in the STEP07 distribution.

## Identity of this input

This is **not merely an unmodified standard distribution**. `Src/msdft.F`,
`diagg_msdft.F`, `diagk_msdft.F`, optional Mulliken/Lowdin diagnostics and
generalized finite-difference Poisson modules are present. `diagon.F` branches
on `lmsdft` and calls the modified diagonalizers. This directly establishes an
MS extension in this particular input. It does not establish authorship,
publication provenance or validation of the extra experimental options.
The standard ground-state and TranSIESTA paths remain distinguishable.
`research/siesta_archive_inventory.json` and `siesta_symbol_map.json` preserve
the observed identity. No claim is based on an assumed filename alone.

## Actual paths and numerical counterparts

Paths below are relative to the supplied `siesta-5.4.2/`, not proposed names.
The Python column describes independent EMA analogues, not ports of AO DFT.

| Function | Inspected source path / routine | STEP07 counterpart |
|---|---|---|
| Bare Hamiltonian and overlap | `Src/hsparse.F`: `hsparse` | `modespace/hamiltonian.py`: `ModeBuilder.build` |
| Density-dependent H / Hartree / XC | `Src/dhscf.F`: `dhscf` | `poisson/device.py`: `Poisson.solve`; EMA has no XC |
| DM solver dispatch | `Src/compute_dm.F`: `compute_dm`; diagonalization vs `transiesta` | `solver.py`: backend dispatch |
| Finite eigenstates | `Src/diagon.F`: `diagon`; `diagg_msdft.F`: `diagg_msdft`; `diagk_msdft.F`: `diagk_msdft` | `multispace/occupation.py`: Hermitian eigensolve |
| Ordinary ground-state Fermi occupation | `Src/diagg.F`, `Src/diagk.F` call `Src/fermid.F`: `fermid`, `fermispin` in module `m_fermid` | Equilibrium `fermi` and scalar number solve; no AO DFT claim |
| Finite biased occupation | `Src/msdft.F`: `msdft`, `msfermid`, `msfermispin`, `stepf` | `multispace/constraints.py`, `occupation.py` |
| Projection in nonorthogonal AO basis | `Src/m_msdft_mulliken_partition.F90`: `classify_real_state`, `classify_complex_state`, `partition_region` | `multispace/projection.py`: orthogonal spatial projectors |
| SCF mixing | `Src/m_mixing_scf.F90`; `Src/mixer.F` | `coupling/mixing.py`: volume-weighted Anderson |
| Convergence checks | `Src/scfconvergence_test.F`: `scfconvergence_test` | `solver.py`: physical unmixed residual checks |
| Electrode layout / voltage / μ | `Src/m_ts_global_vars.f90`, `m_ts_options.F90`, `m_ts_voltage.F90`, `m_ts_chem_pot.F90` | `device/regions.py`, `device/contacts.py` |
| Surface GF | `Src/m_ts_electrode.F90`: `SSR_sGreen_DOS`, `SSR_sGreen_NoDOS`, `create_Green`; `m_ts_gf.F90`: `do_Green`, `read_next_GS` | `transport/negf_reference.py`: analytic `surface_sigma` |
| Self-energy construction | `Src/m_ts_elec_se.F90` | `transport/negf_reference.py`: `lead_matrices` |
| Open DM inside SCF | `Src/m_transiesta.F90`: `transiesta`; `m_ts_dm_update.F90`; equilibrium / nonequilibrium contours in `m_ts_contour_eq.f90`, `m_ts_contour_neq.f90` | `transport/negf.py`: contour + resonance-resolved injectivity |
| Open GF and spectral terms | `Src/m_ts_fullg.F90`, `m_ts_full_scat.F90`, `m_ts_trig.F90`, `m_ts_tri_scat.F90` | `transport/negf_reference.py`: block `rgf` |
| Transport postprocessing | `Util/TS/TBtrans/m_tbt_tri_scat.F90`: `A_Gamma`, `orb_current`; `m_tbt_save.F90`: `state_save`, current output; `Docs/tbtrans.tex` | `transport/ms_postprocess.py` |

## Traced finite MS flow

`compute_dm` dispatches to `diagon`; when MS is enabled, `diagon` selects the
Gamma or k-sampled MS diagonalizer. These assemble the actual generalized
Hermitian H/S problem and call `rdiag`/the k solver. The optional partition
path uses the saved overlap before the eigensolver overwrites work arrays.
The occupation routine receives eigenvalues, eigenvectors, k weights,
temperature and total valence charge. The weighted eigenvectors reconstruct
the DM, which returns to the usual density/Hartree/SCF machinery.

The legacy partition in `msfermid` sums |coefficient|² by atom-index region,
then divides by the number of orbitals in each region. This is **not** the
integrated projector probability of the paper on an orthonormal grid. The
optional Mulliken implementation forms Re[c* (S c)], sums by region without
orbital-count averaging, and caches the state classifications. Mulliken
populations are not generally nonnegative probability projectors. STEP07's
orthonormal EMA space permits genuine positive projectors and does not copy
either nonorthogonal approximation.

For occupation, the source initially finds an equilibrium Fermi level by
bisection against `qtot`, then bisects `ef_mod` while maintaining a prescribed
left/right split. With `MSDFT.SD.BiasReference=left`, μL=ef_mod and
μR=ef_mod−VSD after unit conversion; the symmetric option shifts both sides.
The sum of weighted occupations is checked against total charge. This is a
practical number constraint, not independent local reservoir neutrality.
Only states within a finite energy window are reoccupied in this code;
STEP07 evaluates the Fermi tails globally and records their uncertainty.

The supplied code also contains `partitioned`, `half`, `weighted`, eta-half
and fixed-state options. Their existence is not evidence that all are the
published physically validated prescription. STEP07 adopts the supported
origin rule with explicit rejection; it does not transfer half occupation or
fixed arbitrary state fillings. Ties in the supplied partition routine have
deterministic branch fallbacks; this is not a physical resolution of a
degenerate subspace. STEP07 tests the restricted projector spectrum instead.

## Poisson and open transport

`dhscf` constructs density-dependent potentials. In this supplied variant it
also branches to `modified_cg_poisson_solver` when generalized FD Poisson is
enabled, with epsilon, padding and order options under `Src/m_fd/`. Standard
Poisson and TranSIESTA Hartree fixing remain separate. STEP07 retains the
validated STEP06 conformal tetrahedral P1 operator expressed through DEVSIM;
it does not equate that operator with native DEVSIM FVM or with this FFT/FD
code. Potential continuity is enforced by one shared potential DOF.

In the TranSIESTA branch, `transiesta` manages contour integrations, electrode
surface GF, self-energies, scattering functions and DM updates on each SCF
iteration. The lesser density uses both reservoir populations; it is not
ordinary Fermi filling of the finite device eigenstates. Surface GF routines
implement recursive electrode algorithms and support generalized H/S blocks.
TBtrans under `Util/TS/TBtrans/` is a separate postprocessor of saved H/S and
lead information. STEP07 has orthonormal S=I and analytic uniform FD leads,
so its corresponding algorithms are much smaller and have narrower scope.

Official reference checked:
https://docs.siesta-project.org/projects/siesta/en/stable/reference/siesta.html
(TranSIESTA, electrode chemical potentials, contours and TBtrans documentation).
The current online manual is supporting context; version-specific source
claims above come from the actual supplied archive. SIESTA was not compiled
or run for numerical comparison in this work. A successful Python benchmark
is not reported as reproduction of a published SIESTA device result.
