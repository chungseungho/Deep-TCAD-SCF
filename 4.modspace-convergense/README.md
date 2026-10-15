# STEP 04 — Conservative FEM/FVM–FD transfer, lead self-energy, and SCF convergence

**Directory name (requested):** `4.modspace-convergense`. This is a self-contained
copy of the STEP03 driver with numerics upgraded. The frozen `1.example_build`
folder must remain beside it. It imports STEP01 `config.py` and `DevsimPoisson`.
EMA effective masses, Si material parameters, ballistic NEGF and the 3D
DEVSIM Poisson PDE were **not reparameterized**.

## Why STEP03 failed to converge

At Vgs=0.2 V, Vds=0.1 V, `SCF 35/100` through `SCF 100/100` displays a
persistent **two-cycle**:

- Current alternates between about `9.43951e-9 A` and `9.50326e-9 A`.
- `dPhi` stalls at about `4.859e-3 V` (exceeding 0.002 V original tolerance).
- Charge residual stays near `1.05e-2`.
- Inner DEVSIM Poisson solver reaches very small residuals each iteration.

That establishes the *outer fixed-point loop* is not converging, **not**
that the tetra mesh is broken. The Gmsh log says `No ill-shaped tets`; there
is no mesh-quality proof at every point but no reported mesh crash.

Potential causes addressed here: moving energy quadrature nodes with potential,
non-conservative density remap, contact transverse-mode diagonal approximation,
coarse real Poisson mesh, and density-mixing gain. We cannot attribute the
observed two-cycle to a single cause without running the real device again.

## Meshes and coupling

1. Gmsh: original **tetrahedral** 3D Si/SiO2 mesh, smaller target sizes in
   the channel (0.5 nm) and gate-edge junction region (0.35 nm) via physical
   Gmsh Box background fields; outside 1.0 nm. Still DEVSIM 3D Poisson.
2. Quantum: uniform Cartesian mesh (default x=345, interior y=31, z=25).
   Gate channel lies between x=12 nm and x=36 nm; NEGF includes the S/D too
   because the open-boundary leads must attach to both ends.
3. `transfer.py`: read **original physical Silicon tetra connectivity** from
   `mono_nanosheet_nm.msh` (not a Delaunay retriangulation).
   Find containing tetra for each FD site; construct sparse barycentric map P.
   Potential `phi_fd = P @ phi_devsim`.
4. Electron charge: `n_devsim = M_tet^-1 @ P.T @ M_fd @ n_fd`.
   `M_tet` is **geometric nodal lumped tetra volume**, `M_fd` rectangular
   grid quadrature volume. `P.T` is the exact algebraic adjoint to potential
   interpolation under these weighted inner products; it preserves the
   corresponding discrete total charge.

**Important:** This geometric conservation does not yet prove conservation
against DEVSIM's **internal circumcentric/dense control-volume weights**.
For publication-quality results, compare `M_tet` against DEVSIM's
edge/elements-based volume model and, if they differ, adapt the projection
to the exact DEVSIM weights. Also, interpolating finer does **not** add
physics lost by a coarse Poisson discretization; this version refines the
source tetra mesh as well.

### Why NOT a spline between tetrahedra?

Cubic splines across unstructured tetrahedra can overshoot, violate interface
conditions, and break conservation. Piecewise-affine barycentric FE
interpolation reproduces any affine potential **exactly** on original elements.
The target Cartesian mesh can be refined independently by subdividing the
opposite edges/corners, as requested; this is the regular FD discretization.

## Self-energy physical model

For each source/drain, the endpoint **full mode-space onsite matrix** `H_lead`
(including electrostatic off-diagonal mode mixing) is diagonalized:

`H_lead = U diag(lambda) U†`

With longitudinal coupling `-t I`, each eigenchannel has a semi-infinite
1D retarded surface self-energy `sigma_m^r(E)` and

`Sigma^r = U diag(sigma_m^r) U†`.

This is exact **for the assumed homogeneous periodic semi-infinite discrete
lead**, and fixes the STEP03 approximation that discarded the off-diagonal
contact mode coupling. Gamma = `i (Sigma-Sigma†)` must be nonnegative.

We directly compare RGF Transmission **and** local injectivity matrices
`G_i1 GammaS G_i1†`, `G_iN GammaD G_iN†` to a full dense inverted finite
Hamiltonian including the same self-energies. However the validity of an
ideal homogeneous, ballistic contact for doped silicon is still a **physical
assumption** needing contact doping, Fermi level, and band alignment checks.
Bound-state occupation / eta dependence and scattering are not yet validated.

## Numerical changes

- Stable **energy-grid anchoring** to global zero (multiples of `--energy-step`)
  instead of re-centering `linspace` every SCF iteration. Avoids stepwise
  quadrature jitter when potential moves by sub-meV values.
- Density mixing default 0.05, maximum iteration count 250; tighter phi,
  current, charge tolerances. **Reduced oscillation has to be demonstrated in
  full device SCF, not assumed.**
- Gmsh input mesh copied only inside STEP04 `output/`; nothing overwrites
  STEP01 or STEP03. Long jobs save checkpoints per bias.
- Report geometric volume/charge and affine-potential identity before SCF.

## Locally validated (no external software)

```bash
python -m unittest discover -s tests -v
python diagnostics.py --out figures
```

`diagnostics.py` generates **synthetic** tetra/FD and dense-RGF reference
figures; these tests do NOT mean actual 3D DEVSIM SCF is converged.

## SLURM X2: real-device one-bias smoke test

```bash
cd /home/class/2026_SeungHo/TCAD/Deep-TCAD-SCF/4.modspace-convergense
sbatch job.sh --phase idvgs --max-points 2 --out output_smoke
```

Default parameters: `--mesh-h 1.0 --channel-h 0.5 --junction-h 0.35`,
`--nx 345 --ny 31 --nz 25 --modes 6`,
`--max-iterations 250 --mix 0.05 --phi-tol 0.001`
`--charge-tol 0.01 --current-tol 0.005 --energy-step 0.005`.

If that fails to converge, try **energy convergence** at a smaller step and
smaller mixing using a **new output folder**, e.g.

```bash
sbatch job.sh --phase idvgs --max-points 2 --out output_energyfine \
  --energy-step 0.0025 --mix 0.025 --max-iterations 300
```

For production, run the same options without `--max-points`. Existing
`run_settings.json` rejects changing physics/numerics options on resume;
start a **new output directory** for changed settings.

After Id–Vgs converges, choose and scan Id–Vds:

```bash
sbatch job.sh --phase idvds --out output_smoke
```

Note: **do not** use `--phase idvds` until the full 31-point Id–Vgs run
completes; a two-point smoke test intentionally does not extract Vth.

## Reviewable results

- `figures/00_step03_scf_stagnation.png`: analyzed **real STEP03 logfile** 2-cycle
- `figures/01_tetra_to_fd.png`: affine field exact interpolation (synthetic)
- `figures/02_charge_conservation.png`: geometric charge (synthetic)
- `figures/03_selfenergy_dense_reference.png`: RGF vs full inversion (synthetic)
- `figures/04_lead_selfenergy.png`: self energy real/imag and broadening (synthetic)
- `output_smoke/transfer_diagnostics.json`: ACTUAL-device original tets
- `output_smoke/figures/09_transfer_poisson.png`: ACTUAL DEVSIM potential / FD mapping
- `output_smoke/scf_history.csv`: ACTUAL SCF residual history per bias
- `output_smoke/idvgs.csv`, `output_smoke/states/`: ACTUAL IV and state data
- `STATUS.md`: updated by actual calculations as they occur.

Do not interpret a pretty interpolated image as physical mesh convergence:
run full-SCF 3D mesh-size refinement, FD grid, mode-count and energy-step
sweeps to reach final quantitative confidence.
