# STEP 02 — 2D Transverse Schrödinger + Mode-Space basis

**Scope**: Verify the existing EMA transverse solver mathematically, inspect
mode energies and wavefunctions, and prepare the mode-space Hamiltonian basis.
No DD, no AI, no new EMA effective masses, no claim of a converged NEGF current.

## Directory assumption

```text
Deep-TCAD-SCF/
  1.example_build/        # previous working version; untouched
    config.py
    negf/mode_space.py
    build/result/poisson_only.npz   # if retained (optional)
  2.modespace/            # this step
    run.py
    tests/test_transverse.py
    results/              # auto-generated and gitignored
```

`2.modespace` imports `Config`, `HB2_OVER_2M0`,
`transverse_hamiltonian` and `build_mode_slices` directly from
`1.example_build`. The old EMA equations and mass parameters are not edited.

## Run

Activate the existing `.venv`, then:

```bash
cd /home/class/2026_SeungHo/TCAD/Deep-TCAD-SCF/2.modespace
python -m unittest discover -s tests -v
python run.py --poisson none
python run.py --poisson ../1.example_build/build/result/poisson_only.npz
```

The third command requires that STEP01 saved `poisson_only.npz`, containing
`phi_v`, `x_nm`, `y_nm`, `z_nm`. It interpolates to the original 7x5 transverse
mesh. If missing, use `--poisson none`; this is an explicitly **flat-potential
numerical benchmark, not a biased FET result**. With no flag, `auto` uses the
Poisson file if present, otherwise uses the flat reference. You can point to
any other saved `poisson_only.npz` explicitly.

## Reviewable deliverables (results/)

| Artifact | What it proves/what to look for |
|---|---|
| `01_transverse_modes.png` | 2D modal probability densities; normalized, physically shaped eigenmodes |
| `02_subband_vs_x.png` | Local subband energies along transport; potential/barrier responds to Poisson input |
| `03_grid_convergence.png` | Ground state converges to the continuum infinite-well benchmark as grid refines |
| `04_mode_overlap.png` | Onsite `U^T U = I`; neighbouring slice overlap can change with potential |
| `05_input_potential.png` | Verifies which electrostatic potential was supplied (flat vs DEVSIM) |
| `grid_convergence.csv` | Numerical vs exact *discrete* FD eigenenergies and continuum limit |
| `validation.json` | Automatic PASS/FAIL numerical checks; unit tests run separately |
| `modes.npz` | Mode energies, eigenvectors and longitudinal coupling `hop`, reusable for STEP03 |

## GitHub figures: `figures/`

Every successful `python run.py` updates `figures/*.png`, `STATUS.md`,
`validation.json`, and `grid_convergence.csv`. This directory is **tracked**
by Git and can be reviewed as images in the GitHub UI (open `figures/STATUS.md`).
Raw `results/` output and `modes.npz` remain gitignored.

## Automatic numerical pass criteria

1. Exact discrete-Dirichlet FD eigenvalues agree with independently evaluated
   analytic FD formula within `2e-8 eV`.
2. Maximum `|U^H U-I|` below `1e-10` over all x slices and valleys.
3. Ground-state continuum confinement error below 1% at the finest sampled
   transverse grid (`47 x 39` interior nodes).
4. Every supplied and computed value is finite.
5. Additional unit tests: uniform potential shifts all eigenvalues by
   `-delta_phi` eV; thinner sheet increases confinement energy; flat-slice
   overlaps and hopping are consistent; 2D probability normalizes to one.

These are **not** claims of physical correctness of the full nanosheet device.
For STEP03, need mode-count convergence of **transport current** and 1D NEGF
contact self-energy validation. For STEP04, need a conservative real-space ↔
DEVSIM mapping and fully self-consistent charge conservation tests.

### Scientific cautions

- Silicon is treated as a hard-wall 2D rectangular well within the *existing*
  EMA model; no oxide wavefunction penetration yet.
- The default transverse `ny=7, nz=5` grid is coarse. Refinement above checks
  the *flat-well numerical scheme*. It does **not** refine the Poisson input.
- Valley masses/degeneracy are unchanged from `1.example_build/config.py`.
- `phi_v` from STEP01 is a Poisson-only result with assigned charge, **not**
  a self-consistent quantum transport potential.
- `U_i^H U_{i+1}` is not guaranteed to equal identity when basis changes with x.
  Degenerate modes may rotate; mode tracking/gauge control will need review
  before using these as accurate coupled-mode transport elements.
- Do not report an NEGF current from STEP02. Do not interpret numerical
  precision as agreement with calibrated experiment.

## GitHub

Commit source, README, tests and `figures/` for review.

```bash
git add 2.modespace/run.py 2.modespace/tests 2.modespace/README.md 2.modespace/.gitignore 2.modespace/figures
git commit -m "Add STEP02 transverse mode-space validation"
git push origin <YOUR_BRANCH>
```

`results/` (large numerical artifacts) is ignored by default.
