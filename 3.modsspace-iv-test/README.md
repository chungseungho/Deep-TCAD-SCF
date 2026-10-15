# STEP 03 — Mono-Nanosheet GAA: Full Device I–V Validation

**Project:** `Deep-TCAD-SCF/3.modsspace-iv-test/`  
**State:** Prototype implementation with quantum **unit tests passed**; end-to-end Gmsh/DEVSIM/NEGF device sweep must still be run and validated on the Linux server.

This is **not** a DD MOSFET and does **not** use learned charge. This step builds reference fully-self-consistent 3D Poisson + ballistic coupled mode-space NEGF data for a **single silicon nanosheet**. Existing STEP01 EMA parameters (Si anisotropic valley masses, spin/valley factors) remain unchanged.

## Geometry / physics / mesh

- Silicon sheet: **W = 12 nm**, **T = 5 nm**, **Lg = 24 nm**; source/drain 12 nm each; oxide 2 nm, wrap-around GAA contact.
- Source/Drain donor `1e19 cm^-3`; nominal channel `0 cm^-3`; `T=300 K`.
- Si EMA masses and constants read from `../1.example_build/config.py`. **Do not modify EMA to force desired I–V results.**
- 3D Poisson: Gmsh + DEVSIM. Silicon density is from **NEGF's lesser GF**, no DD involved.
- Transverse (y,z): **31 × 25 internal FD nodes**, 6 fixed EMA modes per valley; this is a **fixed-basis, coupled** mode-space Hamiltonian with projected *off-diagonal* electrostatic mode coupling.
- Longitudinal x: **345 nodes**, Δx ≈ 0.140 nm. Required because `0–6 V` scans otherwise create **artificial lead-band cutoff** in the discrete EMA Hamiltonian. For high-bias checks the minimum valley lead bandwidth must exceed `Vds_max + 1 eV`.
- 3D Gmsh mesh **1 nm target size** (new mesh in `output/mesh/`; the frozen STEP01 mesh is not overwritten).
- Ballistic, coherent, ideal leads, hard-wall transverse confinement, no scattering, no leakage, no self-consistent hole transport.
- `--gate-flatband 0.5`: explicitly assumed *effective gate offset*, **NOT measured/calibrated**. In the Poisson boundary condition `Vgate_eff = Vgs - gate_flatband`. This value is provisional; extracted Vth is model-dependent. Effective contact Fermi level likewise inherits the illustrative STEP01 setting.

### CRITICAL PHYSICAL LIMITATIONS

**6 V is a requested numerical stress scan, not a safe/realistic operating voltage for 2 nm SiO₂.** Real oxide breakdown, tunneling, gate current and inelastic scattering are absent. **Do not present 0–6 V as an experimentally validated GAA operational window.** An apparent plateau or roll-off is not sufficient proof of saturation: verify source/drain spectra, mesh, mode count and energy integration convergence.

### High drain bias electron-density integration

Current uses source-and-drain transmission. Charge uses the *union* of the occupied propagating source **and drain** spectra, because under high Vds the drain-only states also contribute to density, even when their contribution to current transmission is zero. Integrating only the lead spectral overlap gives a wrong Poisson charge.

## Directory prerequisites

```text
TCAD/
├── .venv/
└── Deep-TCAD-SCF/
    ├── 1.example_build/
    │   ├── config.py
    │   ├── geometry/make_mesh.py
    │   └── poisson/devsim_poisson.py   # previous working NetDoping fix
    ├── 2.modespace/
    └── 3.modsspace-iv-test/
        ├── job.sh
        ├── run.py
        ├── quantum.py
        ├── quality.py
        ├── analysis.py
        ├── plotting.py
        ├── tests/
        └── figures/
```

> `1.example_build` must contain the working v0.1.1 Poisson module with initialized `NetDoping`. **Do not silently substitute the older broken version.**

## Run on Linux cluster

```bash
cd /home/class/2026_SeungHo/TCAD/Deep-TCAD-SCF/3.modsspace-iv-test
source ../../.venv/bin/activate
python -m unittest discover -s tests -v
sbatch job.sh
```

To inspect queue / job logs:

```bash
squeue -u "$USER"
tail -f slurm-<JOB_ID>.out
```

### Staged smoke run on your server (recommended before spending an entire 24h allocation)

This uses real DEVSIM+NEGF, not a synthetic I–V curve; it only runs the first **one bias**.

```bash
sbatch job.sh --phase idvgs --max-points 1
```

Then resume into **the same `output/` directory**:

```bash
sbatch job.sh
```

**IMPORTANT:** `--max-points 1` is a partial run and does **not** establish Vth or Id–Vds. The full job runs **31** Id–Vgs points (0–6V, step 0.2V at Vds=0.1V) then extracts threshold and picks three Vgs values below/just above/well above threshold. It subsequently runs **3 × 31 = 93** Id–Vds points (0–6V, step 0.2V).

Default 3 gate biases, rounded to 0.2 V, are:

- Below Vth: `Vth - 0.4 V`
- Just above Vth: `Vth + 0.4 V`
- Strong on: `Vth + 1.4 V`

**These are determined after the computed Id–Vgs curve, never hard-coded to 0.4/0.6/0.8 V.** If Vth is outside the accessible window, the simulation stops with an explanation rather than inventing three conditions. A manually justified Vth can be supplied with `--manual-vth <value>`; the override is explicitly logged.

### Resume / options

- `--out output`: stable checkpoint folder across separate SLURM jobs, `output/` is `.gitignore`d.
- `--phase idvgs`, `--phase idvds`, or `--phase all`.
- `--max-points 1`: smoke test and stop after first new point.
- `--energy-step 0.005`: finer real-axis energy integration.
- `--modes 6 --ny 31 --nz 25 --nx 345`: validated basic resolution; check convergence figures before accepting physical output.
- `--regenerate-mesh`: recreate Poisson mesh under this stage's output. (New `--out` is required if numerical settings change.)
- `--mix 0.20 --max-iterations 35`: conservative initial self-consistency parameters.
- `--allow-unconverged`: diagnostics only; an unconverged point is never silently promoted into the final *validated* Id–Vgs dataset.
- `--gate-flatband 0.5`: model assumption, not a fit. **Changing settings requires a new output directory** to avoid mixing incompatible data.

Example of a separate debugging run:

```bash
sbatch job.sh --out output_debug --vg-max 0.6 --vd-max 0.6 --step 0.2 --nx 345
```

On a shared cluster do not run extra long sweeps interactively on login nodes. Use `sbatch`.

## Expected files — GitHub review workflow

Large checkpoint files are in `output/states/*.npz`, with other data in `output/*.csv` and `output/*.json`. For GitHub, completed plots and status are **automatically copied** into the code directory after each finished bias:

| GitHub file | Physical evaluation |
|---|---|
| `figures/01_transverse_modes.png` | Transverse \|χ\|² of actual 31×25 EMA eigenbasis |
| `figures/02_idvgs.png` | Id–Vgs linear/log and gm, extracted Vth |
| `figures/03_idvds.png` | 3 Id–Vds curves for below/near/strong ON Vgs |
| `figures/04_device_maps.png` | Actual SCF potential and log electron density (x–y, x–z) |
| `figures/05_scf_convergence.png` | Max Δphi and raw NEGF-charge fixed-point residual |
| `figures/06_energy_convergence.png` | Actual SCF potential: one-shot NEGF vs energy integration step |
| `figures/07_mode_convergence.png` | Actual SCF potential: mode count sensitivity |
| `figures/08_subbands_under_bias.png` | Local subbands from actual SCF potential |
| `figures/quality_checks.json` | Energy and mode test values / flags |
| `STATUS.md` | Current step status and figures in Markdown |

After your SLURM run completes, GitHub can display `3.modsspace-iv-test/STATUS.md` and all figures directly:

```bash
cd ..
git add 3.modsspace-iv-test/
git commit -m "STEP03: Mono nanosheet GAA coupled MS-NEGF IV validation"
git push origin HEAD
```

## Definition of PASS

1. All **8 quantum unit tests** pass (this alone is **not** a device validation).
2. At least 31/31 fully **SCF-CONVERGED** Id–Vgs points, followed by 93/93 fully converged Id–Vds points; `NOT_CONVERGED` is never represented as successful.
3. Gated **Id–Vgs** has a credible OFF→ON transition and a defensible constant-current Vth crossing. 0.2 V step is **not adequate by itself for precise SS extraction** — refine around threshold before publication.
4. Three Id–Vds curves exist at the selected below/above/strong-ON Vgs, with correct near-zero Id at Vds=0, monotonic bias behavior in the physically reasonable range, and an understood high-bias response. Full 6 V stress curve is not a calibrated operating range.
5. Transport spectral integration: last two steps (5 meV vs 2.5 meV) agree within a proposed **3%** target on a real converged snapshot; 6 vs 8 mode currents within a proposed **5%** target. These are *fixed-potential tests*, not proof of full self-consistent convergence for different retained modes.
6. Before a paper: systematic **Poisson mesh convergence**, transverse-grid convergence, longitudinal-grid convergence, carrier-mapping charge conservation, proper contact equilibrium and workfunction calibration, and benchmarking against real-space quantum transport must also be completed. This code is a starting framework, not a finished validated industrial-grade device simulator.

## Important implementation notes

- A `CONVERGED` point requires **max potential change**, **charge fixed-point residual** and **relative current change** tolerances; it does not mean physically calibrated.
- The fixed transverse basis allows inter-mode electrostatic coupling without eigenvector phase ambiguities between slices. It is **not** the slice-local basis of STEP02; both describe the same EMA transverse kinetic operator, but their truncation errors differ.
- No carrier transport uses the DEVSIM drift-diffusion equations.
- `DEVSIM` and `gmsh` are not installed in the development runtime that generated this bundle, so the *full coupled 3D device run has not been executed here*. It must be tested in your Linux server; this is explicitly separate from the 8 quantum unit tests we did execute.
