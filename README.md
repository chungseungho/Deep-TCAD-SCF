# Deep-TCAD-SCF — Mono-Nanosheet GAA reference v0.1

**목표:** Mono-nanosheet GAA의 **Gmsh 3D geometry → DEVSIM 3D Poisson → 2D transverse eigenmodes → ballistic coupled mode-space NEGF → Poisson–NEGF SCF** 구현. 추후 DeepSCF-inspired AI에 비교 기준으로 사용하기 위한 prototype입니다.


## 설치 및 실행 (Linux)

```bash
cd /home/class/2026_SeungHo/TCAD/Deep-TCAD-SCF
source ../.venv/bin/activate  # .venv가 TCAD 디렉터리에 있을 경우
python -m pip install numpy scipy matplotlib devsim gmsh

python tests/test_negf.py
python run.py --mode mesh
python run.py --mode poisson --vg 0.4 --vd 0.05 --reuse-mesh
python run.py --mode negf-only --vg 0.4 --vd 0.05
python run.py --mode coupled --vg 0.4 --vd 0.05 --iterations 8 --reuse-mesh
```

**메시 파일 확인:** `gmsh build/mono_nanosheet_nm.msh` (geometry는 nm, DEVSIM용은 `*_cm.msh`). CLI `gmsh`가 설치되어 있지 않더라도 Python gmsh API만으로 mesh 생성 가능.

## 구조

```text
Deep-TCAD-SCF/
├── config.py                     # nm geometry, initial EMA parameters, energies
├── run.py                        # build/poisson/NEGF/self-consistency
├── geometry/make_mesh.py         # 3D conformal Gmsh physical groups
├── poisson/devsim_poisson.py     # 3D Poisson, fixed injected electron density
├── negf/mode_space.py            # 2D eigenmodes + RGF / G^< + density/current
├── tests/test_negf.py            # pure NumPy/SciPy transport tests
└── build/                        # meshes, field/current outputs
```

## 기본 형상 및 물리 가정

- x = transport, y = nanosheet width 12 nm, z = sheet thickness 5 nm.
- Source 12 nm / channel (gate) 24 nm / drain 12 nm. Gate oxide 2 nm is a rectangular GAA annulus; no metal volume (gate boundary imposed on outer oxide surfaces).
- Silicon is one DEVSIM region with source/channel/drain spatial doping; gate oxide separate region. Source and drain contact boundaries at opposite x planes.
- No DEVSIM DD physics. Charge is exported from ballistic NEGF to DEVSIM Poisson as a fixed node-solution.
- Dirichlet hard wall in transverse Schrödinger; piecewise uniform effective masses per Si Delta valley, `+/-` multiplicity 2, spin 2. **Values in config.py are illustrative EMA inputs; not calibrated and not copied from any previous EMA source**. 사용자가 검증한 EMA parameter와 orientation으로 변경해야 함.
- Longitudinal coupled mode space uses basis overlap `U_i^† U_{i+1}`; **mode truncation errors must be checked**. Electron density uses local cross-sectional mode basis and injected spectral functions (lesser-GF equivalent).
- Ballistic transport only; no phonons, traps, nonparabolic bands, workfunction differences, finite-bandwidth correction, or degeneracy-adjusted bulk carrier statistics. Gate is an ideal Dirichlet contact, metal workfunction not modeled.
- Reservoir reference `mu_S = 0.15 eV` relative to flat-band conduction band edge and `mu_D = mu_S - Vd`. These values are explicit assumptions, not fitted to doping equilibrium.
- Poisson model uses `rho=q(N_D - n)` and ignores holes/ionization physics. Contact neutrality is **not guaranteed** by the specified reservoir chemical potentials. Proper source/drain built-in band alignment and neutrality calibration are required before quantitative transistor claims.

## 출력

- `build/result/solution.npz`: x/y/z mesh, potential, NEGF charge, transmission.
- `build/result/history.csv`: iteration, Id, max field difference, elapsed time.
- `build/result/result.txt`: simulation conditions and current.
- `build/result/overview.png`: potential slice and transmission plot.
- `build/result/potential_3d.vtk`: ParaView visualization (if supported).

## 중요한 검증 TODO (차례대로)

1. `--mode mesh`: physical groups and 3D tetra connectivity inspect; physically closed oxide gate boundary and Si/Ox interface.
2. `--mode poisson`: continuity of normal displacement across Si/Ox; Dirichlet values; zero charge Laplace sanity.
3. `tests/test_negf.py`: perfect 1D chain transmission and zero-bias I=0 (implemented).
4. Increasing transverse mesh resolution and mode number: eigenenergies/current convergence; real-space NEGF reference comparison.
5. Fixed 3D↔mode-space mapping: charge conservation; contact neutrality at zero bias; stable SCF residual; adaptive energy integration.
6. Bias-dependent Id–Vgs / Id–Vds; match fully SCF and evaluate reproducibility before any AI training.

## AI 분리 원칙
먼저 fully-SCF reference data 생성과 공간/수송 수렴성을 검증합니다. 이후 NEGF-derived `n_SCF(x,y,z)`를 density target으로 저장해 one-shot charge prediction 모델을 학습합니다.
