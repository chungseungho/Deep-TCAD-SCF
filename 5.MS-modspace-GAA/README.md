# STEP 05: Raised-S/D GAA finite-electrode EMA-Poisson

독립적인 **연구용 reference simulator**다. STEP01-04를 import하거나 수정하지 않는다.
기존 Si EMA masses와 valley/spin multiplicity는 그대로 사용한다.
실험에 맞춘 mass/doping fitting이나 임의 I-V 함수는 없다.

가장 중요한 구분은 다음과 같다.

* **MS-DFT-inspired occupation 연구:** `--backend ms`. 유한 고유상태의 L/C/R
  weight와 논문의 rule I을 검사한다. 모호한 standing/bound state를 임의로
  L/R에 할당하지 않고 중단한다. 엄밀한 MS constrained-search DFT 구현이 아니다.
* **실행 가능한 default transport:** `--backend spectral`. 유한 L/C/R 전체와
  명시적인 외부 reservoir battery를 사용한다. complex eigenstates + analytic
  Fermi integrals로 density를 구하므로 SCF에 **energy-resolved GF, lesser G,
  RGF, semi-infinite surface-GF 반복이 없다**. bath embedding 자체가 사라진
  것은 아니다. 이 density/transport 모델은 원 논문의 MS occupation과 구분한다.
* **비교 연구용:** `--backend lindblad`. isolated electrode eigenmodes에 대한
  Markovian filling과 Lyapunov equation. 좁은 공명에서 ill-conditioned solution,
  reservoir discretization/thermal broadening 문제를 감지하면 실패한다.

`SCF_CONVERGED_PROTOTYPE`는 지정한 수치 fixed point의 수렴만 의미한다.
`quantitative_validation_pass`는 현재 항상 false이다. S/D neutrality, mesh,
mode, finite-reservoir length, broadening 수렴이 검증되기 전에는 학습 정답이나
측정 수준의 TCAD로 인증하지 않는다. 실제 실행 범위는 [STATUS.md](STATUS.md)를 본다.

## A-E 설계 결정

| 과제 | 구현 | 검증과 한계 |
|---|---|---|
| A. 실제 3D GAA | 20 x 15 nm raised Si S/D, 12 x 5 nm nanosheet bridge, 24 nm gate, 2 nm oxide annulus, 양쪽 3 nm SiN spacer | OCC boolean fragment, 물질별 체적, 공유 면, 네 방향 gate surface, 양의 tetra quality 검사. metal volume 대신 연속 외측 oxide Dirichlet boundary |
| B. L-C-R coupling | 공통 cell-centred 공간 격자 위 각 단면 basis를 zero-extend. `H_i,i+1=-tx U_i^T U_i+1`; raised face의 missing-neighbor hard wall도 projection | full basis를 사용한 작은 3D real-space Hamiltonian과 고유값/행렬 동등성 검사. 같은 6x6 identity coupling으로 연결하지 않음 |
| C. MS occupation | L/C/R weight, L/R 분리 chemical potential, machine-degenerate subspace의 region-projector diagonalization | PNAS rule I을 참조하지만 EMA variational theorem을 주장하지 않음. 모호한 비평형 origin은 실패. default spectral backend는 명시적 bath에서 occupation을 결정 |
| D. GF-loop-free SCF density | finite `H_eff=H-i Gamma/2`의 eigensystem과 digamma Fermi integral | 독립적인 직접 GF energy integral과 일치 검증. 모든 retained state 사용. sparse transverse eigensolver, finite density matrix는 dense complex eigensolve |
| E. DC current | 동일 density matrix의 off-diagonal coherence에서 채널 bond flux 추출. T(E)는 수렴 후 별도 postprocessing | undriven channel 양끝/내부 current continuity, zero bias, gauge, reversal, 독립 Landauer 적분 검사. external bath 없는 폐쇄계 DC current로 해석하지 않음 |

상세 식과 논문에서 지원되는 내용/새 가정은 [THEORY.md](THEORY.md)에 정리했다.

## 실행

```bash
cd /home/class/2026_SeungHo/TCAD/Deep-TCAD-SCF/5.MS-modspace-GAA
source ../../.venv/bin/activate
python -m pip install -r requirements.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

python -m pytest -q tests
python run.py --phase mesh --smoke --out output/geometry
python run.py --phase poisson --smoke --out output/poisson

# 빠른 전체 소자 smoke: 매우 coarse quantum grid. production 물성 수렴이 아님.
python run.py --phase point --smoke --vg .8 --vd .1 --out output/smoke

# 실행시간을 줄인 실제 I-V 계산 설정. 아직 mesh/mode convergence 인증 전.
python run.py --phase idvgs --smoke --dx 1 --modes 3 --sd-modes 6 \
  --preconditioner none --out output/iv
python run.py --phase idvds --smoke --dx 1 --modes 3 --sd-modes 6 \
  --preconditioner none --out output/iv --resume

# 기본 전압: Id-Vgs 0:1.2:0.05 at Vds=.1; Id-Vds 0:1.2:0.05 at Vgs=.4,.8,1.2
# 원래 default의 더 많은 modes/더 세밀한 Poisson mesh:
sbatch job.sh --phase idvgs --out output/reference
sbatch job.sh --phase idvds --out output/reference --resume

# energy-resolved transmission은 SCF 완료 후에만 계산
python run.py --phase point --smoke --dx 1 --modes 3 --sd-modes 6 \
  --vg .4 --vd .1 --spectra --out output/spectral_check

# 수렴된 동일 state에서 energy quadrature만 독립적으로 정밀화
python postprocess.py output/spectral_check vg0.400000_vd0.100000

# 별도 단일 바이어스 수렴성 campaign: 실제 SCF를 모두 새로 실행
python diagnostics.py --out output/convergence
python diagnostics.py --out output/convergence --cases \
baseline,mesh_medium,mesh_fine,transverse_fine,longitudinal_fine,modes_low,modes_high,source_drain_24nm,source_drain_36nm,gamma_5mev,gamma_20mev,basis_local,mixing_unscreened --resume
```

실행 때 실제 사용한 버전은 `requirements-lock.txt`에 있다. Linux native dependencies로
Gmsh의 `libXft.so.2`와 DEVSIM의 BLAS/LAPACK가 필요하다. DEVSIM이 math library를
찾지 못하면 설치된 OpenBLAS 경로를 `DEVSIM_MATH_LIBS`에 지정한다
(예: `export DEVSIM_MATH_LIBS=libopenblas.so.0`). 이 의존성은 Python wheel만으로
항상 설치되는 것이 아니다. 기존 HPC의 작동하는 TCAD environment를 사용할 수 있다.

`--basis fixed`가 기본이다. x별 형상에 맞는 fixed transverse kinetic basis에서
전위의 **full projected matrix**를 매 iteration 갱신한다. `--basis local`은
slice별 eigenmodes/overlap을 다시 구한다. 둘 사이의 density/current 차이를
검사해야 한다. retained/omitted mode 사이 crossing은 local truncation을
불연속으로 만들 수 있으므로 mode count 수렴이 필요하다.

`--preconditioner screened`가 기본이다. local `n/kBT` susceptibility는
quasi-Newton **preconditioner에만** 사용한다. 실제 fixed-point residual은
screening이 없는 원래 Poisson-quantum map으로 계산한다. `none`과 비교 가능하다.
Anderson의 부적절한 방향은 positivity guard로 거부하고 adaptive damping한다.
실제 residual이 정체되면 원래 Poisson-quantum potential map에 대한 matrix-free
Newton-Krylov와 Armijo line search로 전환한다. 이때도 동일한 phi/charge/current
기준을 적용한다. Jacobian-vector products도 실제 quantum 계산으로 평가한다.

**Longitudinal bandwidth:** finite-difference EMA의 각 valley bandwidth는
`4 hbar^2/(2 mx dx^2)`이다. `dx=1 nm`에서 최소값은 약 0.166 eV로, 높은 Vds
곡선의 인공적인 current suppression/NDR을 신뢰할 수 없다. default `dx=.25 nm`는
이 최소 bandwidth를 약 2.66 eV로 늘린다. 이것만으로 공간 수렴이 증명되지는 않는다.
위 빠른 I-V 명령의 `--dx 1`은 수치 실습/성능용이며 publication용이 아니다.
`bias_inside_minimum_bandwidth`는 보수적 점검 항목이며 이 값 하나로 인증하지 않는다.

`--resume`는 per-bias checkpoint를 복원한다. physics/numerics configuration이
다르면 기존 output에 섞어 넣지 않는다. 실패 지점은 traceback/checkpoint/history를
남기며, 기본적으로 sweep을 중단한다. `--keep-going`을 지정해 실패를 유지한 채
다음 바이어스를 시도할 수 있다. plot은 수렴한 지점만 사용한다. 따라서 실패한
지점이 있는 sweep은 일부 곡선만 나타날 수 있다. `--max-points`는 smoke 범위다.
새 실행은 package versions와 simulation source SHA256도 저장한다. 코드가 바뀐
출력 폴더에서 resume하면 거부한다. 개발 중 생성된 legacy 결과는 개별 settings와
history를 보존했으며 일부는 source digest 도입 전 결과이다.

## Poisson와 두 격자 사이 mapping

이 형상의 coarse tetra mesh에서 DEVSIM native `NodeVolume` 합은 CAD 체적과
다를 수 있었다. 따라서 이 버전은 native absolute circumcentric model을
그대로 사용하지 않고 다음 **명시적인 DEVSIM discretization**을 설정한다.

* Original tetra P1 gradient에서 signed integrated edge flux를 조립한다.
  `edge_couple_model=SignedElementCouple`.
* 양의 barycentric nodal volume `sum(V_tet/4)`를 실제 equation integration에
  지정한다. `node_volume_model=PhysicalControlVolume`.
* 이는 linear tetra Galerkin stiffness와 동등한 conservative nodal edge balance다.
  native circumcentric FVM과 같은 scheme이라고 주장하지 않는다. signed edge
  conductance가 존재할 수 있어 discrete maximum principle도 보장하지 않는다.
* default/native volume 차이는 `mapping.json`에 물질별로 기록한다. geometric
  volume만 따로 적분하고 실제 Poisson은 다른 volume을 쓰는 방식이 아니다.
* potential은 **원래 silicon tetrahedra**의 barycentric interpolation으로 보낸다.
* density는 local positive kernel의 두 marginals를 control volume과 quantum
  cell volume에 맞춘 balanced map으로 보낸다. constants, positivity, total
  charge를 함께 보존한다. 길이 단위는 quantum/CAD nm, DEVSIM cm이다.
* 이 density map은 affine potential map의 정확한 adjoint가 아니다. work defect를
  기록하며 에너지 변분 원리를 주장하지 않는다. kernel의 유한한 공간 smoothing도
  mapping/mesh 수렴에서 평가해야 한다.

FD grid refinement로 coarse Poisson mesh 오류가 사라진다고 해석하지 않는다.

## 결과 파일

| 경로 | 내용 |
|---|---|
| `mesh/gaa_nm.msh`, `gaa_cm.msh`, `gaa.brep` | 실제 mesh/CAD, nm 표시용과 cm DEVSIM 입력 |
| `mesh/geometry.json` | 물질 체적, physical group/surface, mesh quality |
| `mapping.json`, `contact_calibration.json` | 전하/상수/선형 mapping 검사, 적분 체적, μ calibration |
| `idvgs.csv`, `idvds.csv`, `point.csv` | 상태/잔차/neutrality/current 포함 bias 결과 |
| `scf_history.csv` | **unmixed** phi/charge residual, N, occupation/eigenvalue 변화, runtime |
| `states/vg*_vd*.npz` | 3D potential/n, node data, subbands, modes, state weights, L/R-resolved n, flux |
| `states/*.json` | bath width, density-matrix positivity, current continuity, reservoir spacing |
| `electrostatics.vtk*` | ParaView용 모든 material의 실제 DEVSIM solution |
| `figures/` | mesh cut/material/doping, fields, modes, occupations, SCF, I-V |
| `transfer_metrics.json` | 명시적 constant-current Vth, SS, gm, Ion/Ioff; prototype 지표 |

Gate-perimeter normalization은 `2(width+thickness)`이며 Vth default criterion은
`1e-7 A/um`, 즉 이 geometry에서는 `3.4e-9 A`이다. SS는 threshold 이하,
`1e-16 A` 이상에서 adjacent positive logarithmic slope로 계산한 최소값이다.
실험값으로 calibration한 수치가 아니다.

## DeepSCF로 연결

현재 저장소에는 연구실 DeepSCF network/weights가 없으므로 그 신경망을
이식하거나 학습이 완료됐다고 주장하지 않는다. 먼저 신뢰 가능한 teacher가 필요하다.

```bash
# 현재 결과는 prototype이므로 명시적으로 허용해야 export된다.
python data.py output/iv output/training.npz --include-prototypes
python surrogate.py output/training.npz output/rbf_baseline.npz
```

Dataset은 geometry/grid/physics가 같은 **전체 bias point 단위**로 train/holdout을
분리한다. target은 channel density, one-shot input 후보는 bias와 donor-neutral
non-SCF Poisson field이다. 수렴한 potential은 입력으로 섞지 않는다. 함께 저장하는
S/D-origin density는 thermodynamic local QFL이라는 뜻이 아니다.

`surrogate.py`는 bias 두 개를 사용하는 RBF-ridge log-density baseline일 뿐,
DeepSCF가 아니다. train/holdout error를 내고 geometry/parameter mismatch와
bias extrapolation을 거부한다. `--warm-start-model`로 사용해도 reference SCF를
끝까지 실행한다. 아직 one-shot 정확도/전류 보존이 검증되지 않았으므로 quantum
solver를 완전히 끄는 production 경로는 제공하지 않는다.

다음 학습 단계는 수렴성 campaign으로 teacher 오차를 확정한 뒤 lab DeepSCF에
`(geometry, doping, boundary biases, initial phi) -> n_channel`을 연결하고,
heldout 전체 bias/device에서 Poisson residual과 reference current를 검증하는 것이다.
