# STEP 06 — 3D GAA EMA–Poisson quantum TCAD

실제로 실행되는 독립 소자 simulator다. `STEP 01–05`를 import하거나 수정하지 않는다.
신경망, surrogate, AI training은 포함하지 않는다. 현재 달성 범위와 실제 계산 수치는
[STATUS.md](STATUS.md), 검증 근거는 [VALIDATION.md](VALIDATION.md)를 먼저 확인한다.
이전 단계의 발전 과정, 실제 STEP05 CSV 대조와 코드 결함/모델 한계 구분은 [AUDIT.md](AUDIT.md)에 정리했다.

## 이번 단계의 핵심 수정

1. **삼중 물질 접합의 전위 연속성:** STEP05의 다중 interface 방식에서 실제 공유 노드의 전위 점프를 발견했다. STEP06은 Gmsh의 물질 구분을 유지하면서 **공통 P1 전위 자유도**를 DEVSIM에 조립한다. 각 tetrahedron의 유전율을 stiffness에 반영하므로 Si/oxide/spacer를 같은 유전율로 취급하는 방식이 아니다.
2. **정확한 두 격자 교환 적분:** tetrahedron과 quantum cell의 교집합에서 shape function을 적분한다. 양의 밀도, 두 적분 체적, 균일 밀도, 선형 전위, 전위–전하 work adjoint를 함께 검사한다. nonequilibrium variational 원리를 입증했다는 의미는 아니다.
3. **두 접촉 모델:** default `spectral`은 finite S/D 전체에 명시적인 wide-band reservoir를 연결한다. complex poles/analytic Fermi correlation으로 SCF 내부 에너지별 GF 반복을 사용하지 않는다. `negf`는 semi-infinite wire lead와 contour/real-axis GF 적분을 사용하는 독립 비교 경로다.
4. **수렴 판정:** unmixed potential/charge residual, current change, channel flux와 전체 electrostatic Gauss balance를 검사한다. 실패는 CSV/traceback에 남기고 성공 곡선에서 제외한다.
5. **에너지 유효 범위:** 단순히 Vds와 bandwidth만 비교하지 않고 실제 local subband와 occupied energy window의 valley별 FD dispersion 오차를 기록한다.

**현재 finite-bath 결과는 정량 검증된 AI teacher가 아니다.** Contact neutrality와 mode/grid/bath 수렴이 남으면 `quantitative_validation_pass=False`로 저장한다. Semi-infinite reference의 실제 비평형 적분 실패도 결과에서 숨기지 않는다.

## Geometry / physics

- 12 × 5 nm 단일 Si nanosheet, Lg=24 nm, Ls=Ld=12 nm, spacer=3 nm, tox=2 nm.
- Raised S/D: 20 × 15 nm. x 방향 수송. 위/아래/양옆 4면 gate.
- eps(Si/SiO2/SiN)=11.7/3.9/7.5, T=300 K.
- ND(S/D)=1e19 cm^-3, channel intrinsic, transition=1 nm. Fully ionized donors.
- Si Delta valleys: ml=.916 m0, mt=.19 m0. 각 +/- pair 2, spin 2를 한 번씩만 적용한다.
- Ballistic electron-only Hartree EMA; hard-wall Si surfaces. Hole, phonon, trap, recombination, oxide tunnelling, XC/nonparabolicity는 모델에 없다.
- Flatband=.5 V는 inherited assumption이며 실험 fitting 값이 아니다.

## Quick start

```bash
cd 6.MS-modspace-model
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m pytest -q tests
python run.py --phase point --smoke --dx .5 --modes 3 --sd-modes 6 --vg .8 --vd .1 --out output/point
```

위 quick point는 실행 점검용이다. `.5 nm` longitudinal grid의 정량 수렴을 의미하지 않는다.
실제 sweep과 HPC 명령, 비용은 [HOW_TO_RUN.md](HOW_TO_RUN.md)를 따른다.

## 계산 모듈

| 모듈 | 역할 |
|---|---|
| geometry/build.py, validate.py | OCC Boolean conformal mesh, 물질 체적/공유면/4면 gate 검사 |
| poisson/device.py | 공유 전위 공간의 실제 DEVSIM Poisson, material-weighted P1 stiffness, Gauss 검사 |
| poisson/discretization.py | 독립 tetra gradient/stiffness 적분; native FVM과 구분 |
| coupling/cell_integrals.py, mapping.py | 실제 cell/tetra 교집합 적분과 adjoint charge/potential map |
| modespace/hamiltonian.py | anisotropic transverse eigenmodes, full potential projection, variable-section overlap |
| multispace/occupation.py | closed finite L/C/R spatial-weight 진단; 모호한 active origin 거부 |
| transport/finite_spectral.py | finite wide-band bath의 analytic density matrix |
| transport/negf_reference.py, quadrature.py | semi-infinite lead reference, contour와 적응형 injectivity 적분 |
| solver.py, coupling/mixing.py | 실제 fully SCF, Anderson/Newton–Krylov, unchanged fixed point |
| run.py, provenance.py, resources.py | CLI, 재개 검사, configuration/source identity, 비용 추정 |
| plotting.py, export_fields.py | CSV/NPZ 기반 그림, DEVSIM/quantum ParaView VTU |
| validation/ | 이전 단계 감사, reduced benchmark, 실제 SCF 수렴 campaign |

## 출력

`point.csv`, `idvgs.csv`, `idvds.csv`의 status가 기준이다. `ERROR`/`NOT_CONVERGED`는 검증된 소자 결과가 아니다.
`states/*.npz`에는 geometry/grid/doping/bias/T/masses를 복원할 config JSON, initial/converged phi, quantum n, Poisson input n, subbands, current, residual과 fingerprint/SHA256이 들어간다.
`states/*.json`에는 Pauli/current/neutrality/Gauss/dispersion 진단을 저장한다.
Default finite-spectral 경로의 L/R-origin n은 reservoir injectivity이며 local thermodynamic QFL이라는 뜻이 아니다.
NEGF contour 경로의 equilibrium/bias-injected density는 full L/R-origin split과 구분한다.

`electrostatics.vtk*`는 DEVSIM 전체 물질 전위이며 **좌표는 cm**다. `quantum_fields_nm.vtu`는 마지막 accepted bias의 quantum cell fields이며 **좌표는 nm**다. ParaView에서 둘을 겹치려면 native DEVSIM geometry에 1e7 배율을 적용한다. Potential V, Ec eV, density cm^-3, current A를 사용한다.
원격 push, 학교 HPC 제출, 공용 Python/MPI/CUDA 변경은 수행하지 않는다.
원시 NPZ/VTK/mesh는 전달 ZIP에 포함하며, 생성된 대형 파일과 private 환경·캐시는 기본 source-control 대상에서 제외하도록 구성했다. CSV·보고서·재현 설정은 별도로 확인할 수 있다.
