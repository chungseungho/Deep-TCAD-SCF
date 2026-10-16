# STEP06 — actual implementation and validation status

**Level 3: 실제 3D GAA fully SCF를 검증했다. Level 4/5 정량 reference 또는 AI teacher는 아직 아니다.**

| 판정 | 실제 근거 |
|---|---|
| Level 1 — Code implemented | 독립 Gmsh–DEVSIM–EMA solver, finite spectral / NEGF comparison |
| Level 2 — Numerical tests passed | 44 passed in 1.29s |
| Level 3 — Device SCF converged | 실제 GAA에서 unmixed phi/n/current, Gauss 검사 |
| Level 4 — Numerical reference | 미달: contact/grid/mode/basis 등 전 영역 plateau 미확립 |
| Level 5 — AI teacher ready | 미달: 모든 state의 quantitative_validation_pass=False |

| 실제 sweep | 계산 / 목표 | accepted SCF | 실패 |
|---|---:|---:|---:|
| Id–Vgs, Vds=.1 V | 25/25 | 25 | 0 |
| Id–Vds, Vgs=.4,.8,1.2 V | 75/75 | 75 | 0 |

현재 sweep 설정은 dx=.25 nm, dy=dz=1 nm, modes=2/4, gamma=10 meV, coarse Poisson mesh다.
이는 **실제로 계산된 provisional I–V**이며 mode/contact/continuum convergence를 인증한 소자 사양이 아니다.
별도 campaign에서 dx, gamma, mode, length, basis, Poisson/transverse mesh, SCF tolerance를 실제 coupled SCF로 비교했다. 현재 SCF tolerance 외의 축은 지정한 plateau 기준을 통과하지 못했다. Fine dx의 세 번째 coupled 계산은 미완료다.

## 검출하고 수정한 오류

- STEP05 원본을 saved bias Vgs=.4/Vds=.1 V에서 다시 읽어 실행했을 때 Si/oxide 공유 전위 점프가 0.0263087 V, oxide/spacer는 0.124627 V였다. 삼중 접합의 겹치는 interface 조건을 공유 P1 자유도로 교체했다. 물질별 epsilon은 유지한다.
- Tetra/cell 교집합 적분으로 potential/charge exchange를 통일했다. Uniform density, affine potential, total charge, work adjoint를 실제 GAA mesh에서 검사했다.
- Exact saved field의 physical fixed-point residual과 Gauss balance를 최종 판정에 추가했다.

## 독립 전류 검증

Vgs=.8/Vds=.1 V, dx=.5 nm, modes=3/6의 실제 GAA에서 finite-spectral bond current와 같은 finite-bath Hamiltonian의 독립 Schur/RGF transmission integral을 비교했다. 전체 I–V의 dx=.25 nm, modes=2/4 설정과 구분한다.
초기 energy panel .02/.01/.005 eV 각각에서 최대 상대 차이 1.25e-08였다.
이 비교는 동일 접촉 모델의 계산 일치이며 contact model 자체의 실험적/정량 검증은 아니다.

## 남은 문제

- Broadening, mode truncation, finite S/D와 basis/grid 의존성이 실제 수렴 표에 남는다. 한 쌍의 agreement로 plateau를 선언하지 않는다.
- Connected S/D neutrality는 isolated source calibration과 별개다. Finite-bath 결과의 오차가 남는다.
- Abrupt raised-S/D 접합의 독립적인 mode truncation이 dx 의존성을 키운다. Nested trace 보강 비교는 dx=.5/.25 nm 두 계산에서 전류 차이 0.474%였으며, 두 수준만으로 plateau를 인증하지 않는다. baseline I–V와 섞지 않았다.
- Semi-infinite NEGF의 실제 equilibrium SCF는 수렴했고 I=0이다. 실제 biased GAA는 좁은 공명 적분 기준을 통과하지 못했다. `results/negf_lowbias/`의 ERROR/traceback을 보존한다.
- Fine quantum grid만으로 EMA dispersion, transverse truncation, Poisson mesh 또는 contact 모델의 검증이 끝나지 않는다.
- Hole/scattering/trap/tunnelling/XC/nonparabolicity 및 calibrated workfunction은 포함하지 않는다.

원시 CSV/NPZ/설정/잔차/VTK와 세부 [VALIDATION.md](VALIDATION.md), [HOW_TO_RUN.md](HOW_TO_RUN.md)를 함께 제공한다.
STEP01–05는 보존했고 원격 push와 학교 SLURM 제출은 수행하지 않았다.
