# STEP11 실행 결과

**판정: Release Candidate. RELEASE_READY=false, TEACHER_READY=false.**

실제 코드·Native 계산·회귀 검사를 수행했다. 원래 Raised S/D의 정량 I–V 출시 목표는 아직 충족하지 못했다. 기존 STEP01–10 소스는 보존했다.

## 변경 사항

- Full-basis transverse Schrödinger–Poisson 전극 기준해를 Poisson boundary와 self-energy에서 함께 사용한다.
- 실제 basis overlap을 사용한 Σ=V g V†, block Hamiltonian, banded resonance guidance, 전체 오차 기준 constant-memory 적분을 구현했다.
- Backend-neutral domain 예외, Backend 소유 model-validity hook, 미지원 관측량 None/UNSUPPORTED, 엄격한 재시작 검사를 구현했다.
- 환경/단일점/sweep/refine/export/certify CLI와 자동 필수 게이트 certificate를 추가했다.

## 필수 게이트

| 항목 | 결과 |
|---|---|
| native_integration | PASS |
| backend_replacement | PASS |
| independent_quantum_reference | PASS |
| raised_nonequilibrium_scf | FAIL |
| gate_drain_iv | NOT_RUN |
| lead_device_matching_and_cut | FAIL |
| bound_state_completeness | INCONCLUSIVE |
| mode_convergence | FAIL |
| grid_convergence | FAIL |
| electrostatic_buffer_convergence | FAIL |
| energy_convergence | NOT_RUN |
| clean_install_and_resume | PASS |
| extracted_archive_execution | PASS |

## 실제 수치와 한계

- 최종 pytest: 92 passed, 0 skipped. Native Gmsh/DEVSIM 및 mock backend 포함. 독립 Quantum Reference 13개 사례 통과.
- 20×15 nm 전극의 constant-profile vs SP 밀도 차이 33.93%, μ 차이 8.094 meV. 상수 전위 근사는 채택하지 않았다.
- 전극 grid 1→0.5 nm의 μ 차이는 1.886 meV로 0.1 meV 기준에 미달한다.
- Fully-SCF 8/32→16/64: channel 밀도 0.485%, junction/extension 5.037%, reservoir 3.066% 차이. 평형 결과이며 수송 모드 수렴 인증이 아니다.
- Poisson mesh refinement: channel 1.127%, junction/extension 5.190%, reservoir 1.673% 차이.
- 동일 목표 영역을 유지한 6 nm buffer 확장: channel 0.436%, junction/extension 1.215%, original reservoir 2.207% 차이.

`results/assessment/computed_bias_points.csv`에 실제 SCF 상태를, `failed_bias_attempts.csv`에 실패를 구분했다. CSV/PNG에는 보간·피팅·실패점 대체값이 없다. 원래 요청 sweep 125개 항목의 상태는 `requested_sweep_status.json`에 있다.

Bound-state spectral residue와 continuum+bound 합은 독립 소형 문제에서 검증했다. 실제 소자에서 발견한 pole의 물리적 기원, 모든 dark state, 비평형 점유는 완결되지 않았다. Δx pole의 energy/localization/Gamma 및 mode/mesh/buffer 변화를 원자료로 제공하며 이를 임의로 점유하지 않았다.

## 시간·메모리 및 재현

최종 비평형 시도는 모델 판별 및 정밀화 관문 미통과 후 중단했다. Raised는 최소 1390.8초, Uniform은 최소 1243.2초 실행했으나 첫 비평형 SCF iteration이 완료되지 않았다. 각각 관측 peak RSS 330.6 / 255.3 MiB. 임의 TIMEOUT이나 물리적 해 부재로 표시하지 않았으며 `interrupted_attempts.csv`와 interruption 원기록에 구분했다.
측정 자원은 CPU quota 8 cores, RAM limit 8 GiB. 실제 시간과 peak RSS는 `performance.csv`, 단계별 비용은 profile 로그에 있다. 이는 이 환경의 측정값이며 다른 시스템의 walltime 보장이 아니다. 첫 고차 적분 시도의 메모리 결함과 보수적인 국소 오차 배분을 수정했다. 중단된 시도는 별도 interruption/resource_stop 기록으로 남겼으며 성공 또는 TIMEOUT 없는 수렴으로 꾸미지 않았다.

환경 설치, 실행, checkpoint, 모드 refinement, 전체 목표 sweep 명령은 `RUN_INSTRUCTIONS.md`를 따른다. 모든 필수 항목이 PASS일 때만 certificate가 RELEASE_READY를 true로 계산한다. 현재 결과에서 SS/Ioff/Vth/gm 또는 1% 전류 정확도를 보장하지 않는다.

AI/DFT/phonon/DD는 OUT_OF_SCOPE이다. GitHub push와 학교 HPC 제출을 하지 않았고 SIESTA/private 원본은 ZIP에 포함하지 않았다.
