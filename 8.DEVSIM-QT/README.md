# DEVSIM-QT — STEP08

3D GAA quantum transport의 root-cause audit와 독립 검증 코드.

**이번 단계에서 정량 검증된 비평형 GAA I–V는 완성되지 않았다.** 64개 테스트와 작은 열린 시스템의 독립 density/current 검증을 통과했다. 실제 GAA는 36/72 modes까지 검사했지만 전체 density/contact plateau가 아직 없다. 이 데이터를 AI Teacher로 사용하면 안 된다.

수정한 핵심은 미분류 점유의 전자수 bounds, overlap/subspace tracking, 물리적으로 유효한 trial만 받는 Newton backtracking, current-change acceptance, GF energy/pole batching이다. 예외 무시, 임의 broadening, mass/doping fitting, threshold 완화는 하지 않았다.

| 읽을 파일 | 내용 |
|---|---|
| `ROOT_CAUSE_REPORT.md` | 구현 오류 / 수치 오차 / 모델 한계와 최종 6개 질문 |
| `MSDFT_SOURCE_COMPARISON.md` | 복원한 SIESTA 원본 호출 관계, 식, 줄 위치 |
| `MODE_CONVERGENCE.md` | M1–M5 및 36/72까지 실제 결과, fixed-mu/fixed-N |
| `CONTACT_CONVERGENCE.md` | uniform finite/open, length/mode/contact 실험 |
| `MS_NEGF_COMPARISON.md` | region별 같은-mu density 차이, 미분류 상태 |
| `VALIDATION.md` | 실제 통과한 것과 미완료 항목 |
| `HOW_TO_RUN.md` | 재현 명령 |
| `HPC_PLAN.md` | 전제 조건을 통과한 뒤 수행할 최소 순서 |
| `STATUS.json` | 기계가 읽을 수 있는 실패/인증 상태 |

A=`ms`: closed Multi-Space-inspired EMA/Hartree density, 유효한 경우에만 사용.
B=`negf`: semi-infinite lead self-energy를 쓰는 GF density, independent reference 경로.
C=`spectral`: finite WBL bath analytic density, 별개 접촉 모델.

`results/`에는 완료/거부/중단/OOM 실행을 모두 구분해 보관했다. 그림은 `results/figures/`. `baseline/`은 STEP07 production snapshot과 warm-start seed다. 초기 seed를 쓴 재현을 cold-start 검증으로 부르지 않는다. 최종 provenance는 `PACKAGE_MANIFEST.json`, core snapshots는 `results/source_snapshots/`.

최신 `siesta-5.4.2.tar(1).gz`는 전달되지 않았다. 이전에 제공된 다른 첨부를 복원해 분석했으며 최신 파일과 동일하다고 주장하지 않는다. 연구실 Fortran 원문이나 압축파일은 배포물에 없다. GitHub 변경과 학교 서버 접속/작업 제출은 하지 않았다.
