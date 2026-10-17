> Historical baseline document retained for research provenance. Its results and paths are not STEP11 validation evidence. See RESULTS.md and RELEASE_CERTIFICATE.json for this candidate.

# STEP08 — Validation status

**Unit/native tests 64 passed; GAA physical certification NOT passed.**

| 검증 | 실행 범위 | 판정 |
|---|---|---|
| Quantum/numerical tests | 60 tests | pass |
| 실제 DEVSIM/Gmsh | 4 tests, 별도 process | pass |
| Number bounds / degenerate tracking / invalid trial recovery | 신규 tests | pass |
| Batch-size invariance | 1 vs 96 poles, 1 vs 200 energies | pass |
| Uniform wire | 4 lengths × 10 mode levels × 2 number prescriptions | completed; 내부 수렴 확인 |
| Small constriction | 3 dx × 3 lengths × 5 bases | full-basis identity pass |
| Actual GAA frozen field | M2–M36, two number prescriptions | 전체 plateau 미달성 |
| Same-mu frozen GAA | M2–M20 | closed/open density 차이 남음 |
| Actual MS ambiguous state | 186 records + 3 modes × 7 frozen biases | origin unresolved |
| Open chain density/current | Vd=0/.001/.01/.05/.1 × 3 energy refinements | independent dense reference pass |
| GAA domain Newton | valid half-step recovery 기록 | interrupted, biased solution 미인증 |
| GAA equilibrium / C reference | warm-start 재현 | numerical convergence only |
| Grid/contact coupled GAA convergence | 계획 생성 | 미실행 |
| 정상 GAA I–V / AI teacher | prerequisites failed | 미제출 / 사용 불가 |

Open chain independent density max absolute error=1.50819e-09, current absolute error=1.09184e-13 A per degeneracy. Current conservation과 full dense rho PSD도 확인했다. 이 작은 시스템 검증을 GAA certification으로 확장하지 않는다.

Test logs는 `results/tests_batched/`, 이전 core version 검증은 `results/test_report/`에 있다. Native/quantum 그룹은 분리 실행하며 JUnit completion XML을 요구한다. 예전 monolithic native 종료 문제를 성공으로 세지 않았다.

Final production digest: `fe3c540c88e57070da79b59a8083e2a0bbf947e9d7faea752dac402fb6ff5fdc`. Snapshot은 `results/source_snapshots/`. 원본 source snapshot, 초기 seed, source/config identities와 실패 trace를 포함했다. 검증 script 자체도 최종 package manifest로 고정한다.

![chain](results/figures/open_chain_reference.png)
