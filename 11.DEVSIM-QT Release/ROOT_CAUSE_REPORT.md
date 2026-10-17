> Historical baseline document retained for research provenance. Its results and paths are not STEP11 validation evidence. See RESULTS.md and RELEASE_CERTIFICATE.json for this candidate.

# STEP08 — Root-cause report

프로젝트: **DEVSIM-QT**. 이 패키지는 원인 분석과 검증 가능한 수정을 포함하지만, **정량 검증된 비평형 GAA I–V / Fully-SCF Teacher는 아직 완성되지 않았다.**

## 원인과 실제 수정

| 문제 | 분류 | 근거 | 이번 조치 |
|---|---|---|---|
| unknown origin을 number equation에서 L로 계산 | Algorithmic error | STEP07 `multispace/constraints.py`의 `origin==2` 외 모두 fL | N_low/high, mu interval을 먼저 계산하고 material ambiguity는 density 생성 전에 거부 |
| index 기반 상태 비교 | Algorithmic weakness | eigenvalue crossing/degenerate rotation 가능 | Hungarian overlap matching, degenerate-subspace singular values; origin 결정과 분리 |
| invalid Newton trial이 전체 SCF를 즉시 abort | Algorithmic error | STEP07 failure trace, 실제 GAA bound audit | 유효한 domain에서만 Armijo 감소 스텝 허용; half-step 성공, invalid derivative probe와 stalled search 기록 |
| occupied transverse modes 부족 | Numerical approximation | uniform 4-mode missing charge 45.22%; actual GAA 36-mode density도 변화 | M1–M5 및 M14/M20/M28/M36, fixed-mu/fixed-N 분리 실험 |
| finite wall와 열린 전극 차이 | Boundary/model + finite-size approximation | 중앙은 길이에 수렴하나 outer 3 nm depletion 유지 | uniform mode/length/mu 결합 실험; region별 density comparison |
| standing eigenstate의 injection origin 부재 | Physical applicability limit | symmetric wire: T=1인데 active 19개 모두 ambiguous | threshold 완화나 임의 filling 없이 거부; open path 별도 유지 |
| isolated bound occupation | Physical model incompleteness | Vg=.8,Vd=.001 initial trial의 pole E=0.0052550675 eV, fL−fR=0.00961798 | audit 보존; final solution에서도 유효해야 함. 준비/relaxation 모델 임의 추가하지 않음 |
| GF pole/energy 배열의 과도한 메모리 | Implementation resource problem | 8 GiB에서 M28 unbatched OOM | energy/pole batching; 동일 수치 확인, M36까지 완료 |
| current tolerance가 acceptance에 미사용 | Algorithmic omission | STEP07 acceptance 조건 | A의 후처리는 별도, B/C는 연속 current change도 기존 tolerance로 확인 |

## 수정 전후 재현

Original STEP07 source digest `c0fbe6dc...`를 별도 복원해 실행했다. 이전 density를 초기값으로 사용한 **warm-start 재현**이며 independent cold-start 수렴이라고 주장하지 않는다. Vg=.8 equilibrium MS/NEGF가 다시 수렴했고, MS Vd=.01/.05/.1은 다시 거부됐다. .05의 첫 valley에서 기존 11개 ambiguous state도 재현했다. Finite spectral Vg=.8,Vd=.1은 약 4.49181 uA로 재현됐다.

전자수 수정 후 .05 V의 feasible mu interval은 [-0.04742291141,-0.04686515073] eV이고 audit mu에서 전체 N 범위는 [72.04560806,73.16232108] electrons이다. 1.11671302 electrons의 불확실성을 먼저 드러낸다. 이를 제거했다고 주장하거나 arbitrary L로 맞추지 않는다.

NEGF domain Newton의 실제 GAA 로그에서는 invalid full step 뒤 valid half-step이 채택됐다. 두 audit run은 계산 비용 때문에 중단했고 `INTERRUPTED.json`으로 표시했다. 이전의 즉시-abort 구현 오류는 수정됐지만, 그 실행들을 최종 fixed point 성공/실패로 판정하지 않았다. .001/.01의 초기 trial에서 pole이 발견된 경우 역시 모든 해의 비존재 증명이 아니다.

## 최종 연구 질문에 대한 답

1. **실패는 복합 원인이다.** 구현 오류, 부족한 mode/contact convergence, closed-origin approximation의 물리적 한계를 각각 확인했다.
2. **Uniform wire 내부는 수렴한다. 실제 GAA 전체는 미인증이다.** finite boundary와 open reservoir를 무조건 동일시하면 안 된다.
3. **검사한 GAA 비평형 MS 점유는 안정적으로 정의되지 않았다.** 모드 증가와 tracking만으로 해결되지 않는다.
4. **정당한 이식 범위는 projector/number/gauge/density bookkeeping이다.** SIESTA의 orbital averaging/tie-break/XC/variational constrained search까지 EMA가 구현했다고 할 수 없다.
5. **Open NEGF는 small-system 독립 reference를 통과한다. GAA reference 완성은 남았다.** bound preparation, embedded dark states, joint convergence가 아직 필요하다. C는 별개의 finite bath model이다.
6. **검증된 정상 GAA I–V는 이번에도 제출할 수 없다.** 전제인 density/mode/contact plateau가 통과하지 않았으므로 100개 bias sweep으로 이를 숨기지 않았다. 지원 코드와 조건부 실행 명령은 제공한다.

## 남은 구체적 한계

Bound audit는 모든 propagating lead band 아래의 poles를 검사한다. Embedded dark states의 전역적인 배제, above-band occupied poles, 전체 GAA density-matrix PSD는 완전히 검증하지 않았다. Finite real-space quadrature의 매우 좁은 resonance는 별도 integration budget 문제다. 원래 donor/charge source discretization 차이와 GAA mesh refinement도 미해결이다. 기존 좁은 mode configuration에서 얻은 수렴값은 정량적인 reference가 아니다.

원본 SIESTA 새 첨부는 전달 실패했다. 이전 버전만 분석했으며 소스 비교 범위는 `MSDFT_SOURCE_COMPARISON.md`에 명시했다. 모든 단계가 완료됐다는 표현이나 AI training reference 인증은 하지 않는다.

![domain](results/figures/domain_newton_trials.png)
