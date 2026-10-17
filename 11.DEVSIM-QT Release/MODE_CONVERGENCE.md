> Historical baseline document retained for research provenance. Its results and paths are not STEP11 validation evidence. See RESULTS.md and RELEASE_CERTIFICATE.json for this candidate.

# STEP08 — Mode convergence

**GAA 전체 수렴은 미달성이다.** M1–M5를 넘어 14/28, 20/40, 28/56, 36/72까지 실제 frozen-field 계산을 수행했다. 마지막 단계에서도 전체 density L2 변화는 2.98%, transmission sampled L2 변화는 5.66%다. 채널 charge만 안정되는 것을 전체 검증으로 취급하지 않았다.

## Finite-N STEP07 평형장을 고정한 결과

모든 비교는 같은 xyz 격자, 동일한 frozen phi, 같은 온도다. fixed-mu 열의 mu=-0.057577610836 eV는 해당 finite-N 평형장의 저장값이며 neutral open-lead mu와 다르다. fixed-N은 72.6 electrons를 유지하도록 mu를 다시 푼다.

| Channel / S-D | N (fixed mu) | mu (fixed N, eV) | density refinement (fixed mu) | channel charge refinement (fixed mu) |
| --- | --- | --- | --- | --- |
| 2/4 | 72.6 | -0.0575776 | — | — |
| 4/8 | 117.69 | -0.0737514 | 43.49% | 32.31% |
| 6/12 | 147.186 | -0.0810991 | 29.54% | 11.84% |
| 8/16 | 162.529 | -0.0843066 | 14.83% | 3.792% |
| 10/20 | 173.015 | -0.0863123 | 11.27% | 2.578% |
| 14/28 | 185.939 | -0.088569 | 13.3% | 2.688% |
| 20/40 | 193.983 | -0.0898937 | 9.218% | 1.052% |
| 28/56 | 197.85 | -0.0904736 | 4.92% | 0.3894% |
| 36/72 | 199.711 | -0.0907419 | 2.98% | 0.05398% |

실제 S/D rank는 요청한 sd_modes에 narrow trace를 보강하므로 더 클 수 있다. 정확한 valley별 dimension은 각 summary.json에 있다. 고정 mu에서 N이 증가하는 것 자체는 오류가 아니다. 같은 H의 더 넓은 공간을 채우는 결과이고, 별도로 mu를 움직이는 fixed-N 비교와 구분해야 한다.

`linear_response_conductance_fixed_contact_mu_s`는 저장한 221개 T(E) sample의 thermal integral 진단값이다. 좁은 resonance에 대한 독립 adaptive 검증을 통과한 비평형 전류가 아니며, fixed-N 행에도 비교용 같은 contact mu의 값이 들어 있다. 이를 실제 GAA I–V로 사용하지 않는다. Lowest local subbands, occupied eigenstate counts, regional charges, T(E), closed/open densities도 저장했다. 공간 전체 LDOS(E,r)와 충분한 GAA occupied-subband plateau는 아직 인증하지 않았다.

## Uniform contact DOS

같은 full-basis neutral mu=-0.00962398213094 eV에서 4 modes는 45.22%, 20 modes는 4.59%의 bulk line charge를 놓친다. 64/100/300 modes는 각각 0.05782%, 0.0004437%, 0%. 이것은 rectangular 20x15 nm uniform lead의 결과이며 GAA의 contact-mode 기준을 자동 결정하지 않는다.

## Independent real-space limit

`results/constriction`은 dx=1,0.5,0.25 nm × electrode length=2,4,8 nm × 5 bases = 45 cases이다. 독립 3D stencil로 만든 H와 `U†HU`, 고유값, closed density, reconstruction, dense real-space open inversion의 T/current를 비교했다. Full-basis identity 통과=True. Truncated-basis 오차는 삭제하지 않았다. Nested trace 보강은 projection completeness와 occupied-density convergence를 같은 것으로 만들지 않는다.

## 메모리 수정

Unbatched M28 실험은 8 GiB cgroup에서 OOM(exit -9)로 중단됐다. 결과를 성공으로 처리하지 않았다. Contour pole 8개, transmission energy 16개씩 처리하도록 수정한 뒤 M28과 M36 계산을 완료했다. Pole 수, self-energy, tolerance는 그대로다. Batch 1과 batch 96/200 비교를 테스트했다. 고차 GAA 전체 고유벡터와 real-axis adaptive panel 배열은 여전히 비용이 크다.

![mode](results/figures/mode_refinement.png)

Raw: `results/gaa_modes`, `results/gaa_modes_batched`, `results/matched_mu_modes`, `results/mode_summary.csv`. `gaa_modes_high`는 실패한 unbatched 실행이다. GAA dx/transverse/Poisson joint convergence는 미실행이며 `HPC_PLAN.md`를 따른다.
