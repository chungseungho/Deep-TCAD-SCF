# STEP08 — MS versus open NEGF

**동일 mu의 실제 GAA density 차이는 해결되지 않았다.** 원래 STEP07의 independently coupled equilibrium 차이 24.0657%와 frozen-field 차이 24.7878%는 서로 다른 비교다. 이번에는 저장된 fixed-mu MS 평형장에 두 density solver를 적용해 frozen 24.7878%를 재현했다. M20/40에서도 17.82%가 남는다.

| Region | M2 closed/open difference | M20 closed/open difference |
| --- | --- | --- |
| channel | 6.244e-06% | 2.969e-06% |
| extension | 0.176% | 0.1137% |
| left | 42.6% | 33.36% |
| right | 42.53% | 32.98% |

1. Uniform wire: 같은 H dispersion, full modes, T, mu에서 중앙 밀도의 finite-length limit을 확인했다. 경계부 전체를 같다고 주장하지 않는다.
2. Small constriction: 같은 potential, full real-space reference에서 projection/density/T/current identities를 확인했다. Length와 dx를 분리한 45 case를 계산했다.
3. Actual frozen GAA: fixed-N field와 fixed-contact-mu field를 별도로 사용했다. mode sweep에서 각각 fixed mu와 fixed N을 분리했다.
4. Fully coupled GAA: STEP07의 낮은 모드 평형들을 warm-start 재현했다. 고차 모드의 독립 SCF 비교와 비평형 인증은 미완료다.

## Occupation and tracking

실패 파일의 모든 valley를 저장하도록 바꿨고, 0.01/0.05/0.1 V에서 active ambiguous state 186건의 energies, weights, profile, localization, nearest gap, reach, contrast, boundary stability를 추출했다. `state_audit/active_states.csv`와 NPZ가 근거다. 이것은 실패한 iterate의 상태이지 SCF 해가 아니다.

추가로 `phi(vd)=phi_eq+vd*x/L`의 **명시적 frozen family**에서 M2/M4/M10과 7 bias를 계산했다. 이 controlled family는 actual Poisson bias continuation이 아니다. overlap tracking과 degenerate projector 비교를 수행해도 origin은 별도 미결정으로 남긴다. 

대칭 uniform chain은 active 19개가 모두 ambiguous인데 open T=1이다. 이는 closed standing-state probability만으로 일반적인 injection origin을 결정할 수 없다는 반례다. 원래 MS-DFT 전체 이론이 틀렸다는 주장은 아니다.

![comparison](results/figures/matched_mu_density.png)
![states](results/figures/ambiguous_states.png)
