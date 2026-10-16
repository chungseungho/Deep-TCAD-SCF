# STEP08 — Contact and reservoir physics

**전극 길이를 늘리는 것과 hard wall 옆 reservoir가 bulk-like해지는 것은 다르다.** 같은 uniform wire의 중앙은 열린 전극 밀도로 수렴하지만, 경계 바로 옆의 고정된 구간은 depletion을 유지한다.

| L (nm) | center relative error | outer 3 nm neutrality error |
| --- | --- | --- |
| 12.0 | 0.0531336 | -60% |
| 24.0 | 4.88431e-06 | -60% |
| 48.0 | 7.62455e-11 | -60% |
| 96.0 | 1.23235e-14 | -60% |

이 실험은 fixed uniform potential이며 gradient=0은 입력 조건이다. SCF screening을 증명한 값이 아니다. L=12/24/48/96 nm와 mode=4/8/12/16/20/32/48/64/100/300, fixed mu/fixed N을 모두 계산했다. LDOS는 `-df/dE`로 thermal convolution한 delta spectrum이다. GF에 인위적인 eta를 넣지 않았다. Transverse spacing 1/0.5/0.25 nm에서 full-basis bulk density와 neutrality mu도 비교했다.

STEP07의 약 -48% reservoir charge 오차를 단일 원인으로 단정할 수 없다. 불충분한 occupied modes, 끝 hard wall, finite-N mu와 neutral-contact mu 차이, 충분하지 않은 screening이 함께 작용한다. 전체 N을 맞춰도 지역별 전하가 맞는 것은 아니다. 또한 nodal donor quadrature 72.840918878과 quantum-cell donor 72.6의 차이는 서로 다른 source discretization의 오차이며, 이를 fitting으로 renormalize하지 않았다.

Finite WBL path C는 4 electrode lengths × 4 gamma 값을 독립 small chain에서 비교했다. L=64의 current error는 gamma .0025/.005/.01/.02 eV에서 각각 47.62%, 18.87%, 2.757%, 0.7703%. 특정 gamma의 일치를 열린 전극의 물리적 인증으로 삼지 않는다. 채널 density 오차도 동시에 기록했다.

Energy-dependent embedding 대안은 positive-residue finite-chain pole expansion을 평가했다. 8/16/32/64/128 poles의 complex-axis self-energy 오차는 `results/uniform/summary.json`에 있다. 이 결과는 equilibrium contour representation의 가능성만 보여 준다. Biased GF-free density, initial preparation, long-time/continuum limit을 해결한 구현은 아직 없다.

![boundary](results/figures/reservoir_boundary.png)
![finite bath](results/figures/finite_bath_convergence.png)

실제 GAA electrode length × mode × mu의 fully coupled plateau는 미달성이다. Uniform/reduced-wire 결과로 이를 대신 주장하지 않는다.
