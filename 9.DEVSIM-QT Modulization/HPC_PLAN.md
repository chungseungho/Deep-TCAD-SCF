# STEP08 — 최소 후속 계산 순서

**17개 설정을 생성했지만 학교 서버 실행/제출은 0건이다.** `hpc_plan/plan.json`의 모든 executed 값은 false다. 먼저 다음 관문을 통과한 뒤 다음 계산을 시작한다. 100-point sweep부터 제출하지 않는다.

1. **Frozen mode gate:** 같은 field/동일 mu에서 36/72, 48/96, 60/120을 비교한다. 직교 기저 full limit는 현재 transverse grid에서 channel 60, S/D 300이다. 관심 energy window의 occupied density/region charge/LDOS, adaptive T/current가 세 연속 단계에서 안정되는지 확인한다. 채널 charge 하나의 <2%는 충분하지 않다. 필요하면 S/D modes를 더 늘려야 한다. 현재 36/72도 전체 density 2.98%, sampled T 5.66% 변화가 남는다.
2. **Equilibrium contact gate:** 모드가 안정된 수준에서 source/drain length=12/24/36 nm를 사용한다. 먼저 동일한 contact mu를 고정한 비교와, full lead DOS로 중성을 맞춘 비교를 구분한다. 자동 calibration은 retained basis마다 mu를 재계산하므로 이를 같은-mu 비교라고 부르면 안 된다. 모든 경우 outer wall 바로 옆이 아니라 접합/끝에서 떨어진 bulk-like plateau를 찾아 projector 경계를 정한다. 서로 다른 projector 경계에서의 안정성도 검사한다.
3. **Open bound-state gate:** reference equilibrium을 확보한 뒤 Vd=0,.001,.01,.05,.1 V 순으로 진행한다. 새 유효 seed에서 계속하고, domain trial 실패와 final root failure를 구분한다. 최종 전위에 isolated state가 있으면 preparation history/relaxation 모델을 먼저 명시하고 검증해야 한다. 현재 closed-origin 또는 임의 평균 filling으로 대체할 수 없다. Embedded dark-state audit도 필요하다.
4. **Joint grids:** dx=.5/.25/.125 nm; dy=dz=1/.5/.25 nm; 세 Poisson meshes를 순차 비교한다. Transverse grid를 바꿀 때 modes를 자동으로 충분하게 늘려야 한다. Mode 기준을 하나도 확인하지 않은 채 grid만 8배 키우면 비용만 증가한다. 이 단계부터 regional n, barrier, current, Gauss, current continuity와 stricter SCF/energy tolerances를 함께 비교한다.
5. **I–V:** 위 관문을 통과한 한 경로와 설정을 `certified.json`에 기록한 뒤 transfer 25개, output 75개를 실행한다. 성공한 bias만 표시하고 실패를 보간하지 않는다. 0–6 V는 물리 적용 범위를 별도 기록하는 stress experiment로 다룬다.

먼저 실행할 명령 예:

```bash
./job.sh test results/server_tests
./job.sh modes baseline/ms_fixed_mu_seed.npz results/server_frozen 36,48,60
./job.sh point --config hpc_plan/negf_m36_sd72_L12.json --vg .8 --vd 0 --out results/server_eq_m36
```

SLURM 제출은 사용자가 자신의 partition/account와 확보한 memory/time으로 수행한다.

```bash
sbatch --partition=<사용자_partition> --cpus-per-task=<CPU수> --mem=<메모리> --time=<시간> job.sh modes baseline/ms_fixed_mu_seed.npz results/server_frozen 36,48,60
```

`STEP08_PYTHON`으로 전용 Python 실행 파일을 지정할 수 있다. job.sh 자체에는 자동 sbatch, ssh, partition/account 추정이 없다.

## 비용과 필요한 구현 개선

이 구현은 block RGF를 사용하지만 finite H/eigenvectors와 bound audit 일부를 dense로 유지한다. 8 GiB 환경의 unbatched M28에서 OOM이 확인됐고, energy batching으로 M36은 완료됐다. 다음 단계에서는 valley를 순차 처리하고 occupied eigenpairs/contour density를 활용하여 full eigensystem 메모리를 줄이는 것이 좋다. Full GAA basis에서는 한 real dense matrix만 수 GiB이므로 이론적 최소 배열 크기보다 여유 있는 자원이 필요하다. 실행 전 dimension과 배열 수로 memory를 산정한다. 64 GiB와 같은 숫자를 서버 권한이나 충분성 보증으로 가정하지 않는다.

비평형 adaptive quadrature는 energy-node별 matrix output과 panel error 배열도 크다. Energy batching만으로 모든 메모리 병목이 해결된 것은 아니다. 전역 componentwise error budget을 유지하는 streaming accumulation을 추가 검증해야 한다. Runtime budget 때문에 중단된 STEP08 GAA Newton run은 완료로 세지 않는다.

## 서버에서 새로 필요한 최소 물리 판단

큰 서버만 사용한다고 injection origin 또는 bound-state preparation이 자동 결정되지는 않는다. 충분한 mode/contact에도 ambiguity가 남으면 A의 적용 영역을 제한한다. B의 bound state는 equilibrium에서 시작한 adiabatic preparation 등 **구체적인 물리 모델과 상태 생성/소멸 규칙**을 먼저 정해야 한다. 이 패키지는 그 물리를 임의로 골라 density를 반환하지 않는다.
