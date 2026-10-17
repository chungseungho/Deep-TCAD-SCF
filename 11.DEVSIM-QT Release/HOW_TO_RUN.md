# STEP08 — Reproduction

모든 명령은 `8.DEVSIM-QT/` 안에서 실행한다. 다른 STEP 폴더는 필요하지 않다. Python 3.12.14, requirements의 고정 버전으로 검사했다. Gmsh의 GLU와 DEVSIM의 BLAS 등 native system libraries는 실행 환경에 있어야 한다. 시스템 환경을 덮어쓰지 않고 전용 환경을 사용한다.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
.venv/bin/python validation/run_tests.py results/user_tests
```

DEVSIM 설치가 BLAS library 탐색을 요구하면 해당 시스템에서 사용할 실제 BLAS에 맞춰 `DEVSIM_MATH_LIBS`와 library search path를 설정한다. 패키지에 특정 학교 서버 경로나 shared environment 변경 명령은 넣지 않았다.

## 빠른 독립 검증

```bash
.venv/bin/python validation/step08_uniform.py results/user_uniform
.venv/bin/python validation/step08_constriction.py results/user_constriction
.venv/bin/python validation/step08_open_chain.py results/user_open_chain
```

Uniform 실험은 mode/length/fixed-mu/fixed-N을 분리한다. Constriction은 independent full 3D stencil과 dense inversion을 비교한다. Open-chain은 Vd=0,.001,.01,.05,.1 V에서 dense GF integration과 생산용 contour/injectivity 계산을 비교한다. 이 chain 결과는 GAA I–V가 아니다.

## GAA frozen mode study

```bash
.venv/bin/python validation/step08_modes.py baseline/ms_seed.npz results/user_modes 2,4,6,8,10,14,20
.venv/bin/python validation/step08_modes.py baseline/ms_seed.npz results/user_modes_high 20,28,36
.venv/bin/python validation/step08_modes.py baseline/ms_fixed_mu_seed.npz results/user_matched_mu 2,4,6,8,10,14,20
.venv/bin/python validation/step08_origin_sensitivity.py baseline/ms_fixed_mu_seed.npz results/user_origin
```

`ms_seed.npz`는 finite-N STEP07 field(mu≈−.057578 eV), `ms_fixed_mu_seed.npz`는 fixed-contact-mu field(mu≈+.009559 eV)다. 비교에서 바꿔 쓰면 안 된다. 출력에 요청/실제 dimensions, mu, N, 지역별 charge, density, T(E)가 있다. High-mode 실험은 GF batching을 써도 dense eigenproblem 비용이 크다.

## 실제 low-mode GAA 재현

```bash
.venv/bin/python run.py idvds --config configs/step08_ms.json --gates .8 --drain-range 0,.01,.05,.1 --initial-state baseline/ms_seed.npz --out results/user_ms
.venv/bin/python run.py idvds --config configs/step08_negf.json --gates .8 --drain-range 0,.001,.01,.05,.1 --initial-state baseline/negf_seed.npz --out results/user_negf
.venv/bin/python run.py point --config baseline/spectral_config.json --vg .8 --vd .1 --initial-state baseline/spectral_seed.npz --out results/user_spectral
```

MS bias rejection과 NEGF bound rejection은 예상되는 검증 결과이며 exit 2를 무시해 성공 처리하지 않는다. NEGF nonlinear search는 길게 걸릴 수 있다. `*_newton_trials.jsonl`에서 invalid trial / accepted step / physical image를 구분한다. `trials/*.npz`의 유효한 trial도 `numerically_converged=false`이며 학습용 해가 아니다. 미분류 점유의 `*_number_interval.json`과 valley별 NPZ를 확인한다.

각 출력 폴더는 config/source identity가 고정된다. source나 config를 바꾸면 새 폴더를 사용한다. `--resume`은 동일 identity의 수렴한 state를 재사용한다. 결과 CSV의 NaN은 실패이며 보간하지 않는다.

## I–V와 0–6 V 지원

다음 명령은 구현되어 있지만 **이 패키지에는 아직 certified.json이라는 검증된 configuration이 없다.** 검증 완료 후 사용자가 그 설정을 저장해 사용한다.

```bash
.venv/bin/python run.py idvgs --config certified.json --vd .1 --gate-range 0:1.2:.05 --out results/user_transfer
.venv/bin/python run.py idvds --config certified.json --gates .4,.8,1.2 --drain-range 0:1.2:.05 --out results/user_output
.venv/bin/python run.py stress --config certified.json --axis drain --vg .8 --drain-range 0:6:.05 --out results/user_stress
```

0–6 V parsing은 `--dry-run`으로 확인할 수 있다. 실제 high-field numerical completion과 EMA physical validity는 별도다. Oxide breakdown, impact ionization, scattering/heating, nonparabolicity는 현재 모델에 없다. 기존 `resonance_order`, `resonance_max_iterations` 필드는 legacy 호환용이며 현재 integration 제어에 연결되어 있지 않다. 실제 제어는 `energy_*`, `fermi_poles`, `contour_batch_size`, `transmission_batch_size`다.

## 보고서와 package

```bash
.venv/bin/python validation/step08_reports.py
.venv/bin/python validation/step08_package_audit.py
```

보고서 generator는 동봉된 전체 결과 폴더의 존재를 전제로 한다. 소스-only 복사에서 모든 연구 실험을 재실행하기 전에는 보고서를 새로 생성하지 않는다.
