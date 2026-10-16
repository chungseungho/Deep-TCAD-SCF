# 설치, 재현 및 HPC 실행

## 환경

Linux Python 3.12를 권장한다. 실제 계산은 Python 3.12.14에서 수행했으며, 테스트 환경은 `results/environment.json`과 `requirements-lock.txt`에 기록한다. 아래 `python3`가 해당 전용 Python을 가리키도록 환경/module을 먼저 선택한다.

```bash
cd 6.MS-modspace-model
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m pytest -q tests
```

Gmsh는 libXft.so.2 등 OS library가 필요할 수 있다. DEVSIM이 BLAS를 찾지 못하면 이미 설치된 OpenBLAS에 맞춰 `DEVSIM_MATH_LIBS=libopenblas.so.0`을 설정한다. 기존 학교 공용 환경을 변경하지 말고 작동하는 native library와 STEP06 전용 venv를 사용한다. 자동 OS 설치나 MPI/CUDA 변경은 없다.

## 단계별 실행

```bash
python run.py --phase mesh --smoke --out output/mesh
python run.py --phase poisson --smoke --out output/poisson
python run.py --phase point --smoke --dx .5 --modes 3 --sd-modes 6 --vg .8 --vd .1 --out output/point
python run.py --phase point --smoke --backend negf --dx .5 --modes 3 --sd-modes 6 --vg .8 --vd 0 --out output/equilibrium
```

`--smoke`의 dx=3 nm 그대로 수송 곡선을 해석하지 않는다. 위 point는 dx를 명시적으로 override한다.
실제 실행된 transfer/output의 provisional 설정은 dx=.25 nm, modes=2/4, dy=dz=1 nm, gamma=.01 eV, coarse Poisson mesh이다. 적은 mode 수와 contact 수렴은 검증 전이다.

```bash
python run.py --phase idvgs --smoke --dx .25 --modes 2 --sd-modes 4 --keep-going --out output/transfer
python run.py --phase idvds --smoke --dx .25 --modes 2 --sd-modes 4 --keep-going --out output/output
```

기본 sweep은 Vgs=0:1.2:.05 at Vds=.1, output은 Vgs=.4,.8,1.2에서 Vds=0:1.2:.05다. 상태별 status/dispersion/neutrality를 함께 확인한다. 같은 폴더에서 두 계산을 동시에 실행하지 않는다.

중단 후 동일 physics/numerics/source인 경우에만 `--resume`로 재개한다. 다른 설정은 새 output 폴더를 사용한다. 실패한 이전 bias는 더 작은 허용오차로 바꾸지 않으며, accepted 인접 상태를 initial guess로 사용해 재시도할 수 있다. `--initial-state path/to/state.npz`는 같은 Poisson nodes의 순수 수치 initial guess다. 최종 quantum–Poisson SCF를 생략하지 않는다.

`--resume`에서 해당 bias의 checkpoint가 있으면 그 checkpoint가 `--initial-state`보다 우선한다. 정체된 계산을 다른 accepted 초기값에서 다시 시작하려면 **새 output 폴더**에 `--initial-state`를 지정한다. 이전 실패/중간 기록은 그대로 보존한다.

```bash
python run.py --phase idvds --smoke --dx .25 --modes 2 --sd-modes 4 --keep-going --out output/output --resume
python run.py --phase plot --smoke --dx .25 --modes 2 --sd-modes 4 --out output/transfer
```

## 수렴성

```bash
python validation/convergence.py --groups mesh,dx,transverse,modes,length,gamma,basis,scf --out output/campaign
# plan.json / 각 config / 메모리 추정 검토 후 실제 local 실행
python validation/convergence.py --run --groups gamma,dx --out output/campaign
python validation/convergence.py --summary-only --out output/campaign
python validation/benchmarks.py --out output/benchmarks
```

기본 campaign bias는 Vgs=.8, Vds=.1이다. `--vg`, `--vd`로 equilibrium/off/threshold/on 및 .5/1 V를 각각 별도 폴더에서 실행한다. 세 단계 이상과 연속 refinement의 current/charge/barrier 차이를 검사한다. 한 bias, 한 축의 agreement는 전체 검증이 아니다. NEGF campaign에서는 gamma가 적용되지 않는다.

```bash
python validation/convergence.py --backend negf --groups mesh,dx,modes,length,energy,poles,basis,scf --out output/negf_plan
```

Semi-infinite 실제 비평형 GAA에서 adaptive integral이 에러를 내는 경우는 실패 결과다. `energy_max_nodes`, `energy_step`, `energy_order`, `energy_rtol`을 독립적으로 검증하고, 임의 broadening이나 느슨한 SCF tolerance로 바꾸지 않는다. 같은 finite-bath model의 post-SCF transmission은 `--spectra`로 계산할 수 있지만, sampled spectrum은 energy refinement 없이 독립 current 정확도 인증이 아니다.

실제 GAA finite-bath 전류의 별도 Schur/RGF 검증은 `python validation/post_scf_current.py results/device_lowbias vg0.800000_vd0.100000`으로 재현한다. `python validation/frozen_contact_comparison.py results/device_lowbias vg0.800000_vd0.100000`은 같은 저장 전위에서 semi-infinite 접촉으로 바꾼 **frozen-field 전류만** 비교한다. 이 결과는 별도 fully-SCF biased reference가 아니다.

## 자원과 SLURM

```bash
python run.py --dry-run --dx .25 --modes 6 --sd-modes 10
```

Finite spectral은 valley별 dense complex diagonalization이고 O(N^3)이다. 기본 fine preset은 valley당 1680차원이며 추가 refinement는 빠르게 비싸진다. 실제 한 bias 시간은 raw CSV에 기록되며 하드웨어가 바뀌면 다시 측정한다. `resources.json`의 메모리 추정에는 상수가 보수적으로 들어가지만 OS/Python 전체 peak를 보장하지 않는다. NEGF는 block/energy 비용과 매우 좁은 공명의 adaptive 비용을 추가로 고려한다.

`job.sh`에는 학교 자원 권한을 가정하는 partition/memory/time을 고정하지 않았다. 학교 정책을 확인해 `sbatch` 옵션으로 지정한다. 이전 X2/2CPU/32GB 사례가 현재 허용된다는 뜻은 아니다.

```bash
# 실제 partition 이름과 허용 자원을 직접 선택한다.
sbatch --partition=PARTITION --cpus-per-task=CPUS --mem=MEMORY --time=TIME job.sh --phase point --vg .8 --vd .1 --out output/hpc_point
```

이 작업에서 학교 SLURM 제출이나 원격 GitHub push는 수행하지 않았다. ZIP의 `6.MS-modspace-model/`을 저장소 옆에 복사한 뒤 사용자가 검토/commit/push한다.

## 접합 basis 보강 비교 실험

```bash
python validation/junction_overlap.py
python validation/junction_compatible.py --phase point --smoke --backend negf --dx .5 --modes 2 --sd-modes 4 --vg .8 --vd .1 --out output/nested_contact
```

보강 경로는 low S/D modes에 zero-extended channel modes를 추가하고 같은 실제 FD operator를 전부 재투영한다. 추가 rank, 확장 코드 SHA256과 configuration identity를 별도로 기록한다. 기존 I–V 결과와 섞지 않는다. 이 실험의 성공 여부는 해당 point.csv로 판단하며, 단위 검증이나 Gram overlap=I만으로 실제 수렴을 주장하지 않는다.

계산 종료 후 `python validation/enrich_states.py --root results`로 static material mesh와 mass/valley metadata, teacher eligibility=false catalog를 묶을 수 있다. 물리적 field 값은 바꾸지 않는다.

## 결과 정리 및 전달 무결성

분리 실행한 bias들은 source/config fingerprint가 같고 서로 중복되지 않을 때만 합친다. 각 partition solver를 먼저 종료한 뒤 `validation/merge_output_partitions.py --target TARGET --parts PART1 PART2 --require-complete`를 사용한다. 실제 CSV/NPZ를 복사할 뿐, 누락 전류를 보간하지 않는다. `partition_execution.json`에 원래 partition과 CSV hash를 남긴다.

최종 그림/보고서는 `python validation/report.py`, 저장 상태와 provenance 검사는 `python validation/package_audit.py`로 생성한다. Audit는 STEP01–05가 옆에 있을 때만 이전 파일의 checksum도 검사하며, STEP06 단독 추출에서는 해당 검사를 skip했다고 표시한다. 동일 config 결과를 다른 수렴 축의 공통 baseline으로 재사용한 경우에는 `reuse_identical_configuration.json`을 확인한다.

원시 계산의 core SHA256과 전달 코드의 SHA256이 다르면 `results/source_snapshots/index.json`을 확인한다. 실제 계산에 사용한 원본 core ZIP을 함께 보존했다. 최종 변경은 resume 이후 저장 전위의 native export 복원과 venv를 source hash에서 제외하는 처리다. Hamiltonian, 재료값, 점유 또는 SCF 기준을 바꾸지 않았으며 기존 NPZ의 source identity도 덮어쓰지 않았다. 원래 폴더를 `--resume`하려면 별도 복사본에 해당 core archive를 추출하고, 전달 코드를 사용할 때는 새 output 폴더를 선택한다.

원본 core를 사용할 때는 당시처럼 venv를 그 복사본 **밖에** 두고 활성화한다. 원본에는 이전 hash 규칙과 export 동작이 포함되며, 현재의 delivery regression 두 테스트는 수정된 전달 코드용이다. 일반적인 새 계산은 전달 코드와 새 output 폴더를 사용하는 방법을 권장한다.

미완료 checkpoint를 새 코드의 순수 초기값으로 사용할 수도 있다. 예를 들어 `python validation/checkpoint_seed.py results/convergence/dx_2 --out output/dx_seed.npz` 뒤 `python run.py --phase point --vg .8 --vd .1 --config results/convergence/dx_2.json --initial-state output/dx_seed.npz --out output/dx_fresh`로 다시 fully-SCF를 수행한다. 이 seed는 `converged=False`, `seed_only=True`이며 학습 가능한 물리 state가 아니다.

`python validation/reexport_saved.py --root results`는 저장된 accepted 입력 전하로 선형 Poisson만 복원해 native field를 내보낸다. 저장된 Si 전위와 1e-10 V 이내로 같은지 검사하고 `field_export.json`과 `field_export_audit.json`에 바이어스·좌표 단위·오차를 남긴다. Quantum SCF/current/NPZ를 새 결과로 대체하는 명령이 아니다.

모든 solver가 종료된 뒤 `python validation/make_delivery.py --out /absolute/path/Deep-TCAD-SCF_STEP06.zip`으로 source와 실제 결과, 전체 파일별 SHA256 manifest 및 CRC 검사를 포함한 ZIP을 만든다. 캐시, private venv와 수정 전 개발 pilot은 제외하고, 실제 실패한 reference와 미완료 checkpoint는 포함한다.
