#!/usr/bin/env bash
# User submits this script. No partition/account/resource entitlement is assumed.
# Example: sbatch --partition=<yours> --cpus-per-task=4 --mem=8G --time=24:00:00 job.sh stage-a-ms
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export PYTHONDONTWRITEBYTECODE=1
TASK="${1:-test}"
PYTHON_EXE="${STEP07_PYTHON:-$PWD/.venv/bin/python}"
case "$TASK" in
  test) "$PYTHON_EXE" validation/run_tests.py ;;
  stage-a-ms) "$PYTHON_EXE" run.py idvgs --config configs/gaa_ms.json --out results/user_ms_transfer --vd .1
              "$PYTHON_EXE" run.py idvds --config configs/gaa_ms.json --out results/user_ms_output ;;
  stage-a-negf) "$PYTHON_EXE" run.py idvgs --config configs/gaa_negf.json --out results/user_negf_transfer --vd .1
                "$PYTHON_EXE" run.py idvds --config configs/gaa_negf.json --out results/user_negf_output ;;
  stress-gate) "$PYTHON_EXE" run.py stress --config configs/gaa_ms.json --out results/user_gate_stress --axis gate --gate-range 0:6:.25 --vd .1 ;;
  stress-drain) "$PYTHON_EXE" run.py stress --config configs/gaa_ms.json --out results/user_drain_stress --axis drain --drain-range 0:6:.25 --vg .8 ;;
  *) echo "Unknown task: $TASK" >&2; exit 2 ;;
esac
