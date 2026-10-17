#!/usr/bin/env bash
# STEP08 DEVSIM-QT. The user supplies partition/account/time/memory to sbatch.
# No job is submitted by this script, and no cluster entitlement is assumed.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export PYTHONDONTWRITEBYTECODE=1
TASK="${1:-test}"
if (($#)); then shift; fi
PYTHON_EXE="${STEP08_PYTHON:-$PWD/.venv/bin/python}"
case "$TASK" in
  test) exec "$PYTHON_EXE" validation/run_tests.py "${1:-results/user_tests}" ;;
  uniform) exec "$PYTHON_EXE" validation/step08_uniform.py "${1:-results/user_uniform}" ;;
  constriction) exec "$PYTHON_EXE" validation/step08_constriction.py "${1:-results/user_constriction}" ;;
  modes) exec "$PYTHON_EXE" validation/step08_modes.py "$@" ;;
  chain) exec "$PYTHON_EXE" validation/step08_open_chain.py "${1:-results/user_chain}" ;;
  plan) exec "$PYTHON_EXE" validation/step08_plan.py "${1:-hpc_plan_user}" ;;
  point|idvgs|idvds|stress) exec "$PYTHON_EXE" run.py "$TASK" "$@" ;;
  *) echo "Tasks: test uniform constriction modes chain plan point idvgs idvds stress" >&2; exit 2 ;;
esac
