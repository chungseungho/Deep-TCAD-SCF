#!/bin/bash
#SBATCH --job-name=GAA-STEP06
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=step06-%j.out
#SBATCH --error=step06-%j.err
# Select partition, CPU, memory and walltime on sbatch's command line according
# to your school's policy; no partition/resource entitlement is presumed.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "${TCAD_VENV:-.venv}/bin/activate"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1
python -u run.py "$@"
