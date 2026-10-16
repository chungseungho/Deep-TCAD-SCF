#!/bin/bash
#SBATCH --job-name=GAA-STEP05
#SBATCH --partition=X2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
# STEP04 uses X2 and ../../.venv; preserve the actual repository server convention.
source "${TCAD_VENV:-../../.venv}/bin/activate"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export GIT_COMMIT="$(git rev-parse HEAD)"
python -u run.py "$@"
