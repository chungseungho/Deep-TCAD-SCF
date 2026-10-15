#!/bin/bash
#SBATCH --job-name=GAA-MSIV
#SBATCH --partition=X2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --export=ALL
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source "../../.venv/bin/activate"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1
# Stable output directory intentionally allows safe restart on a new SLURM job.
OUT_DIR="${OUT_DIR:-output}"
python -u run.py --out "${OUT_DIR}" "$@"
echo "GAA mode-space IV simulation terminated normally"
