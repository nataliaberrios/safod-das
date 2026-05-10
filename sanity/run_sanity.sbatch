#!/bin/bash
#SBATCH --job-name=safod_sanity
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Single SLURM script that runs one of {preflight, cc, plot} based on STAGE env var.
# Conventions match run_daily_array.sbatch in the parent repo.
#
# Examples:
#   sbatch --export=ALL,STAGE=preflight run_sanity.sbatch
#   sbatch --export=ALL,STAGE=cc        run_sanity.sbatch
#   sbatch --export=ALL,STAGE=plot      run_sanity.sbatch
#   sbatch --export=ALL,STAGE=all       run_sanity.sbatch     # default
#
# Override target date, channel range, and virtual source by exporting before sbatch:
#   export SANITY_DATE=2024-10-23
#   export SANITY_DATES=2024-10-22:2024-10-28
#   export SANITY_CH_START=150
#   export SANITY_CH_END=800
#   export SANITY_SOURCE_CH=150
#   sbatch --export=ALL run_sanity.sbatch

set -euo pipefail

STAGE="${STAGE:-all}"

# Match parent repo conventions
module --force purge
module load devel math
module load cuda/11.7.1
module load py-pytorch/2.0.0_py39
module unload py-scipy || true

export PYTHONNOUSERSITE=1

REPO_DIR="${REPO_DIR:-/home/groups/ettore88/nberrios/safod-das}"
SANITY_DIR="${REPO_DIR}/sanity"
PYTHON="${PYTHON:-/home/users/nberrios/miniconda3/envs/das/bin/python}"

export PYTHONPATH=/home/users/nberrios/miniconda3/envs/das/lib/python3.9/site-packages:/share/software/user/open/py-pytorch/2.0.0_py39/lib/python3.9/site-packages:${REPO_DIR}/DAS-utilities/python:${REPO_DIR}

cd "${REPO_DIR}"
mkdir -p logs
cd "${SANITY_DIR}"

run_preflight() {
    echo "=== Running preflight.py ==="
    "${PYTHON}" preflight.py
}

run_cc() {
    echo "=== Running sanity_cc.py ==="
    "${PYTHON}" sanity_cc.py
}

run_plot() {
    echo "=== Running sanity_plot.py ==="
    "${PYTHON}" sanity_plot.py
}

case "${STAGE}" in
    preflight) run_preflight ;;
    cc)        run_cc ;;
    plot)      run_plot ;;
    all)
        run_preflight
        run_cc
        run_plot
        ;;
    *)
        echo "Unknown STAGE='${STAGE}'. Use one of: preflight | cc | plot | all"
        exit 1
        ;;
esac

echo "Done."
