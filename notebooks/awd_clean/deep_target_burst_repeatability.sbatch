#!/bin/bash
#SBATCH --job-name=awd_deep_burst
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=awd_clean/deep_burst_%j.out
#SBATCH --error=awd_clean/deep_burst_%j.err

set -euo pipefail
ml gcc/12.4.0
source /home/users/nberrios/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
python awd_clean/deep_target_burst_repeatability.py
