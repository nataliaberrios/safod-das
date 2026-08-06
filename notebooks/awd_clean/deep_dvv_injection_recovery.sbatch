#!/bin/bash
#SBATCH --job-name=awd_deep_dvv
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=awd_clean/deep_dvv_%j.out
#SBATCH --error=awd_clean/deep_dvv_%j.err

set -euo pipefail
ml gcc/12.4.0
source /home/users/nberrios/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
python -u awd_clean/deep_dvv_injection_recovery.py --stage all
