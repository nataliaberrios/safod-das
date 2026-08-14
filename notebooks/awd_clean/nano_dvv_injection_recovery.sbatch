#!/bin/bash
#SBATCH --job-name=nano_dvv
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=awd_clean/nano_dvv_%j.out
#SBATCH --error=awd_clean/nano_dvv_%j.err

set -euo pipefail
source /home/users/nberrios/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
python awd_clean/nano_dvv_injection_recovery.py --stage all
