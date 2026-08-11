#!/bin/bash
#SBATCH --job-name=fullscan
#SBATCH --partition=serc
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48GB
#SBATCH --array=0-199%40
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
# 574,301 files over 200 tasks = ~2,870 files each, ~1.2 h at 1.5 s/file.
# %40 caps concurrency so the Lustre read rate stays reasonable (~40 x 70 MB/s).
# Only per-day summaries and the top 500 peaks are retained, so output is small.
ml gcc/12.4.0
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/faultzone/repeaters
# NO PIPE. Piping python -u through grep re-buffers the output, so the smoke
# test's log sat at 0 bytes for 27 minutes with the job running fine and no way
# to tell. Write straight to the log; filter when reading it.
/home/users/nberrios/miniconda3/envs/das/bin/python -u full_scan.py
