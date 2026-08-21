#!/bin/bash
#SBATCH -J burst_ts
#SBATCH -p serc
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH -o ../../logs/burst_ts_%j.out
#SBATCH -e ../../logs/burst_ts_%j.err

# Extract one burst as a continuous Nano time series for the drop-marker figure.
# BURST selects which burst; defaults to 0.
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/drop_catalog
ml gcc/12.4.0
/home/users/nberrios/miniconda3/envs/das/bin/python -u build_burst_timeseries.py
