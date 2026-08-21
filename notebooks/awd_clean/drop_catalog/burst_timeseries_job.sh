#!/bin/bash
#SBATCH -J burst_ts
#SBATCH -p serc
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=20GB
#SBATCH -o ../../logs/burst_ts_%j.out
#SBATCH -e ../../logs/burst_ts_%j.err

# Extract one burst as a continuous Nano time series for the drop-marker figure.
# BURSTS selects which bursts; defaults to 0,10,20,30,40,48.
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/drop_catalog
ml gcc/12.4.0
/home/users/nberrios/miniconda3/envs/das/bin/python -u build_burst_timeseries.py
