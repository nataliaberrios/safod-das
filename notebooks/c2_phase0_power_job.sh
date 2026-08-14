#!/bin/bash
#SBATCH -J c2_pow
#SBATCH -p serc
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH -o logs/c2_pow_%j.out
#SBATCH -e logs/c2_pow_%j.err

# C2 Phase 0 -- is the tube-wave amplitude candidate list significant at all?
# Reads figures/awd_2026/epoch_stacks_paired_deep_all.npz (2.7 GB), so the
# memory request covers loading the sliced epoch array at float64 plus the
# strain-rate gradient and the bandpass working copies.

cd /home/groups/ettore88/nberrios/safod_das_git/notebooks
ml gcc/12.4.0
/home/users/nberrios/miniconda3/envs/das/bin/python -u c2_phase0_power.py
