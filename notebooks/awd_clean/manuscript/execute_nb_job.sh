#!/bin/bash
#SBATCH -J exec_nb
#SBATCH -p serc
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128GB
#SBATCH -o ../../logs/exec_nb_%j.out
#SBATCH -e ../../logs/exec_nb_%j.err

# Execute AWD_reproduce_analysis.ipynb with outputs saved, so the committed
# notebook shows its figures without the reader running it first.  Memory is
# sized for canonical_epoch_stacks_paired_deep_all.npz (2.7 GB) plus working
# copies.
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/manuscript
ml gcc/12.4.0
/home/users/nberrios/miniconda3/envs/das/bin/python -u execute_notebook.py \
    AWD_reproduce_analysis.ipynb --kernel das --timeout 2400 --dpi 110
