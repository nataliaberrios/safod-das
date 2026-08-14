#!/bin/bash
#SBATCH -J nb_rebuild
#SBATCH -p serc
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128GB
#SBATCH -o ../../logs/nb_rebuild_%j.out
#SBATCH -e ../../logs/nb_rebuild_%j.err

cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/manuscript
ml gcc/12.4.0
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
echo "=== regenerate notebook from build_notebook.py ==="
$PY -u build_notebook.py
echo
echo "=== execute with outputs saved ==="
$PY -u execute_notebook.py AWD_reproduce_analysis.ipynb --kernel das --timeout 2400 --dpi 110
