#!/bin/bash
#SBATCH -J drop_nb
#SBATCH -p serc
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH -o ../../logs/drop_nb_%j.out
#SBATCH -e ../../logs/drop_nb_%j.err

# Regenerate AWD_drop_catalog.ipynb from its build script, then execute it with
# the das kernel so the committed notebook carries its figures.
#
# The build step needs the SYSTEM python3 (nbformat); execution needs the das
# env (the das kernel, and no nbformat). That split is the reason for two
# interpreters here -- see CLAUDE.md.

cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean/drop_catalog
ml gcc/12.4.0
PY=/home/users/nberrios/miniconda3/envs/das/bin/python

echo "=== regenerate the notebook (system python3: has nbformat) ==="
python3 -u build_drop_notebook.py || exit 1

echo
echo "=== execute with outputs saved (das kernel) ==="
$PY -u ../manuscript/execute_notebook.py AWD_drop_catalog.ipynb \
    --kernel das --timeout 900 --dpi 110
