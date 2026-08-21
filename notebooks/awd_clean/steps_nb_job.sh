#!/bin/bash
#SBATCH -J stepsnb
#SBATCH -p serc
#SBATCH -t 04:00:00
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
# Regenerate the ambient-CC walkthrough end to end.
#   1. deep_cc_steps.npz was built with the +-0.35 s lag window and a source
#      389 m down the hole. Rebuild at the derived window, from the WELLHEAD.
#   2. rebuild the notebook (system python3: the das env has no nbformat)
#   3. EXECUTE it, so the committed notebook actually carries its figures.
#      run_nb_cells.py is not enough -- it writes nothing back, so it cannot
#      tell you the notebook renders. execute_notebook.py exits non-zero if a
#      cell that draws a figure captured none.
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
$PY -u deep_cc_steps.py --nfiles 300
python3 -u build_ambient_cc_steps_notebook.py
$PY -u manuscript/execute_notebook.py Ambient_CC_Steps.ipynb || \
  python3 -u manuscript/execute_notebook.py Ambient_CC_Steps.ipynb
