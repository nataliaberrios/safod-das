#!/bin/bash
#SBATCH -J census_band
#SBATCH -p serc
#SBATCH -t 01:30:00
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -euo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python
# The resolution prediction: the 3200 m/s target is inside the k=0 main lobe
# below ~12 Hz and outside it above.  If the mechanism is aperture resolution,
# the upper half should separate the fan from k=0 and the lower half should not.
for band in "5 12" "12 20"; do
  set -- $band
  echo "=============== ${1}-${2} Hz, k0 removed ==============="
  CENSUS_FMIN=$1 CENSUS_FMAX=$2 K0_REMOVE=1 $PY -u ambient_apparent_velocity_census.py
done
