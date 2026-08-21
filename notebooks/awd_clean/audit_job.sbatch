#!/bin/bash
# Numerical audit of 2026-08-20: run every audit_test_*.py, probe the real
# $OAK record geometry that the memory findings depend on, and regenerate the
# fk_downgoing products with the corrected (Hermitian) quadrant mask.
#SBATCH -J audit_num
#SBATCH -p serc
#SBATCH -t 02:00:00
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
PY=/home/users/nberrios/miniconda3/envs/das/bin/python

echo "############ 1. synthetic audit suites ############"
for t in audit_test_fk.py audit_test_correlator.py audit_test_geometry.py \
         audit_test_readrecords.py test_block_spectra_chunking.py; do
  echo; echo "===== $t ====="
  $PY -u "$t"; echo "  exit=$?"
done

echo; echo "############ 2. real record geometry (item F inputs) ############"
$PY -u - <<'EOF'
import numpy as np, h5py, sys
sys.path.insert(0, ".")
import ambient_lellouch2019_exact_stack as ex
for arm in ("deepA", "deepC"):
    try:
        rows = ex.deep_rows(arm)
    except Exception as e:
        print("%s: %s" % (arm, e)); continue
    p = ex.corrected_path(rows.file.iloc[0])
    with h5py.File(p, "r") as h:
        g = h["Acquisition/Raw[0]"]
        d = g["RawData"]
        fs = float(g.attrs.get("OutputDataRate", float("nan")))
        dx = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", float("nan")))
        first = d[:64, :4]
        print("%s: %d records, RawData %s %s, fs %.1f Hz, dx %.4f m, %.1f s/record"
              % (arm, len(rows), d.shape, d.dtype, fs, dx, d.shape[0] / fs))
        print("       value range in a 64-sample probe: %g .. %g"
              % (float(first.min()), float(first.max())))
        # what one 4 h block of 343 channels costs at float64
        n = d.shape[0] * 240
        print("       4 h x 343 ch = %.2f GB float64, %.2f GB float32"
              % (n * 343 * 8 / 2**30, n * 343 * 4 / 2**30))
EOF

echo; echo "############ 3. regenerate fk_downgoing with the fixed mask ############"
for spec in "deep_record_section.npz 1675 deep_dense" \
            "deep_record_section_cm.npz 1675 deep_cm" \
            "deep_record_section_pubcm.npz 1675 deep_pub" \
            "deep_record_section_nano.npz 2950 nano"; do
  set -- $spec
  [ -f "$1" ] || { echo "missing $1, skipped"; continue; }
  echo; echo "===== $3 : $1 at $2 m/s ====="
  $PY -u fk_downgoing.py --input "$1" --v-ref "$2" --tag "$3"
done
echo; echo "############ done ############"
