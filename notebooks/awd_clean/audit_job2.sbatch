#!/bin/bash
# Follow-up probe: does the float32 conversion before np.diff lose the signal?
# The Deep RawData is int32 with |values| ~ 1e6, where float32 resolution is
# ~0.06-0.12 counts. If the sample-to-sample difference is only a few counts,
# converting to float32 BEFORE differencing quantises the strain rate.
# Every script in this directory except deep_cc_steps.py does exactly that.
#SBATCH -J audit_f32
#SBATCH -p serc
#SBATCH -t 00:40:00
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
set -uo pipefail
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
/home/users/nberrios/miniconda3/envs/das/bin/python -u - <<'EOF'
import sys
import numpy as np, h5py
sys.path.insert(0, ".")
import ambient_lellouch2019_exact_stack as ex
from scipy.ndimage import uniform_filter1d

fs = 1000.0
n_ram = ex.odd_ram_samples(0.1, fs)
for arm in ("deepA", "deepC"):
    rows = ex.deep_rows(arm)
    p = ex.corrected_path(rows.file.iloc[0])
    ch = list(range(211, 211 + 64))              # in-hole, near-vertical
    with h5py.File(p, "r") as h:
        raw_i = h["Acquisition/Raw[0]/RawData"][:, ch].T.astype(np.int64)
    print("\n===== %s : %s  %s =====" % (arm, p.name, raw_i.shape))
    print("  raw int32 range        : %d .. %d" % (raw_i.min(), raw_i.max()))
    print("  |raw| median           : %.3e" % np.median(np.abs(raw_i)))
    eps32 = np.spacing(np.float32(np.median(np.abs(raw_i))))
    print("  float32 spacing there  : %.4f counts" % eps32)

    d_exact = np.diff(raw_i, axis=1).astype(np.float64)
    print("  |sample-to-sample diff|: median %.2f, p05 %.2f, p95 %.2f counts"
          % (np.median(np.abs(d_exact)), np.percentile(np.abs(d_exact), 5),
             np.percentile(np.abs(d_exact), 95)))
    print("  -> float32 quantisation is %.2f%% of the median difference"
          % (100.0 * eps32 / max(np.median(np.abs(d_exact)), 1e-30)))

    # the two pipelines, on identical input
    r64 = np.diff(raw_i.astype(np.float64), axis=1) * fs
    r32 = np.diff(raw_i.astype(np.float32), axis=1) * np.float32(fs)
    err = r32.astype(np.float64) - r64
    print("  strain rate: rel RMS error of the float32 path = %.3e"
          % (np.sqrt(np.mean(err ** 2)) / np.sqrt(np.mean(r64 ** 2))))

    # does it survive the running-abs-mean and the correlation?
    def ram(x):
        x = np.asarray(x, dtype=np.float64)
        w = uniform_filter1d(np.abs(x), size=n_ram, axis=1, mode="nearest")
        sc = float(np.nanmedian(w)) or 1.0
        return x / np.maximum(w, np.finfo(np.float32).eps * sc)
    n_win = int(30 * fs); n_fft = 1 << int(np.ceil(np.log2(2 * n_win - 1)))
    out = {}
    for lab, r in (("f64", r64), ("f32", r32)):
        x = ram(r)
        acc = np.zeros((x.shape[0] - 1, n_fft // 2 + 1), dtype=np.complex128)
        nw = 0
        for s in range(0, x.shape[1] - n_win + 1, int(15 * fs)):
            W = np.fft.rfft(x[:, s:s + n_win], n=n_fft, axis=1)
            acc += np.conj(W[0]) * W[1:]
            nw += 1
        out[lab] = acc / nw
    rel = float(np.linalg.norm(out["f32"] - out["f64"]) /
                np.linalg.norm(out["f64"]))
    print("  cross-spectra (%d windows): rel diff float32 vs float64 = %.3e" % (nw, rel))
    print("  VERDICT: %s" % ("float32 conversion is SAFE here" if rel < 1e-3
                             else "float32 conversion CHANGES the result"))
EOF
