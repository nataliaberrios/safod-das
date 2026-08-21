#!/usr/bin/env python3
"""Does streaming `block_spectra` over sub-chunks equal reading the block whole?

`deep_timeseries.block_spectra` was changed on 2026-08-20 to read a block in
sub-chunks, because reading 4 h whole cost ~40 GB and OOM-killed 12 of 16 array
tasks. Chunking a windowed accumulation is the kind of change that quietly loses
data: any window straddling a chunk join is dropped unless a tail is carried
across, and a 1-in-80 window loss would never show up as an error -- only as a
slightly weaker stack that looks like an honest result.

The first implementation split in TIME. This test rejected it: window counts
matched (79/79/79, so the tail carry worked) but the cross-spectra differed by
1.96e-3, because the 0.1 s running-abs-mean saw a chunk edge instead of real
data for 0.05 s either side of every join -- 20 joins over 1200 s predicts
1.7e-3, which is what was measured. That difference is physically negligible,
and that is exactly why it was unacceptable: the only tolerance that admits it
is one loose enough to also admit a real bug.

`block_spectra` now splits in CHANNEL, which has no joins at all, so the two
properties are both EXACT:

  WINDOW COUNT equal, and equal to the analytic count for the full record.
  This is the property that catches a dropped window.

  CROSS-SPECTRA equal to 1e-12, i.e. bit-identical up to summation order.
  If this ever loosens, the fix is to find the bug, not to raise the tolerance.

Synthetic data only -- no $OAK reads -- so this runs in seconds and can gate a
job before it spends node hours.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import deep_timeseries as dts

FS = 500.0
REC_N = int(60 * FS)          # 60 s records, as the Deep fibre writes them
N_REC = 20                    # 20 min
N_CH = 8


def main() -> int:
    rng = np.random.default_rng(20260820)
    src = rng.standard_normal(N_REC * REC_N).astype(np.float32)
    master = np.empty((N_CH, N_REC * REC_N), dtype=np.float32)
    for c in range(N_CH):
        lag = 7 * c                                   # a clean moveout to correlate
        master[c] = np.roll(src, lag) + 0.3 * rng.standard_normal(src.size)
    master = np.cumsum(master, axis=1)                # so np.diff has something to undo

    def fake_read(paths, channels):
        idx = [int(Path(p).stem) for p in paths]
        cols = np.concatenate([np.arange(i * REC_N, (i + 1) * REC_N) for i in idx])
        return master[np.asarray(channels)][:, cols], FS, 1.0

    dts.steps.read_records = fake_read
    dts.SOURCE_CH = 0
    paths = ["%d" % i for i in range(N_REC)]
    rx = list(range(1, N_CH))

    whole, nw_whole, fs, n_fft = dts.block_spectra(paths, rx, (5.0, 20.0),
                                                   target_bytes=1e12)
    chunk, nw_chunk, _, _ = dts.block_spectra(paths, rx, (5.0, 20.0),
                                              target_bytes=6.0 * REC_N * N_REC * 4 * 3)

    n_win, n_step = int(dts.WINDOW_S * FS), int(dts.STEP_S * FS)
    total = N_REC * REC_N
    expect = (total - n_win) // n_step + 1

    rel = float(np.linalg.norm(chunk - whole) / max(np.linalg.norm(whole), 1e-30))
    ok = True

    print("streaming block_spectra vs reading the block whole")
    print("  %d channels, %d x 60 s records = %.0f s at %.0f Hz"
          % (N_CH, N_REC, total / FS, FS))
    print("  window %.0f s, step %.0f s, n_fft %d" % (dts.WINDOW_S, dts.STEP_S, n_fft))
    print()
    print("  windows, one group      : %d" % nw_whole)
    print("  windows, many groups    : %d" % nw_chunk)
    print("  windows, analytic       : %d" % expect)
    if nw_whole == nw_chunk == expect:
        print("  PASS  every group saw the same complete window sequence")
    else:
        print("  FAIL  window counts disagree across channel groups")
        ok = False

    print("  cross-spectra rel diff  : %.3e" % rel)
    if rel < 1e-12:
        print("  PASS  splitting by channel is exact, as it must be")
    else:
        print("  FAIL  channel grouping changed the accumulated cross-spectra.")
        print("        Do NOT raise this tolerance -- channel groups are")
        print("        independent, so any difference at all is a bug.")
        ok = False

    # A dropped window is a ~1 % effect on nw and would pass a loose amplitude
    # check, so confirm the count test can actually see one.
    print()
    print("  (a single dropped window would move the count to %d, which the"
          % (expect - 1))
    print("   exact-equality test above rejects)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
