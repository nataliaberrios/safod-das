#!/usr/bin/env python3
"""Does the zero-lag peak decay with stacking, or does the arrival just grow?

THE QUESTION. The Deep record section shows a strong zero-lag stripe that buries
the 1675 m/s arrival beyond ~450 m. The natural hypothesis is that stacking longer
will average it away. That is testable and is tested here rather than argued.

WHAT STACKING CAN AND CANNOT DO. Averaging suppresses contributions that are
INCOHERENT between windows, as 1/sqrt(N). A zero-lag correlation peak arises from
noise shared between the source channel and its neighbours, which correlates at
zero lag in EVERY window, so it accumulates at full strength exactly as a real
arrival does and the RATIO between them does not improve. Whether the Deep
zero-lag peak behaves that way is an empirical question about this data.

THE MEASUREMENT. Cross-spectra are accumulated batch by batch, and after each
batch the partial stack is turned into a gather and three quantities recorded:

  A_zero     peak envelope amplitude within +-20 ms of zero lag
  A_arrival  peak envelope amplitude in a window around t = x/1675 m/s,
             excluding the zero-lag zone so the two cannot overlap
  ratio      A_arrival / A_zero      <-- the quantity that decides it

If the ratio GROWS as a power of hours, stacking is winning and more of it helps.
If the ratio is FLAT, the zero-lag peak is as coherent as the arrival and no amount
of stacking will separate them; the fix has to be spatial rather than temporal.

Both amplitudes are also reported in absolute terms, because a flat ratio with
both terms growing means something different from a flat ratio with both static.

Measured only over offsets 100-450 m, where the earlier section showed the
arrival is coherent; including the range where it has already died would dilute
A_arrival with noise and bias the test toward "flat".

Output: deep_zerolag_vs_stack.{npz,png,txt}
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, hilbert, sosfiltfilt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ambient_lellouch2019_exact_stack as ex
import arrival_velocities as av
import deep_cc_steps as steps

STEM = HERE / "deep_zerolag_vs_stack"
SRC = 400
V_ARR = av.V_DEEP_ARRIVAL        # retracted from 1675
OFF_LO, OFF_HI = 100.0, 450.0      # where the arrival is coherent
ZERO_HALF_S = 0.020                # +-20 ms counts as "zero lag"
ARR_HALF_S = 0.030                 # gate around the predicted arrival time
WINDOW_S, STEP_S = 30.0, 15.0
BAND = (5.0, 20.0)
# Derived: the arrival gate sits at x/V_ARR and offsets run to ~450 m, which at
# 1350 m/s is 0.33 s -- inside 0.35 s only barely, and outside it for any slower
# trial. Deriving it removes the coincidence.
MAX_LAG_S = round(av.required_lag_s(500.0, av.V_DEEP_ARRIVAL_RANGE[0]), 2)


def gather_from(acc, n_used, n_fft, fs):
    avg = acc / max(1, n_used)
    full = np.fft.fftshift(np.fft.irfft(avg, n=n_fft, axis=-1), axes=-1)
    lags = (np.arange(n_fft) - n_fft // 2) / fs
    sos = butter(4, list(BAND), btype="bandpass", fs=fs, output="sos")
    keep = np.abs(lags) <= MAX_LAG_S
    return sosfiltfilt(sos, full, axis=-1)[:, keep], lags[keep]


def measure(gather, lags, offsets):
    """Zero-lag and arrival amplitudes, from the envelope, per trace then median."""
    env = np.abs(hilbert(gather, axis=1))
    zero = np.abs(lags) <= ZERO_HALF_S
    a_zero, a_arr = [], []
    for row, x in zip(env, offsets):
        t = x / V_ARR
        gate = (np.abs(lags - t) <= ARR_HALF_S) & (~zero)   # never overlap zero lag
        if gate.sum() < 3:
            continue
        a_zero.append(row[zero].max())
        a_arr.append(row[gate].max())
    return float(np.median(a_zero)), float(np.median(a_arr))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nfiles", type=int, default=288, help="60 s records (288 = 4.8 h)")
    ap.add_argument("--batch", type=int, default=12)
    a = ap.parse_args()

    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    rows = ex.deep_rows("deepA").iloc[: a.nfiles]
    paths = [ex.corrected_path(r) for r in rows.file]
    say("Does stacking remove the zero-lag peak?  Deep fibre, arm deepA")
    say("  %.2f h available in this run (%d x 60 s records)"
        % (len(paths) * 60 / 3600.0, len(paths)))
    say("  %s to %s UTC" % (rows.time.iloc[0], rows.time.iloc[-1]))
    say("  arrival gate at t = x/%.0f m/s, offsets %.0f-%.0f m, zero lag = +-%.0f ms"
        % (V_ARR, OFF_LO, OFF_HI, ZERO_HALF_S * 1e3))
    say("")

    acc = None; n_used = 0; n_fft = fs = dx = None; offsets = None
    curve = []
    for i in range(0, len(paths), a.batch):
        chunk = paths[i:i + a.batch]
        if fs is None:
            with __import__("h5py").File(chunk[0], "r") as h:
                dx = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", 2.0419))
            lo_c, hi_c = int(OFF_LO / dx), int(OFF_HI / dx)
            channels = [SRC] + [SRC + j for j in range(lo_c, hi_c + 1)]
            offsets = np.arange(lo_c, hi_c + 1) * dx
        raw, fs, _ = steps.read_records(chunk, channels)
        x = np.asarray(raw, dtype=np.float32); del raw
        x[:, 1:] = np.diff(x, axis=1) * np.float32(fs); x[:, 0] = 0.0
        n_ram = ex.odd_ram_samples(0.1, fs)
        for j in range(0, x.shape[0], 64):
            blk = x[j:j + 64]
            w = uniform_filter1d(np.abs(blk), size=n_ram, axis=1, mode="nearest")
            sc = float(np.nanmedian(w)) or 1.0
            np.divide(blk, np.maximum(w, np.float32(np.finfo(np.float32).eps * sc)),
                      out=blk)
            del w
        n_win, n_step = int(WINDOW_S * fs), int(STEP_S * fs)
        if acc is None:
            n_fft = 1 << int(np.ceil(np.log2(2 * n_win - 1)))
            acc = np.zeros((len(channels) - 1, n_fft // 2 + 1), dtype=np.complex128)
        for s in range(0, x.shape[1] - n_win + 1, n_step):
            S = np.conj(np.fft.rfft(x[0, s:s + n_win], n=n_fft))
            for r in range(1, x.shape[0]):
                acc[r - 1] += S * np.fft.rfft(x[r, s:s + n_win], n=n_fft)
            n_used += 1
        del x
        g, lags = gather_from(acc, n_used, n_fft, fs)
        z, arr = measure(g, lags, offsets)
        hours = (i + len(chunk)) * 60 / 3600.0
        curve.append((hours, z, arr, arr / z if z else np.nan, n_used))
        say("  %5.2f h | %5d windows | zero-lag %.4e | arrival %.4e | ratio %.4f"
            % (hours, n_used, z, arr, arr / z if z else np.nan))

    c = np.array([(r[0], r[1], r[2], r[3]) for r in curve], dtype=float)
    lh = np.log(c[:, 0])
    b_zero = float(np.polyfit(lh, np.log(np.maximum(c[:, 1], 1e-300)), 1)[0])
    b_arr = float(np.polyfit(lh, np.log(np.maximum(c[:, 2], 1e-300)), 1)[0])
    b_ratio = float(np.polyfit(lh, np.log(np.maximum(c[:, 3], 1e-30)), 1)[0])

    say("")
    say("=== power-law fits against hours stacked ===")
    say("  zero-lag amplitude   N^%+.3f" % b_zero)
    say("  arrival amplitude    N^%+.3f" % b_arr)
    say("  RATIO arrival/zero   N^%+.3f   <- the one that decides it" % b_ratio)
    say("")
    if b_ratio > 0.15:
        say("  THE RATIO IS GROWING. Stacking IS separating the arrival from the")
        say("  zero-lag peak, so more hours will help and the earlier sections were")
        say("  under-stacked. Extend the stack.")
    elif b_ratio < 0.05:
        say("  THE RATIO IS FLAT. The zero-lag peak accumulates exactly as fast as")
        say("  the arrival, which is what happens when both are coherent window to")
        say("  window. No amount of additional stacking will separate them; the")
        say("  remedy has to be SPATIAL (it is a low-wavenumber component) rather")
        say("  than temporal. Consistent with the 2024-25 Nano convergence result,")
        say("  where the raw score grew as N^+0.158 while detectability stayed flat")
        say("  at N^+0.042.")
    else:
        say("  INTERMEDIATE (N^%+.3f). Weakly favourable to more stacking; not" % b_ratio)
        say("  decisive. Report the exponent rather than a verdict.")

    INK, C1, C2, C3 = "#444444", "#0072B2", "#D55E00", "#009E73"
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.3), constrained_layout=True)
    fig.suptitle("Deep fibre: does stacking remove the zero-lag peak?   |   "
                 "%s to %s UTC" % (rows.time.iloc[0].strftime("%Y-%m-%d %H:%M"),
                                   rows.time.iloc[-1].strftime("%Y-%m-%d %H:%M")),
                 fontsize=10.5, y=1.02)
    ax[0].loglog(c[:, 0], c[:, 1], "o-", color=C2, ms=5,
                 label="zero-lag peak (N^%+.2f)" % b_zero)
    ax[0].loglog(c[:, 0], c[:, 2], "s-", color=C1, ms=5,
                 label="arrival at %.0f m/s (N^%+.2f)" % (V_ARR, b_arr))
    ax[0].set(xlabel="hours stacked", ylabel="envelope amplitude",
              title="(a) Both terms")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which="both")
    ax[1].semilogx(c[:, 0], c[:, 3], "o-", color=C3, ms=6)
    ax[1].set(xlabel="hours stacked", ylabel="arrival / zero-lag",
              title="(b) The ratio: N^%+.3f" % b_ratio)
    ax[1].grid(alpha=.3, which="both")
    fig.savefig(str(STEM) + ".png", dpi=300, bbox_inches="tight")
    np.savez_compressed(str(STEM) + ".npz", curve=c, b_zero=b_zero, b_arr=b_arr,
                        b_ratio=b_ratio, offsets=offsets, v_arrival=V_ARR)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
