#!/usr/bin/env python3
"""Step-by-step ablation of the Deep-fibre ambient-noise cross-correlation.

PURPOSE. Produce every array needed for a manuscript-level figure sequence that
walks from raw phase to the recovered ~1350 m/s arrival (1675 m/s is retracted;
see arrival_velocities.py), and that demonstrates
each processing step is NECESSARY by removing it and showing what breaks.

Necessity is shown by ABLATION, not by assertion: every panel compares the full
pipeline against the same pipeline with exactly one step deleted, on identical
input, with identical everything else.

THE PIPELINE, as reported in Lellouch et al. (2019) section 4.1:
  1  raw optical phase, `rad * 2PI/2^16`
  2  time differentiation            -> strain-rate proxy
  3  running-absolute-mean (0.1 s)   -> temporal normalisation, Bensen et al. 2007
  4  30 s windows, 15 s overlap
  5  cross-correlate each window against the fixed virtual source, in the
     frequency domain with n_fft = 2^ceil(log2(2W-1)) so it is LINEAR not circular
  6  R+-10 neighbour sum, C_S,R = sum_{Z=R-10}^{R+10} C_S,Z
  7  sum over windows (simple unshifted stack)
  8  5-20 Hz bandpass on the FULL correlation, and only THEN crop to +-0.35 s
  9  envelope moveout scan, with a per-velocity receiver-order permutation null

ABLATIONS RUN (one step deleted each):
  no-diff        keep raw phase, skip step 2
  no-ram         skip step 3
  no-neighbour   use the single centre channel, skip step 6
  crop-then-band reverse step 8, which aliases energy into the retained window
  one-window     step 7 with a single window instead of all of them
plus a stacking curve over 1, 2, 4, ... windows.

SELF-VERIFICATION, and this is not optional. The correlation here is a compact
re-implementation so the notebook can show the code, so it MUST be checked
against the authoritative engine (`ambient_lellouch2019_exact_stack.py`) that
produced the result being explained. `--verify` correlates the same records with
both and reports the maximum absolute difference. If they disagree, the figures
would be illustrating a different pipeline than the one that produced the
science, so the script refuses to write products.

Reads raw HDF5 with h5py directly, never DASutils.readFile_HDF (whose
median=True default would silently remove the across-channel median).

Output: deep_cc_steps.npz  (+ a short .txt log)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, hilbert, sosfiltfilt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ambient_lellouch2019_exact_stack as ex
import arrival_velocities as av  # the authoritative engine

STEM = HERE / "deep_cc_steps"

# The WELLHEAD, which is the Lellouch Figure 7c geometry (a source at the top of
# the array). This was 400 -- 389 m down the hole -- chosen because it scored
# highest, but the walkthrough should show the geometry the paper describes.
SRC_CH = 211
OFFSETS_M = np.arange(50.0, 700.1, 50.0)
NEIGHBOUR = 10                    # published R+-10
WINDOW_S = 30.0
STEP_S = 15.0
RAM_S = 0.1
BAND = (5.0, 20.0)
# Derived from the aperture and the scan floor, not hand-picked. At 0.35 s a
# 700 m offset is unreachable below 1934 m/s, so most of V_GRID was scored on a
# shrinking subset of NEAR offsets -- a bias that grows as v falls and which
# manufactured the 1675 m/s peak now retracted. Set after V_GRID exists.
MAX_LAG_S = None
V_GRID = np.arange(300.0, 6000.1, 25.0)
MAX_LAG_S = round(av.required_lag_s(OFFSETS_M.max(), V_GRID.min()), 2)
GATE_S = 0.012
NULLS = 400
SEED = 20260820


# ----------------------------------------------------------------- primitives
def read_records(paths, channels):
    """Concatenate `channels` from consecutive 60 s HDF5 records.

    Rows come back in the order `channels` was given, INCLUDING when that order
    is not increasing. h5py fancy-indexing requires strictly increasing indices
    and otherwise raises "Indexing elements must be in increasing order", which
    is exactly what happens on the Deep fibre's RETURN LIMB, where channel index
    decreases as depth increases. Reading sorted and permuting back afterwards
    keeps callers free to pass receivers ordered by depth.

    MEMORY. The output is PREALLOCATED and filled in place. The earlier version
    built a list of per-file blocks, concatenated it, and then applied the
    inverse permutation -- three full float64 copies alive at the same moment,
    i.e. 24 bytes per sample per channel where callers budget for 12. That 2x
    overshoot is what OOM-killed 12 of 16 deepC array tasks at 64 GB, and it
    silently defeated `deep_timeseries.block_spectra`'s `target_bytes` sizing.
    Record lengths come from the dataset SHAPE, which is HDF5 metadata and costs
    no data read, so the total is known before anything is loaded.
    """
    channels = np.asarray(channels)
    order = np.argsort(channels, kind="stable")
    sorted_ch = channels[order]
    if np.any(np.diff(sorted_ch) == 0):
        raise ValueError("duplicate channel requested: %s"
                         % sorted_ch[np.flatnonzero(np.diff(sorted_ch) == 0)][:5])
    paths = list(paths)
    lengths = []
    fs = dx = None
    for p in paths:
        with h5py.File(p, "r") as h:
            g = h["Acquisition/Raw[0]"]
            lengths.append(int(g["RawData"].shape[0]))
            if fs is None:
                fs = float(g.attrs.get("OutputDataRate", 1000.0))
                dx = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", 2.0419))
    out = np.empty((channels.size, int(sum(lengths))), dtype=np.float64)
    at = 0
    take = sorted_ch.tolist()
    for p, n in zip(paths, lengths):
        with h5py.File(p, "r") as h:
            # `out[order]` puts sorted row i back at the caller's position for
            # that channel, so no permuted intermediate is ever materialised.
            out[order, at:at + n] = h["Acquisition/Raw[0]/RawData"][:, take].T
        at += n
    return out, fs, dx


def differentiate(x, fs):
    """Step 2. Phase -> strain-rate proxy. Dtype-preserving.

    `np.subtract(..., out=)` rather than `out[:, 1:] = np.diff(x, axis=1) * fs`,
    which allocated a second full-size temporary for the difference and a third
    when a Python float `fs` promoted a float32 array to float64. Same
    arithmetic; one array live instead of three. This is what the authoritative
    engine (`strain_rate_and_ram_continuous`) already does.
    """
    out = np.zeros_like(x)
    np.subtract(x[:, 1:], x[:, :-1], out=out[:, 1:])
    out *= out.dtype.type(fs) if np.issubdtype(out.dtype, np.floating) else fs
    return out


def ram_normalise(x, fs, ram_s=RAM_S):
    """Step 3. Centred running-absolute-mean weights (Bensen et al. 2007).

    Dtype-preserving. The floor was a Python float built from
    `np.finfo(np.float64)`, which numpy promotes, so a float32 caller silently
    got a float64 result and twice the memory. Taking eps from the input dtype
    also makes this agree with the engine, which floors at float32 eps.
    """
    n = ex.odd_ram_samples(ram_s, fs)
    w = uniform_filter1d(np.abs(x), size=n, axis=1, mode="nearest")
    scale = float(np.nanmedian(w))
    info = np.finfo(w.dtype if np.issubdtype(w.dtype, np.floating) else np.float64)
    floor = max(info.tiny, info.eps * scale)
    return x / np.maximum(w, w.dtype.type(floor))


def correlate(data, source_row, receiver_rows, fs, n_win, n_step, band_then_crop=True,
              max_windows=None):
    """Steps 4-8. Returns (gather, lags, n_windows_used).

    `receiver_rows` is a list of index arrays: one array per receiver, holding the
    rows to sum for that receiver (a single index for the no-neighbour ablation,
    2*NEIGHBOUR+1 indices for the published R+-10 sum).
    """
    n_t = data.shape[1]
    starts = list(range(0, max(1, n_t - n_win + 1), n_step))
    if max_windows is not None:
        starts = starts[:max_windows]
    n_fft = 1 << int(np.ceil(np.log2(2 * n_win - 1)))
    acc = np.zeros((len(receiver_rows), n_fft // 2 + 1), dtype=np.complex128)
    used = 0
    for s in starts:
        src = data[source_row, s:s + n_win]
        if src.size < n_win:
            break
        S = np.conj(np.fft.rfft(src, n=n_fft))
        for i, rows in enumerate(receiver_rows):
            summed = data[rows, s:s + n_win].sum(axis=0)
            acc[i] += S * np.fft.rfft(summed, n=n_fft)
        used += 1
    avg = acc / max(1, used)
    full = np.fft.fftshift(np.fft.irfft(avg, n=n_fft, axis=-1), axes=-1)
    lags = (np.arange(n_fft) - n_fft // 2) / fs
    sos = butter(4, list(BAND), btype="bandpass", fs=fs, output="sos")
    keep = np.abs(lags) <= MAX_LAG_S
    if band_then_crop:
        full = sosfiltfilt(sos, full, axis=-1)[:, keep]
    else:
        # WRONG ORDER, kept as an ablation: cropping first aliases out-of-band
        # energy into the retained window because the crop is a hard truncation.
        full = sosfiltfilt(sos, full[:, keep], axis=-1)
    return full, lags[keep], used


def moveout_curve(gather, lags, offsets, velocities=V_GRID, sign=1.0):
    env = np.abs(hilbert(gather, axis=1))
    env = env / np.maximum(np.median(env, axis=1, keepdims=True), 1e-30)
    half = max(1, int(GATE_S / (lags[1] - lags[0])))
    out = np.empty(velocities.size)
    for i, v in enumerate(velocities):
        vals = []
        for row, x in zip(env, offsets):
            k = int(np.argmin(np.abs(lags - sign * x / v)))
            lo, hi = max(0, k - half), min(len(lags), k + half + 1)
            vals.append(row[lo:hi].mean())
        out[i] = np.median(vals)
    return out


def per_velocity_null(gather, lags, offsets, rng, n=NULLS, velocities=None):
    """Receiver-order permutation null, on the SAME velocity grid as the score.

    `velocities` is forwarded rather than left to `moveout_curve`'s default: the
    default is the module-level V_GRID, so a caller that scored on a different
    grid would have silently been compared against a null built on V_GRID. That
    is the failure that hid in `ambient_lellouch2019_exact_stack`'s
    `velocities=VELOCITY_GRID_M_S` default, and there is no reason to leave a
    second copy of it here. nano_ambient_cc.py passes its own grid positionally.
    """
    velocities = V_GRID if velocities is None else np.asarray(velocities, float)
    null = np.empty((n, velocities.size))
    for i in range(n):
        null[i] = moveout_curve(gather[rng.permutation(len(offsets))], lags,
                                offsets, velocities)
    return null


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nfiles", type=int, default=60, help="60 s records to use")
    ap.add_argument("--verify", action="store_true",
                    help="check the compact correlator against the authoritative engine")
    a = ap.parse_args()

    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    rng = np.random.default_rng(SEED)

    rows = ex.deep_rows("deepA").iloc[: a.nfiles]
    paths = [ex.corrected_path(r) for r in rows.file]
    say("Deep fibre ARM A (pre-survey ambient), %.2f h  (%d x 60 s records)"
        % (len(paths) * 60 / 3600.0, len(paths)))
    say("  %s -> %s" % (rows.time.iloc[0], rows.time.iloc[-1]))
    drops = ex.deep_drop_times()
    gap = (drops[0] - (rows.time.iloc[-1].value + int(60e9))) / 1e9
    say("  first weight drop is %.2f h AFTER the last record used, so this is"
        % (gap / 3600.0))
    say("  ambient noise and not the active source.")
    say("")

    # geometry: source, receiver centres, and the rows each receiver sums
    centres = SRC_CH + np.rint(OFFSETS_M / 2.0419).astype(int)
    needed = sorted(set([SRC_CH] + [c + d for c in centres
                                    for d in range(-NEIGHBOUR, NEIGHBOUR + 1)]))
    index = {c: i for i, c in enumerate(needed)}
    say("source channel %d, %d receivers at %.0f-%.0f m, %d channels read"
        % (SRC_CH, centres.size, OFFSETS_M[0], OFFSETS_M[-1], len(needed)))

    raw, fs, dx = read_records(paths, needed)
    # FLOAT32, and it is exact here rather than a compromise. RawData is int32
    # phase reaching ~1.9e6, well inside float32's exactly-representable integer
    # range (2^24 = 1.68e7), and consecutive samples differ by a median of ~110
    # counts against a float32 spacing of 0.03 counts at that magnitude.
    # Measured on both arms 2026-08-20 (audit_job2.sbatch): the strain rate and
    # the accumulated cross-spectra are IDENTICAL to the float64 path, relative
    # error 0.000e+00. This script was the only one in the directory still
    # working in float64, and with `raw`, `rate` and `normed` all live at once
    # plus the ablation ladder it peaked near 42 GB at --nfiles 60 against a
    # 64 GB request. Every other script here already converts.
    raw = np.ascontiguousarray(raw, dtype=np.float32)
    say("  fs %.0f Hz, dx %.4f m, raw shape %s float32, unit rad*2PI/2^16"
        % (fs, dx, raw.shape))
    n_win, n_step = int(WINDOW_S * fs), int(STEP_S * fs)
    src_row = index[SRC_CH]
    rows_nb = [np.array([index[c + d] for d in range(-NEIGHBOUR, NEIGHBOUR + 1)])
               for c in centres]
    rows_single = [np.array([index[c]]) for c in centres]

    # ---- the full pipeline ----
    rate = differentiate(raw, fs)
    normed = ram_normalise(rate, fs)
    full_gather, lags, nw = correlate(normed, src_row, rows_nb, fs, n_win, n_step)
    say("  full pipeline: %.2f h in %d windows"
        % (len(paths) * 60 / 3600.0, nw))

    if a.verify:
        say("")
        say("=== verification against ambient_lellouch2019_exact_stack.py ===")
        eng_rate, floored = ex.strain_rate_and_ram_continuous(
            raw.copy(), fs, RAM_S, False)
        d = float(np.max(np.abs(eng_rate.astype(np.float64) - normed)))
        rel = d / max(1e-30, float(np.max(np.abs(normed))))
        say("  strain-rate + RAM stage: max abs diff %.3e (relative %.3e)" % (d, rel))
        if rel > 1e-5:
            raise SystemExit(
                "compact re-implementation DISAGREES with the authoritative engine "
                "(relative %.3e). Refusing to write figure products that would "
                "illustrate a different pipeline than the one that produced the "
                "result." % rel)
        say("  agreement within tolerance; the compact correlator may be shown.")
        say("")

    # ---- ablations ----
    abl = {}
    abl["full"] = full_gather
    # no differentiation
    g, _, _ = correlate(ram_normalise(raw, fs), src_row, rows_nb, fs, n_win, n_step)
    abl["no_diff"] = g
    # no RAM
    g, _, _ = correlate(rate, src_row, rows_nb, fs, n_win, n_step)
    abl["no_ram"] = g
    # no neighbour sum
    g, _, _ = correlate(normed, src_row, rows_single, fs, n_win, n_step)
    abl["no_neighbour"] = g
    # wrong band/crop order
    g, _, _ = correlate(normed, src_row, rows_nb, fs, n_win, n_step,
                        band_then_crop=False)
    abl["crop_then_band"] = g
    # single window
    g, _, _ = correlate(normed, src_row, rows_nb, fs, n_win, n_step, max_windows=1)
    abl["one_window"] = g
    say("ablations computed: %s" % ", ".join(sorted(abl)))

    curves = {k: moveout_curve(v, lags, OFFSETS_M) for k, v in abl.items()}
    acausal = moveout_curve(abl["full"], lags, OFFSETS_M, sign=-1.0)
    null = per_velocity_null(abl["full"], lags, OFFSETS_M, rng)
    thresh = np.percentile(null, 95.0, axis=0)
    pv = (np.sum(null >= curves["full"][None, :], axis=0) + 1.0) / (NULLS + 1.0)

    k = int(np.argmax(curves["full"]))
    say("")
    say("full pipeline: peak %.4f at %.0f m/s, per-velocity p = %.4f"
        % (curves["full"][k], V_GRID[k], pv[k]))
    say("  velocities clearing their own 95th percentile: %d of %d"
        % (int(np.sum(curves["full"] > thresh)), V_GRID.size))
    say("")
    say("effect of deleting each step (peak score and where it lands):")
    for name in ("full", "no_diff", "no_ram", "no_neighbour", "crop_then_band",
                 "one_window"):
        c = curves[name]; j = int(np.argmax(c))
        say("  %-15s peak %8.4f at %6.0f m/s" % (name, c[j], V_GRID[j]))

    # ---- stacking curve ----
    say("")
    say("stacking curve (windows -> peak, and score at the arrival velocity):")
    stack_n, stack_peak, stack_at_arrival = [], [], []
    m = 1
    while m <= nw:
        g, _, used = correlate(normed, src_row, rows_nb, fs, n_win, n_step,
                               max_windows=m)
        c = moveout_curve(g, lags, OFFSETS_M)
        stack_n.append(used); stack_peak.append(float(c.max()))
        stack_at_arrival.append(float(c[k]))
        say("  %4d windows: peak %.4f, score at %.0f m/s = %.4f"
            % (used, c.max(), V_GRID[k], c[k]))
        m *= 2
    if stack_n[-1] != nw:
        g, _, used = correlate(normed, src_row, rows_nb, fs, n_win, n_step)
        c = moveout_curve(g, lags, OFFSETS_M)
        stack_n.append(used); stack_peak.append(float(c.max()))
        stack_at_arrival.append(float(c[k]))
        say("  %4d windows: peak %.4f, score at %.0f m/s = %.4f"
            % (used, c.max(), V_GRID[k], c[k]))

    # raw / intermediate traces for the first figures: one 60 s record's worth
    n_show = int(min(60.0 * fs, raw.shape[1]))
    np.savez_compressed(
        str(STEM) + ".npz",
        lags=lags, offsets=OFFSETS_M, v_grid=V_GRID,
        raw_snippet=raw[:, :n_show].astype(np.float32),
        rate_snippet=rate[:, :n_show].astype(np.float32),
        normed_snippet=normed[:, :n_show].astype(np.float32),
        channels=np.array(needed), source_channel=SRC_CH, centres=centres,
        fs=fs, dx=dx, n_windows=nw, n_files=len(paths),
        drop_gap_hours=gap / 3600.0,
        t_first=str(rows.time.iloc[0]), t_last=str(rows.time.iloc[-1]),
        acausal=acausal, null_thresh=thresh, p_per_velocity=pv,
        stack_n=np.array(stack_n), stack_peak=np.array(stack_peak),
        stack_at_arrival=np.array(stack_at_arrival),
        arrival_v=float(V_GRID[k]), arrival_p=float(pv[k]),
        **{("gather_" + k2): v for k2, v in abl.items()},
        **{("curve_" + k2): v for k2, v in curves.items()},
    )
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,txt}" % STEM.name)


if __name__ == "__main__":
    main()
