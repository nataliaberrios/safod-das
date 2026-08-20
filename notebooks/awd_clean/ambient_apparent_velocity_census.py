#!/usr/bin/env python3
"""Where does the along-fibre energy sit in apparent velocity, INCLUDING k ~ 0?

THE CLAIM UNDER TEST.  Figure 7d shows our constant-offset correlations peaking at
exactly zero lag at every depth (0.0 ms at 50, 250, 550 and 700 m), while
Lellouch's released 7d peaks at the P travel time and migrates with depth
(20.7 -> 11.5 ms, giving 2416 -> 4357 m/s). A zero-lag peak between channels 50 m
apart means the signals are simultaneous: INFINITE apparent velocity along the
fibre. The reading is therefore that the 2024-25 field is dominated by energy with
no along-fibre moveout, where Lellouch's was dominated by downgoing energy.

This tests that directly, and it tests the one thing the earlier census could not:
`ambient_fk_energy_census.py` excluded the k = 0 column via |k| > 1e-9, which the
audit found discards 22.1 % of the in-band energy. Energy at k ~ 0 IS energy at
infinite apparent velocity -- the quantity this claim is about. The old census
threw away the measurement that bears on the answer.

WHAT IS FIXED RELATIVE TO THE WITHDRAWN CENSUS
  - k = 0 is INCLUDED, in an explicit "no along-fibre moveout" bin.
  - A geometric baseline is computed: the fraction of in-band (f,k) CELLS each
    velocity bin occupies, so shares are reported against chance rather than raw.
    The old census had none, and 98.4 % of cells lie below 1500 m/s by
    construction, which is why its headline could not distinguish SAFOD from
    white noise.
  - Hann tapers in BOTH space and time before the 2-D FFT, so the fan measurement
    is not dominated by Dirichlet leakage from the k ~ 0 peak. Untapered, the
    audit measured a 24 % bias.
  - Energy is summed, not per-window variance-normalised, so the budget is a real
    budget.

IDENTICAL PROCESSING ON BOTH ARMS, asserted at run time. Five of six withdrawn
claims on 2026-08-14 came from two arms of a comparison being processed
differently -- most recently because `extract_all.py` takes DASutils
`median=True` by default and stripped the common mode from one epoch only. Here
both epochs are read RAW (h5py for 2024-25, np.fromfile for the 2017 release, no
DASutils anywhere), pushed through one `analyse()` function, decimated to a common
250 Hz, cut to a common aperture in METRES, and trace-edited by the same rule. The
operations applied to each arm are recorded and compared before any result prints.

Output: ambient_apparent_velocity_census.{npz,png,txt}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, detrend, resample_poly, sosfiltfilt

HERE = Path(__file__).resolve().parent
STEM = HERE / "ambient_apparent_velocity_census"
CSV = Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv")
LEL = Path("/scratch/users/nberrios/lellouch2017")   # raw 2017 binaries (correlograms live in repo lellouch_traces/)

FMIN, FMAX = 5.0, 20.0
FS_COMMON = 250.0            # the 2017 release rate; 2024-25 decimates 500 -> 250
APERTURE_M = 700.0           # common aperture both arms can supply after editing
WIN_S = 2.0                  # the 2017 pre-event window is short; both use 2 s
EDIT_LO, EDIT_HI = 0.2, 5.0
CH_LO_2024 = 23              # wellhead (G0)

# apparent-velocity bins, INCLUDING an explicit no-moveout bin at the top
V_EDGES = np.array([0, 500, 1000, 1500, 2000, 2500, 3000, 4000,
                    6000, 10000, 30000, np.inf])


def corrected_path(p):
    return Path(str(p).replace("/data/SAFODAS1-harddrive-transfer",
                               "/data/SAFOD/SAFODAS1-harddrive-transfer"))


def edit_traces(x):
    """Replace outlier channels by linear interpolation from their neighbours.

    NOT by deletion. A 2-D FFT assumes UNIFORM spatial sampling, so removing
    channels would leave an irregularly sampled array and corrupt the wavenumber
    axis -- the very axis this script measures. Interpolation keeps dx uniform
    while removing the glitches (the 2017 release requires trace editing: M1p33
    carries 16 channels at ~1e6x the median).
    """
    rms = np.sqrt((x ** 2).mean(axis=1))
    med = np.median(rms[rms > 0]) if np.any(rms > 0) else 0.0
    if med <= 0:
        return x, 0
    bad = ~((rms > EDIT_LO * med) & (rms < EDIT_HI * med))
    if bad.all():
        return x, 0
    y = x.copy()
    idx = np.arange(len(rms))
    good = ~bad
    for t in range(x.shape[1]):
        y[bad, t] = np.interp(idx[bad], idx[good], x[good, t])
    return y, int(bad.sum())


def analyse(raw, fs, dx, label, ops, log):
    """One code path. Both epochs go through exactly this, and `ops` records it."""
    steps = []
    x = np.asarray(raw, dtype=float)

    x, dropped = edit_traces(x); steps.append("trace_edit_interp(0.2-5.0 x median RMS)")
    log.append("    %-10s trace edit interpolated %d of %d channels" % (label, dropped, x.shape[0]))

    x = np.diff(x, axis=1, prepend=x[:, :1]); steps.append("diff(strain_rate)")

    if abs(fs - FS_COMMON) > 1e-9:
        from math import gcd
        up, dn = int(FS_COMMON), int(fs); g = gcd(up, dn)
        x = resample_poly(x, up // g, dn // g, axis=1); fs = FS_COMMON
        steps.append("resample_to_%.0fHz" % FS_COMMON)
    else:
        steps.append("resample_to_%.0fHz(noop)" % FS_COMMON)

    x = detrend(x, axis=1, type="linear"); steps.append("detrend")
    x = sosfiltfilt(butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos"),
                    x, axis=1); steps.append("bandpass_%g-%g" % (FMIN, FMAX))

    nch = int(round(APERTURE_M / dx))
    if x.shape[0] < nch:
        raise SystemExit("%s: only %d channels, need %d" % (label, x.shape[0], nch))
    x = x[:nch]; steps.append("aperture_%.0fm(%d ch)" % (APERTURE_M, nch))

    nw = int(WIN_S * fs)
    if x.shape[1] < nw:
        raise SystemExit("%s: record shorter than the %.1f s window" % (label, WIN_S))

    # Hann in space AND time before the 2-D FFT -- controls leakage from k ~ 0
    wt = np.hanning(nw)[None, :]
    wx = np.hanning(nch)[:, None]
    steps.append("hann_space_and_time")

    P = None; nwin = 0
    for s in range(0, x.shape[1] - nw + 1, nw):
        seg = x[:, s:s + nw] * wt * wx
        F = np.fft.fftshift(np.fft.fft2(seg))
        P = np.abs(F) ** 2 if P is None else P + np.abs(F) ** 2
        nwin += 1
    steps.append("sum_power(over %d windows, no per-window normalisation)" % nwin)

    ops[label] = steps
    k = np.fft.fftshift(np.fft.fftfreq(nch, dx))
    f = np.fft.fftshift(np.fft.fftfreq(nw, 1.0 / fs))
    return P, k, f, nwin


def budget(P, k, f):
    """Energy and CELL COUNT per apparent-velocity bin, k = 0 included."""
    K, F = np.meshgrid(k, f, indexing="ij")
    band = (np.abs(F) >= FMIN) & (np.abs(F) <= FMAX)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(np.abs(K) > 0, np.abs(F) / np.abs(K), np.inf)
    tot_e = P[band].sum(); tot_c = band.sum()
    rows = []
    for lo, hi in zip(V_EDGES[:-1], V_EDGES[1:]):
        # The top bin must be CLOSED at infinity. k = 0 gives v = inf, and
        # `inf < inf` is False, so a half-open top bin silently discards the
        # entire k = 0 column -- which is precisely the quantity this script
        # exists to measure, and which the previous census also dropped.
        m = band & (v >= lo) & ((v <= hi) if np.isinf(hi) else (v < hi))
        rows.append((lo, hi, 100 * P[m].sum() / tot_e, 100 * m.sum() / tot_c))
    covered = sum(r[2] for r in rows)
    if not (99.0 <= covered <= 101.0):
        raise SystemExit("binning does not conserve energy: %.2f%% covered" % covered)
    return rows


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    ops = {}

    say("Apparent-velocity budget along the fibre, k = 0 INCLUDED")
    say("  band %g-%g Hz | common rate %.0f Hz | common aperture %.0f m | %.1f s windows"
        % (FMIN, FMAX, FS_COMMON, APERTURE_M, WIN_S))
    say("  Hann taper in space and time; energy summed, not per-window normalised")
    say("")

    # ---- 2024-25, raw HDF5, no DASutils ----
    db = pd.read_csv(CSV, sep=r"\s+"); db = db[db.nSamples > 0].copy()
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    day = db[db.t.dt.strftime("%Y-%m-%d") == "2024-12-20"].sort_values("t").reset_index(drop=True)
    seg = []
    for row in day.iloc[::240].head(6).itertuples(index=False):
        f = corrected_path(row.file)
        if not f.is_file():
            continue
        with h5py.File(f, "r") as h:
            g = h["Acquisition/Raw[0]"]
            d = g["RawData"][:, CH_LO_2024:].astype(np.float32).T
            fs24 = float(g.attrs.get("OutputDataRate", 500.0))
            dx24 = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", 1.0))
        seg.append(d)
    x24 = np.concatenate(seg, axis=1); del seg
    P24, k24, f24, n24 = analyse(x24, fs24, dx24, "2024-25", ops, log)

    # ---- 2017, raw binary release ----
    parts = []
    for name in ("M1p33", "M2p46"):
        p = LEL / name
        if p.is_file():
            v = np.fromfile(p, dtype="<f4").reshape(1250, 800).T
            parts.append(v[:, :int(2.5 * 250.0)])          # pre-event only
    if not parts:
        raise SystemExit("2017 release not found")
    x17 = np.concatenate(parts, axis=1)
    P17, k17, f17, n17 = analyse(x17, 250.0, 1.0, "2017", ops, log)

    say("")
    say("--- processing applied to each arm ---")
    for lab, steps in ops.items():
        say("  %-8s %s" % (lab, " -> ".join(steps)))
    a, b = ops["2024-25"], ops["2017"]
    same = [s.split("(")[0] for s in a] == [s.split("(")[0] for s in b]
    say("  IDENTICAL OPERATION SEQUENCE: %s" % same)
    if not same:
        say("  ABORTING: the arms differ, which is the error behind five withdrawn")
        say("  claims on 2026-08-14. No result is reported.")
        Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
        return

    say("")
    say("--- apparent-velocity budget: energy %% (and %% of available cells) ---")
    r24, r17 = budget(P24, k24, f24), budget(P17, k17, f17)
    say("  %-18s %18s %18s" % ("velocity band", "2017", "2024-25"))
    for (lo, hi, e17, c17), (_, _, e24, c24) in zip(r17, r24):
        lbl = "%5.0f - %5.0f" % (lo, hi) if np.isfinite(hi) else "%5.0f +  (k=0)" % lo
        say("  %-18s %7.2f%% (%5.2f%%) %7.2f%% (%5.2f%%)" % (lbl, e17, c17, e24, c24))

    def frac(rows, lo, hi):
        return sum(e for a_, b_, e, c in rows if a_ >= lo and b_ <= hi)
    nm17, nm24 = frac(r17, 10000, np.inf), frac(r24, 10000, np.inf)
    fan17, fan24 = frac(r17, 2500, 4000), frac(r24, 2500, 4000)
    say("")
    say("  NO ALONG-FIBRE MOVEOUT (>10 km/s, incl. k=0):  2017 %6.2f%%   2024-25 %6.2f%%"
        % (nm17, nm24))
    say("  BODY-WAVE FAN (2.5-4 km/s):                    2017 %6.2f%%   2024-25 %6.2f%%"
        % (fan17, fan24))
    say("  ratio no-moveout / fan:                        2017 %6.2f    2024-25 %6.2f"
        % (nm17 / max(fan17, 1e-9), nm24 / max(fan24, 1e-9)))
    say("")
    if nm24 > 1.5 * nm17 and (nm24 / max(fan24, 1e-9)) > 1.5 * (nm17 / max(fan17, 1e-9)):
        say("  SUPPORTS the Figure 7d reading: the 2024-25 field carries markedly more")
        say("  of its 5-20 Hz energy with NO along-fibre moveout, and markedly less of")
        say("  it in the body-wave fan, relative to 2017 under identical processing.")
    elif nm24 < nm17:
        say("  REFUTES the Figure 7d reading: 2024-25 has LESS no-moveout energy than")
        say("  2017, so the zero-lag correlation peak is not explained by the raw")
        say("  wavefield composition and another mechanism is required.")
    else:
        say("  INCONCLUSIVE: the difference is not large enough to carry the claim.")
    say("")
    say("  LIMITS: 2017 is 2 records and ~5 s total, so its budget is an estimate from")
    say("  very little data; the 2017 pre-event window of a triggered record is not")
    say("  the same population as continuous ambient; and the gauge lengths differ")
    say("  (10 m vs 16.335 m), though at 3200 m/s that is a 0.1-1.1 %% effect.")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
    lbls = ["%.0f-%.0f" % (a_, b_) if np.isfinite(b_) else ">%.0f" % a_
            for a_, b_, _, _ in r17]
    yy = np.arange(len(lbls))
    ax[0].barh(yy - 0.2, [e for _, _, e, _ in r17], 0.4, label="2017", color="crimson")
    ax[0].barh(yy + 0.2, [e for _, _, e, _ in r24], 0.4, label="2024-25", color="steelblue")
    ax[0].plot([c for _, _, _, c in r17], yy, "k.", ms=6, label="chance (cell share)")
    ax[0].set_yticks(yy); ax[0].set_yticklabels(lbls, fontsize=8)
    ax[0].set_xlabel("% of 5-20 Hz energy"); ax[0].legend(fontsize=8)
    ax[0].set_title("apparent-velocity budget, k=0 included")
    for P, k, f, lab, c in ((P17, k17, f17, "2017", "crimson"),
                            (P24, k24, f24, "2024-25", "steelblue")):
        pos = f >= 0
        prof = P[:, pos].sum(axis=1); prof /= prof.max()
        ax[1].semilogy(k, prof, color=c, lw=1.2, label=lab)
    ax[1].axvline(0, color="k", lw=0.6)
    for v, s in ((3200, "-"), (-3200, "-")):
        ax[1].axvline(10.0 / v, color="green", ls=s, lw=0.8)
    ax[1].set_xlim(-0.05, 0.05); ax[1].set_xlabel("wavenumber (cycles/m)")
    ax[1].set_ylabel("normalised power"); ax[1].legend(fontsize=8)
    ax[1].set_title("wavenumber profile (green: 3200 m/s at 10 Hz)")
    fig.savefig(str(STEM) + ".png", dpi=190)
    np.savez(str(STEM) + ".npz", rows_2017=np.array(r17, dtype=object),
             rows_2024=np.array(r24, dtype=object), k24=k24, f24=f24,
             k17=k17, f17=f17, P24=P24, P17=P17)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
