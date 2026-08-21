#!/usr/bin/env python3
"""Ambient-noise cross-correlation on the CEMENTED Nano fibre, June 2026.

WHY THIS IS THE BEST-CONFIGURED ATTEMPT IN THE PROJECT.

Three things have to line up for borehole ambient body-wave interferometry, and
this is the only dataset here where all three do.

1  COUPLING. `manuscript/sections/MANUSCRIPT_DISCUSSION.md` establishes, from the
   AWD active-source survey, that the two June 2026 fibres see different waves:
   the CEMENTED Nano fibre records a coherent mode near 2950 m/s, while the
   WIRELINE Deep fibre records one near 1547 m/s. Its reading: "A cemented fiber
   is mechanically tied to the formation and senses strain transmitted through it;
   a wireline fiber hangs in borehole fluid and couples preferentially to energy
   guided by the fluid column." Lellouch et al. (2019) worked on a fibre "in a
   0.9 mm steel tube cemented between casing strings" -- i.e. the cemented case.
   Ambient CC on the Deep fibre duly recovered ~1525-1800 m/s (the fluid mode) and
   nothing at 3200 m/s. Nano is the fibre that should be able to see a body wave,
   and ambient CC has never been run on it.

2  APERTURE. 732 channels at 1.26606202 m = 927 m instrumented extent, against
   890 m on the 2024-25 main-hole fibre and Lellouch's 800 m of usable fibre. The
   50-700 m receiver span therefore fits, with the source in the top ~180 channels.

3  ILLUMINATION. Nano coverage is 2026-06-15 22:48 to 2026-06-17 23:48 UTC. A
   nodal deployment was underway at SAFOD on 2026-06-15 -- personnel drilling,
   digging and deploying nodes -- which is anthropogenic surface activity directly
   above the hole, the near-vertical illumination that borehole interferometry
   requires and that the unattended 2024-25 archive measurably lacks.

THE PRE-DROP CONSTRAINT. The first AWD weight drop is 2026-06-16T23:47:09Z. Every
record used here ENDS before that, so this is ambient noise and not the active
source. That leaves roughly 25 h of pre-drop Nano ambient.

Because no depth registration for Nano is assumed, the source channel is scanned
over the top of the fibre and every position is reported.

Processing is the published Lellouch et al. (2019) section 4.1 sequence, and the
correlator is the one verified against `ambient_lellouch2019_exact_stack.py` in
`deep_cc_steps.py`: differentiate to strain rate, running-absolute-mean 0.1 s,
30 s windows with 15 s overlap, R+-10 neighbour sum, simple stack, 5-20 Hz
bandpass on the full correlation and only then crop to +-0.35 s.

Note R+-10 spans +-12.7 m at Nano's 1.266 m spacing, against +-10 m at Lellouch's
1 m spacing. The channel count is kept at the published value and the metric
difference is declared rather than tuned away.

Statistics: per-velocity receiver-order permutation nulls, which unlike a
max-over-grid statistic cannot be moved by the choice of scan limits. That
distinction already invalidated one earlier deep-fibre p-value.

Reads protobuf with DASutils.readFile_protobuf using median=False, so the
across-channel median is NOT removed.

Output: nano_ambient_cc.{npz,png,txt}
"""
from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import resample_poly

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
for candidate in ("/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python",
                  "/home/groups/edunham/nberrios/safod_das_git_link/DAS-utilities/python",
                  "/home/groups/ettore88/nberrios/safod-das/DAS-utilities/python"):
    if Path(candidate).is_dir() and candidate not in sys.path:
        sys.path.insert(0, candidate)
from DASutils import readFile_protobuf  # noqa: E402

import deep_cc_steps as steps  # verified primitives  # noqa: E402

STEM = HERE / "nano_ambient_cc"
NANO_DIR = Path("/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/Nano")
DX_NANO = 1.26606202
RAW_KW = dict(filter=False, median=False, detrend=False, tapering=False)
FIRST_DROP = datetime(2026, 6, 16, 23, 47, 9, tzinfo=timezone.utc)

WELLHEAD_NPZ = HERE / "nano_find_wellhead.npz"

FS_TARGET = 250.0   # = 2.5 * fmax from the reader, and Lellouch's own rate
OFFSETS_M = np.arange(50.0, 700.1, 50.0)
NEIGHBOUR = 10
WINDOW_S, STEP_S = 30.0, 15.0
V_GRID = np.arange(300.0, 6000.1, 25.0)
NULLS = 400
SEED = 20260820


def nano_time(name):
    parts = os.path.basename(name).split("_")
    try:
        return datetime.strptime(parts[1] + " " + parts[2].replace(".", ":"),
                                 "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
    except (IndexError, ValueError):
        return None


def nano_wellhead() -> int | None:
    """First downhole Nano channel, as measured by `nano_find_wellhead.py`.

    None if that diagnostic has not been run, or ran and found no clear step.
    """
    if not WELLHEAD_NPZ.is_file():
        return None
    wh = int(np.load(WELLHEAD_NPZ)["wellhead"])
    return wh if wh > 0 else None


def require_downhole(channel: int, what: str = "source") -> None:
    """Refuse a channel that `nano_find_wellhead` places in surface/air fibre.

    THE TRAP THIS CLOSES. Deep channel 98 was scanned as a virtual source, was
    surface lead-in, and returned nothing -- a null about the cable, not the
    wavefield. `nano_find_wellhead.py` ran on 2026-08-20 and put the Nano
    wellhead at channel 73 ("channels 0-72 look like SURFACE/AIR fibre and must
    not be used"), yet every Nano product written before that guard used
    channel 10. All of those are negative results, which is exactly what a
    source in the air produces, so the guard has to be a refusal and not a
    warning.
    """
    wh = nano_wellhead()
    if wh is None:
        print("WARNING: no Nano wellhead measurement in %s, so %s channel %d "
              "cannot be checked for being downhole. Run nano_find_wellhead.py."
              % (WELLHEAD_NPZ.name, what, channel), flush=True)
        return
    if channel < wh:
        raise SystemExit(
            "REFUSING to run: %s channel %d is SURFACE/AIR fibre.\n"
            "  nano_find_wellhead.py places the Nano wellhead at channel %d "
            "(%.0f m along fibre);\n"
            "  channels 0-%d must not be used. Pass %d or greater."
            % (what, channel, wh, wh * DX_NANO, wh - 1, wh))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nfiles", type=int, default=48, help="5-minute records")
    # 73 is the measured wellhead; the previous default started at 10 and 40,
    # both of which are above it, i.e. in the air.
    ap.add_argument("--sources", type=int, nargs="+", default=[80, 120, 160, 200])
    a = ap.parse_args()
    for s in a.sources:
        require_downhole(s)

    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    rng = np.random.default_rng(SEED)

    files = sorted(NANO_DIR.glob("*.pb"))
    dated = [(p, nano_time(p.name)) for p in files]
    pre = [(p, t) for p, t in dated if t is not None and t < FIRST_DROP]
    say("Ambient CC on the CEMENTED Nano fibre")
    say("  %d Nano records total, %d before the first weight drop (%s)"
        % (len(dated), len(pre), FIRST_DROP))
    if not pre:
        raise SystemExit("no pre-drop Nano records")
    use = pre[: a.nfiles]
    say("  using %.2f h  (%d x 5-minute records): %s -> %s UTC"
        % (len(use) * 300 / 3600.0, len(use), use[0][1], use[-1][1]))
    loc0 = use[0][1].astimezone(__import__("zoneinfo").ZoneInfo("America/Los_Angeles"))
    loc1 = use[-1][1].astimezone(__import__("zoneinfo").ZoneInfo("America/Los_Angeles"))
    say("  local time: %s -> %s (nodal deployment was 2026-06-15)" % (loc0, loc1))
    gap_h = (FIRST_DROP - use[-1][1]).total_seconds() / 3600.0
    say("  last record ends %.2f h before the first weight drop: ambient, not source"
        % gap_h)
    say("")

    # MEMORY. Nano is 10 kHz raw: one 5-minute file is 300 s x 10 kHz x 732 ch x 8 B
    # = 17.6 GB, so reading 48 files undecimated OOM-killed a 96 GB job. The reader
    # desamples to 2.5*fmax, so fmax=100 Hz yields 250 Hz -- comfortably above the
    # 5-20 Hz band, 4 ms lag steps, and the same rate Lellouch et al. worked at.
    # Files are read in small batches and cast to float32 as they arrive.
    batches, fs, dx = [], None, None
    BATCH = 6
    for i in range(0, len(use), BATCH):
        chunk = [str(p) for p, _ in use[i:i + BATCH]]
        arr, info = readFile_protobuf(chunk, fmin=1.0, fmax=100.0,
                                      desampling=True, **RAW_KW)
        if fs is None:
            fs = float(info["fs"]); dx = float(info.get("dx", DX_NANO))
        batches.append(np.asarray(arr, dtype=np.float32))
        del arr
        say("    batch %d/%d: %s at %.0f Hz"
            % (i // BATCH + 1, (len(use) + BATCH - 1) // BATCH,
               batches[-1].shape, fs))
    raw = np.concatenate(batches, axis=1); del batches
    say("  read %s at fs=%.0f Hz, dx=%.5f m (%.0f m extent)"
        % (raw.shape, fs, dx, raw.shape[0] * dx))
    raw = np.asarray(raw, dtype=np.float64)
    n_win, n_step = int(WINDOW_S * fs), int(STEP_S * fs)
    say("  %.2f h of data, %d-sample windows, %d-sample step"
        % (raw.shape[1] / fs / 3600.0, n_win, n_step))
    say("")

    # published pipeline, verified primitives
    rate = steps.differentiate(raw, fs)
    normed = steps.ram_normalise(rate, fs)

    results = {}
    for src in a.sources:
        centres = src + np.rint(OFFSETS_M / dx).astype(int)
        if centres.max() + NEIGHBOUR >= raw.shape[0]:
            keep = centres + NEIGHBOUR < raw.shape[0]
            centres = centres[keep]; offs = OFFSETS_M[keep]
        else:
            offs = OFFSETS_M
        if centres.size < 6:
            say("  src %d: only %d receivers fit, skipped" % (src, centres.size))
            continue
        rows_nb = [np.arange(c - NEIGHBOUR, c + NEIGHBOUR + 1) for c in centres]
        gather, lags, nw = steps.correlate(normed, src, rows_nb, fs, n_win, n_step)
        cs = steps.moveout_curve(gather, lags, offs, V_GRID)
        ac = steps.moveout_curve(gather, lags, offs, V_GRID, sign=-1.0)
        null = np.empty((NULLS, V_GRID.size))
        for i in range(NULLS):
            null[i] = steps.moveout_curve(gather[rng.permutation(len(offs))],
                                          lags, offs, V_GRID)
        thresh = np.percentile(null, 95.0, axis=0)
        pv = (np.sum(null >= cs[None, :], axis=0) + 1.0) / (NULLS + 1.0)
        k = int(np.argmax(cs))
        clears = int(np.sum(cs > thresh))
        interior = 0 < k < V_GRID.size - 1
        say("--- source channel %d (%.2f h stacked, %d receivers %.0f-%.0f m) ---"
            % (src, raw.shape[1] / fs / 3600.0, centres.size, offs[0], offs[-1]))
        say("    peak %.4f at %.0f m/s | per-velocity p = %.4f | %s"
            % (cs[k], V_GRID[k], pv[k], "interior" if interior else "AT A SCAN EDGE"))
        say("    velocities clearing their own 95th pct: %d of %d" % (clears, V_GRID.size))
        if clears:
            good = V_GRID[cs > thresh]
            say("      band %.0f-%.0f m/s, best p %.4f at %.0f m/s"
                % (good.min(), good.max(), pv.min(), V_GRID[int(np.argmin(pv))]))
        say("    causal/acausal at the peak: %.3f" % (cs[k] / ac[k] if ac[k] else np.nan))
        for v_ref in (2950.0, 3200.0):
            say("    per-velocity p at %.0f m/s: %.4f"
                % (v_ref, float(np.interp(v_ref, V_GRID, pv))))
        say("    curve ends: %.4f at %.0f, %.4f at %.0f m/s"
            % (cs[0], V_GRID[0], cs[-1], V_GRID[-1]))
        say("")
        results[src] = dict(cs=cs, ac=ac, thresh=thresh, pv=pv, gather=gather,
                            lags=lags, offsets=offs, n_windows=nw,
                            peak_v=float(V_GRID[k]), peak_p=float(pv[k]),
                            clears=clears, interior=interior)

    if not results:
        raise SystemExit("no source channel produced a gather")

    say("=== reading ===")
    hit = {s: r for s, r in results.items() if r["clears"] > 0 and r["interior"]}
    say("  source channels with an interior peak clearing its own null: %s"
        % (sorted(hit) if hit else "NONE"))
    body = {s: r for s, r in results.items()
            if r["clears"] > 0 and 2500.0 <= r["peak_v"] <= 4000.0}
    say("  ... and peaking inside the 2500-4000 m/s body-wave fan: %s"
        % (sorted(body) if body else "NONE"))
    say("")
    if body:
        say("  A BODY-WAVE-RANGE ARRIVAL IS PRESENT on the cemented Nano fibre from")
        say("  ambient noise alone. This is the configuration Lellouch et al. (2019)")
        say("  worked in -- cemented fibre, ~900 m aperture, surface activity -- and")
        say("  it is the first body-wave-range ambient result in this project.")
        say("  Remaining gates before it is a result: causal dominance, an")
        say("  input-level null, and reproduction in an independent time window.")
    elif hit:
        say("  Coherent arrivals are recovered, but none peaks in the 2500-4000 m/s")
        say("  body-wave fan. Report the velocities found; do not describe this as a")
        say("  Lellouch reproduction.")
    else:
        say("  No source channel yields an interior peak clearing its own")
        say("  per-velocity null. Ambient CC does not recover an arrival on Nano in")
        say("  this window, despite cemented coupling and documented surface")
        say("  activity -- which would make illumination insufficient on its own.")
    say("")
    say("  LIMITS: %d records = %.2f h, one pre-drop window; source channel scanned"
        % (len(use), raw.shape[1] / fs / 3600.0))
    say("  but Nano has no depth registration here; R+-10 spans +-%.1f m at this"
        % (NEIGHBOUR * dx))
    say("  spacing; %d permutations per velocity; the fan is the frozen selection." % NULLS)

    n = len(results)
    from zoneinfo import ZoneInfo as _Z
    _la0 = use[0][1].astimezone(_Z("America/Los_Angeles"))
    _la1 = use[-1][1].astimezone(_Z("America/Los_Angeles"))
    _lab = ("%s %s-%s local  (%s to %s UTC)"
            % (_la0.strftime("%Y-%m-%d"), _la0.strftime("%H:%M"),
               _la1.strftime("%H:%M"),
               use[0][1].strftime("%Y-%m-%d %H:%M"),
               use[-1][1].strftime("%Y-%m-%d %H:%M")))
    fig, ax = plt.subplots(2, n, figsize=(4.6 * n, 7.8), squeeze=False,
                           constrained_layout=True)
    fig.suptitle("Nano fibre (cemented), %.2f h stacked   |   %s"
                 % (raw.shape[1] / fs / 3600.0, _lab), fontsize=10.5, y=1.004)
    for i, (src, r) in enumerate(sorted(results.items())):
        a0 = ax[0][i]
        g = r["gather"]
        gn = g / np.maximum(np.abs(g).max(axis=1, keepdims=True), 1e-30)
        for o, row in zip(r["offsets"], gn):
            a0.plot(r["lags"], -o + row * 40.0, "k-", lw=0.6)
        a0.plot(r["offsets"] / 2950.0, -r["offsets"], "--", color="#009E73", lw=1.4,
                label="2950 m/s (Nano active)")
        a0.plot(r["offsets"] / 3200.0, -r["offsets"], ":", color="#D55E00", lw=1.4,
                label="3200 m/s (Lellouch)")
        a0.set(xlim=(-0.35, 0.35), xlabel="lag (s)",
               ylabel="offset below source (m)" if i == 0 else "",
               title="Nano src %d, %d windows" % (src, r["n_windows"]))
        a0.legend(fontsize=6)
        a1 = ax[1][i]
        a1.plot(V_GRID / 1e3, r["cs"], "k-", lw=1.6, label="causal")
        a1.plot(V_GRID / 1e3, r["ac"], "-", color="0.6", lw=1.0, label="acausal")
        a1.plot(V_GRID / 1e3, r["thresh"], "--", color="crimson", lw=1.2,
                label="per-velocity null 95th")
        a1.axvspan(2.5, 4.0, color="#009E73", alpha=.12)
        a1.set(xlabel="trial velocity (km/s)",
               ylabel="moveout score" if i == 0 else "",
               title="peak %.0f m/s, p=%.4f, %d clear"
                     % (r["peak_v"], r["peak_p"], r["clears"]))
        a1.legend(fontsize=6); a1.grid(alpha=.3)
    fig.savefig(str(STEM) + ".png", dpi=190)
    np.savez_compressed(str(STEM) + ".npz", v_grid=V_GRID, fs=fs, dx=dx,
                        n_files=len(use), hours=raw.shape[1] / fs / 3600.0,
                        t_first=str(use[0][1]), t_last=str(use[-1][1]),
                        drop_gap_hours=gap_h,
                        **{f"src{s}_{k}": v for s, r in results.items()
                           for k, v in r.items() if k != "gather"},
                        **{f"src{s}_gather": r["gather"] for s, r in results.items()})
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{npz,png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
