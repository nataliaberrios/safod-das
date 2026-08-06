"""Where each AWD drop sits inside its raw file, and what the read-time taper did to it.

Manifest-only: reads no DAS data, runs in seconds.

The issue
---------
``readFile_protobuf`` and ``readFile_HDF`` default to ``taper=0.4`` and
``tapering=True``, and apply ``scipy.signal.windows.tukey(ntRawTot, alpha=0.4)``
to the **whole file** before filtering (``DASutils.py`` around line 1910, and the
same block in the HDF reader). Callers pass one file at a time, so ``ntRawTot``
is one file: 300 s for Nano, 62 s for Deep. Tukey alpha=0.4 is flat only between
20% and 80% of the record, so the first and last 60 s of every Nano file and the
first and last 12.4 s of every Deep file are amplitude-scaled, down to zero at
the very edge.

Weight drops are timed by the survey, not by file boundaries, so a drop's taper
weight depends on nothing but where it happened to land inside a file. That is a
multiplicative amplitude error, uncorrelated with any physics, imposed before the
burst mean is taken.

Two separate effects
--------------------
**1. A per-drop scale factor.** Each drop is multiplied by whatever the window
happens to be worth at its offset. Timing and normalised cross-correlation
survive this; per-drop and per-burst amplitude, RMS, SNR and energy observables
do not, and neither does anything regressed against them.

**2. A ramp across the analysis window, which is worse.** The taper is not flat
over the 0.5 s pre / 3.0 s post window that gets cut around each drop. On Deep
that 3.5 s is 5.6% of a 62 s file, and inside the taper ramp the weight can
change by a large factor from the start of the window to the end. That is not a
scale factor -- it is a time-varying gain applied to the waveform, so it
suppresses late arrivals relative to early ones, survives per-trace noise-RMS
normalisation (the noise window sits at one end of the ramp), and mimics
attenuation with time and therefore with distance along the fiber.

Timing, moveout, apparent velocity and F-K ridge positions remain unaffected by
both effects: the weight is real, positive and zero-phase.

``nano_hierarchical_repeatability.py`` passes ``tapering=False`` explicitly and is
clean. ``paired_stack_job.py`` and ``paired_stack_job_deep_all.py`` do not, so the
canonical epoch stacks and everything derived from them inherit it. The ambient
F-K chain reads through ``h5py`` directly and never touches this code path.

Output
------
figures/awd_2026/plain_look/fig10_taper_audit.png and taper_audit.csv
"""
from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal.windows import tukey

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "awd_manifest.csv"
OUT_DIR = Path(
    "/home/groups/ettore88/nberrios/safod_das_git/notebooks/figures/awd_2026/plain_look"
)

NANO_DURATION_S = 300.0     # build_manifest.py
DEEP_DURATION_S = 62.0      # build_manifest.py
ALPHA = 0.4                 # DASutils readFile_* default `taper`
FS = 1000.0
PRE_S, POST_S = 0.5, 3.0    # the window paired_stack_job*.py cuts around each drop


def taper_weight(offset_s: np.ndarray, duration_s: float, alpha: float = ALPHA):
    """The exact Tukey weight the reader applies at each drop's file offset."""
    n = int(round(duration_s * FS))
    w = tukey(n, alpha=alpha)
    idx = np.clip((offset_s * FS).astype(int), 0, n - 1)
    return w[idx]


def load():
    rows = []
    with open(MANIFEST) as fh:
        for r in csv.DictReader(fh):
            rows.append(dict(
                burst=int(r["burst_id"]),
                t=datetime.fromisoformat(r["utc_time"]),
                nano_ok=r["nano_available"] == "1",
                deep_ok=r["deep_available"] == "1",
                nano_off=float(r["nano_offset_s"] or "nan"),
                deep_off=float(r["deep_offset_s"] or "nan"),
            ))
    return rows


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load()
    t0 = min(r["t"] for r in rows)
    hours = np.array([(r["t"] - t0).total_seconds() / 3600.0 for r in rows])
    burst = np.array([r["burst"] for r in rows])

    nano_ok = np.array([r["nano_ok"] for r in rows])
    deep_ok = np.array([r["deep_ok"] for r in rows])
    nano_off = np.array([r["nano_off"] for r in rows])
    deep_off = np.array([r["deep_off"] for r in rows])

    wn = np.where(nano_ok, taper_weight(np.nan_to_num(nano_off), NANO_DURATION_S), np.nan)
    wd = np.where(deep_ok, taper_weight(np.nan_to_num(deep_off), DEEP_DURATION_S), np.nan)

    # the ramp across the cut window: weight at window end / weight at window start
    def ramp(off, dur, ok):
        w0 = taper_weight(np.clip(np.nan_to_num(off) - PRE_S, 0, dur), dur)
        w1 = taper_weight(np.clip(np.nan_to_num(off) + POST_S, 0, dur), dur)
        r = np.divide(w1, w0, out=np.full_like(w1, np.nan), where=w0 > 1e-6)
        return np.where(ok, r, np.nan)

    rn = ramp(nano_off, NANO_DURATION_S, nano_ok)
    rd = ramp(deep_off, DEEP_DURATION_S, deep_ok)

    fig, ax = plt.subplots(4, 2, figsize=(15, 16))

    # A: the window itself, with the drops placed on it
    for col, (dur, off, ok, w, name) in enumerate([
        (NANO_DURATION_S, nano_off, nano_ok, wn, "Nano (300 s files)"),
        (DEEP_DURATION_S, deep_off, deep_ok, wd, "Deep (62 s files)"),
    ]):
        n = int(round(dur * FS))
        tt = np.arange(n) / FS
        a = ax[0, col]
        a.plot(tt, tukey(n, alpha=ALPHA), "k", lw=1.4,
               label=r"Tukey $\alpha$=0.4 applied by the reader")
        a.plot(off[ok], w[ok], "o", ms=3, alpha=0.5, color="C3",
               label="where the drops land")
        a.axvspan(0, 0.2 * dur, color="C3", alpha=0.10)
        a.axvspan(0.8 * dur, dur, color="C3", alpha=0.10)
        a.set_xlabel("seconds into the raw file")
        a.set_ylabel("amplitude weight")
        a.set_ylim(-0.05, 1.1)
        a.set_title(f"A  {name}: {100 * np.mean(w[ok] < 0.99):.0f}% of drops are attenuated",
                    fontsize=11)
        a.legend(fontsize=8, loc="lower center")
        a.grid(alpha=0.3)

        # B: histogram of the weights actually applied
        a = ax[1, col]
        a.hist(w[ok], bins=40, color="C0", alpha=0.85)
        a.axvline(1.0, color="k", lw=1.0, ls="--")
        a.set_xlabel("taper weight applied to the drop")
        a.set_ylabel("number of drops")
        a.set_title(f"B  {name}: median {np.nanmedian(w[ok]):.2f}, "
                    f"{int(np.sum(w[ok] < 0.5))} drops below 0.5", fontsize=11)
        a.grid(alpha=0.3)

    # C: the part that matters -- how it drifts through the survey
    a = ax[2, 0]
    a.plot(hours[nano_ok], wn[nano_ok], ".", ms=3, alpha=0.35, color="C0",
           label="individual drops")
    bmeans_h, bmeans_w = [], []
    for b in np.unique(burst):
        m = nano_ok & (burst == b)
        if m.sum():
            bmeans_h.append(np.mean(hours[m]))
            bmeans_w.append(np.nanmean(wn[m]))
    a.plot(bmeans_h, bmeans_w, "o-", color="C3", lw=1.4, ms=4,
           label="burst mean")
    a.set_xlabel("hours into the survey")
    a.set_ylabel("mean taper weight")
    a.set_ylim(0, 1.08)
    a.set_title("C  Nano: a smooth ~50% amplitude excursion with no physical cause",
                fontsize=11)
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    # D: what an amplitude observable inherits
    a = ax[2, 1]
    bw = np.array(bmeans_w)
    a.plot(bmeans_h, bw / np.nanmedian(bw), "o-", color="C3", lw=1.4, ms=4)
    a.axhline(1.0, color="k", lw=0.8, ls="--")
    a.set_xlabel("hours into the survey")
    a.set_ylabel("burst amplitude bias, relative to median")
    a.set_title(f"D  Multiplicative bias on any burst RMS observable\n"
                f"range {np.nanmin(bw) / np.nanmedian(bw):.2f}x to "
                f"{np.nanmax(bw) / np.nanmedian(bw):.2f}x", fontsize=11)
    a.grid(alpha=0.3)

    # E/F: the ramp across the 3.5 s analysis window -- a time-varying gain
    for col, (r, ok, name, dur) in enumerate([
        (rn, nano_ok, "Nano", NANO_DURATION_S),
        (rd, deep_ok, "Deep", DEEP_DURATION_S),
    ]):
        a = ax[3, col]
        v = r[ok & np.isfinite(r)]
        a.hist(np.log10(v[v > 0]), bins=45, color="C1", alpha=0.85)
        a.axvline(0.0, color="k", lw=1.2, ls="--", label="flat (no ramp)")
        a.set_xlabel(r"$\log_{10}$(gain at window end / gain at window start)")
        a.set_ylabel("number of drops")
        bad = int(np.sum((v > 1.2) | (v < 1 / 1.2)))
        a.set_title(
            f"{'E' if col == 0 else 'F'}  {name}: gain ramp across the "
            f"{PRE_S + POST_S:g} s cut window\n"
            f"{bad} drops ({100 * bad / v.size:.0f}%) ramp by more than 20%; "
            f"extremes {np.nanmin(v):.2f}x to {np.nanmax(v):.2f}x", fontsize=11)
        a.legend(fontsize=8)
        a.grid(alpha=0.3)

    fig.suptitle(
        "10  Read-time Tukey taper: an amplitude modulation set by file boundaries, "
        "not by the Earth", fontsize=13)
    fig.tight_layout()
    out = OUT_DIR / "fig10_taper_audit.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)

    with open(OUT_DIR / "taper_audit.csv", "w", newline="") as fh:
        wtr = csv.writer(fh)
        wtr.writerow(["burst_id", "utc_time", "hours_into_survey",
                      "nano_offset_s", "nano_taper_weight", "nano_window_ramp",
                      "deep_offset_s", "deep_taper_weight", "deep_window_ramp"])
        for i, r in enumerate(rows):
            wtr.writerow([r["burst"], r["t"].isoformat(), f"{hours[i]:.4f}",
                          f"{nano_off[i]:.3f}", f"{wn[i]:.4f}", f"{rn[i]:.4f}",
                          f"{deep_off[i]:.3f}", f"{wd[i]:.4f}", f"{rd[i]:.4f}"])

    print("--- read-time Tukey taper audit ---")
    for name, w, ok, dur in [("Nano", wn, nano_ok, NANO_DURATION_S),
                             ("Deep", wd, deep_ok, DEEP_DURATION_S)]:
        v = w[ok]
        print(f"{name}: {ok.sum()} drops, file {dur:g} s, flat zone "
              f"{0.2 * dur:.0f}-{0.8 * dur:.0f} s")
        print(f"   attenuated (<0.99): {int((v < 0.99).sum())} "
              f"({100 * (v < 0.99).mean():.1f}%)   below 0.5: {int((v < 0.5).sum())}")
        print(f"   weight median {np.nanmedian(v):.3f}  min {np.nanmin(v):.3f}")
    print(f"burst-mean Nano weight: {np.nanmin(bw):.3f} to {np.nanmax(bw):.3f} "
          f"({np.nanmax(bw) / np.nanmin(bw):.2f}x swing across the survey)")
    print(f"\ngain ramp across the {PRE_S + POST_S:g} s cut window "
          f"(end/start; 1.0 = flat):")
    for name, r, ok in [("Nano", rn, nano_ok), ("Deep", rd, deep_ok)]:
        v = r[ok & np.isfinite(r)]
        print(f"   {name}: median {np.nanmedian(v):.3f}, "
              f">20% ramp on {int(np.sum((v > 1.2) | (v < 1 / 1.2)))} drops, "
              f"range {np.nanmin(v):.3f} to {np.nanmax(v):.3f}")
    print("\nwrote", out)


if __name__ == "__main__":
    main()
