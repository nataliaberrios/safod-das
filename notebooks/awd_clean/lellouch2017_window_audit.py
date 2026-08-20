#!/usr/bin/env python3
"""Where is the earthquake in M1p33 / M2p46, and is the census's 2017 arm noise?

WHY THIS MATTERS.  ambient_apparent_velocity_census.py builds its 2017 arm as

    v[:, :int(2.5 * 250.0)]          # "pre-event only"

i.e. the first 2.5 s of each 5 s earthquake record.  The comment asserts that
window is pre-event; nothing verifies it.  The whole cross-epoch mechanism claim
("66 % of 2024-25 band energy at k = 0 versus 0.17 % in 2017") depends on it:

  - if the first 2.5 s really is quiet, the 2017 arm is a legitimate ambient
    proxy and the comparison stands (subject to its other stated limits: 5 s of
    data, and a gauge-length difference of 10 m vs 16.335 m);
  - if the arrival is inside that window, the 2017 arm contains a loud coherent
    propagating wave and the 2024-25 arm does not, so the comparison is
    apples-to-oranges and the mechanism claim must be withdrawn.  A loud arrival
    swamps instrumental common mode, which would produce exactly the near-zero
    k = 0 share observed, with no interrogator difference required.

For a local microearthquake recorded at SAFOD the source is close and triggered
records commonly place the arrival early, so this is not a remote possibility.

WHAT IT MEASURES
  1. Broadband RMS envelope versus time, per record, to locate the arrival.
  2. An automatic arrival pick: the first time the short-term RMS exceeds
     `TRIGGER` x the median RMS of the quietest quarter of the record.
  3. The k = 0 energy share computed SEPARATELY for the pre-arrival and
     post-arrival portions of the same record.  This is the direct, within-2017
     control on the mechanism claim: if a loud arrival suppresses the k = 0
     share within a single record, then any cross-epoch k = 0 comparison between
     a window containing an arrival and one that does not is measuring signal
     presence, not instrument.

No 2024-25 data and no filters are involved; this audits one assumption.

Output: lellouch2017_window_audit.{png,txt}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, detrend, sosfiltfilt

HERE = Path(__file__).resolve().parent
STEM = HERE / "lellouch2017_window_audit"
LEL = Path("/scratch/users/nberrios/lellouch2017")

FS = 250.0
NCH, NT = 800, 1250
DX = 1.0
FMIN, FMAX = 5.0, 20.0
CENSUS_WINDOW_S = 2.5      # what the census calls "pre-event only"
STA_S = 0.10               # short-term average length for the picker
TRIGGER = 3.0              # x quiet-quarter median RMS


def k0_share(x, fs, dx, fmin=FMIN, fmax=FMAX):
    """Fraction of in-band 2-D power at exactly k = 0, same recipe as the census."""
    if x.shape[1] < int(0.5 * fs):
        return np.nan, 0
    nw = int(min(2.0, x.shape[1] / fs) * fs)
    x = detrend(np.diff(x, axis=1), axis=1)
    x = sosfiltfilt(butter(4, [fmin, fmax], btype="bandpass", fs=fs, output="sos"),
                    x, axis=1)
    wt, wx = np.hanning(nw)[None, :], np.hanning(x.shape[0])[:, None]
    P = None
    n = 0
    for s in range(0, x.shape[1] - nw + 1, nw):
        F = np.fft.fftshift(np.fft.fft2(x[:, s:s + nw] * wt * wx))
        P = np.abs(F) ** 2 if P is None else P + np.abs(F) ** 2
        n += 1
    if P is None:
        return np.nan, 0
    k = np.fft.fftshift(np.fft.fftfreq(x.shape[0], dx))
    f = np.fft.fftshift(np.fft.fftfreq(nw, 1.0 / fs))
    K, F2 = np.meshgrid(k, f, indexing="ij")
    band = (np.abs(F2) >= fmin) & (np.abs(F2) <= fmax)
    at_zero = band & (K == 0.0)
    return 100.0 * P[at_zero].sum() / P[band].sum(), n


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)

    say("Audit: is the census's 2017 arm actually pre-event noise?")
    say("  census takes the first %.1f s of each 5 s record and calls it pre-event"
        % CENSUS_WINDOW_S)
    say("")

    fig, ax = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    any_found = False
    verdict_rows = []

    for row, name in enumerate(("M1p33", "M2p46")):
        p = LEL / name
        if not p.is_file():
            say("  missing: %s" % p)
            continue
        any_found = True
        v = np.fromfile(p, dtype="<f4").reshape(NT, NCH).T      # (channels, time)
        t = np.arange(NT) / FS

        # broadband RMS across channels, on the strain-rate-like difference
        d = detrend(np.diff(v, axis=1), axis=1)
        d = sosfiltfilt(butter(4, [FMIN, FMAX], btype="bandpass", fs=FS, output="sos"),
                        d, axis=1)
        rms = np.sqrt(np.mean(d ** 2, axis=0))
        nsta = max(1, int(STA_S * FS))
        sta = np.convolve(rms, np.ones(nsta) / nsta, mode="same")
        quarter = len(sta) // 4
        quiet = float(np.median(np.sort(sta)[:quarter]))
        above = np.flatnonzero(sta > TRIGGER * quiet)
        t_pick = float(t[above[0]]) if above.size else np.nan

        say("--- %s ---" % name)
        say("  record %.2f s, %d channels at %.0f m spacing, %.0f Hz"
            % (NT / FS, NCH, DX, FS))
        say("  quiet-quarter median STA %.4g | trigger at %.1fx" % (quiet, TRIGGER))
        if np.isnan(t_pick):
            say("  NO arrival picked -- record looks quiet throughout")
        else:
            say("  first arrival picked at t = %.3f s" % t_pick)
            inside = t_pick < CENSUS_WINDOW_S
            say("  arrival inside the census's %.1f s window: %s"
                % (CENSUS_WINDOW_S, "YES -- the window is NOT pre-event" if inside
                   else "no -- the window is genuinely pre-event"))

        # k=0 share before vs after the pick: the within-2017 control
        if not np.isnan(t_pick) and t_pick > 0.6:
            i = int(t_pick * FS)
            pre, npre = k0_share(v[:, :i], FS, DX)
            post, npost = k0_share(v[:, i:], FS, DX)
            say("  k=0 share, pre-arrival  (%.2f s, %d win): %s"
                % (i / FS, npre, "n/a" if np.isnan(pre) else "%.2f %%" % pre))
            say("  k=0 share, post-arrival (%.2f s, %d win): %s"
                % ((NT - i) / FS, npost, "n/a" if np.isnan(post) else "%.2f %%" % post))
            if not (np.isnan(pre) or np.isnan(post)) and post > 0:
                say("  ratio pre/post = %.1fx -- a loud arrival %s the k=0 share"
                    % (pre / post, "SUPPRESSES" if pre > post else "does not suppress"))
            verdict_rows.append((name, t_pick, pre, post))
        else:
            cen, ncen = k0_share(v[:, :int(CENSUS_WINDOW_S * FS)], FS, DX)
            say("  k=0 share over the census window (%d win): %s"
                % (ncen, "n/a" if np.isnan(cen) else "%.2f %%" % cen))
            verdict_rows.append((name, t_pick, np.nan, np.nan))
        say("")

        ax[row, 0].plot(t[: len(sta)], sta, "k-", lw=0.9)
        ax[row, 0].axhline(TRIGGER * quiet, color="crimson", ls="--", lw=1,
                           label="%.0fx quiet" % TRIGGER)
        ax[row, 0].axvspan(0, CENSUS_WINDOW_S, color="steelblue", alpha=.15,
                           label="census 2017 window")
        if not np.isnan(t_pick):
            ax[row, 0].axvline(t_pick, color="darkorange", lw=1.6,
                               label="pick %.2f s" % t_pick)
        ax[row, 0].set(xlabel="time (s)", ylabel="STA of cross-channel RMS",
                       title="%s: where is the arrival?" % name)
        ax[row, 0].legend(fontsize=7); ax[row, 0].grid(alpha=.3)

        lim = float(np.percentile(np.abs(d), 99.0))
        ax[row, 1].imshow(d, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                          extent=[t[0], t[-2], NCH, 0], interpolation="nearest")
        ax[row, 1].axvline(CENSUS_WINDOW_S, color="k", lw=1.4)
        if not np.isnan(t_pick):
            ax[row, 1].axvline(t_pick, color="darkorange", lw=1.4)
        ax[row, 1].set(xlabel="time (s)", ylabel="channel",
                       title="%s: 5-20 Hz record (black = census cut)" % name)

    if not any_found:
        raise SystemExit("2017 release not found at %s" % LEL)

    picks = [r[1] for r in verdict_rows if not np.isnan(r[1])]
    say("=== verdict ===")
    if not picks:
        say("  No arrival picked in either record: the census's 2017 arm is")
        say("  plausibly noise, and the cross-epoch comparison is NOT invalidated")
        say("  by signal presence. Its other stated limits still apply (5 s of")
        say("  data; gauge length 10 m vs 16.335 m).")
    elif all(p < CENSUS_WINDOW_S for p in picks):
        say("  Every picked arrival falls INSIDE the census's 2017 window, so that")
        say("  arm contains a loud coherent arrival while the 2024-25 arm is")
        say("  ambient noise. The cross-epoch k=0 mechanism claim is WITHDRAWN:")
        say("  it measures signal presence, not interrogator behaviour.")
    elif not any(p < CENSUS_WINDOW_S for p in picks):
        # The clean case, and the one that actually occurred: picks at 4.100 s
        # and 4.916 s of 5.00 s, both AFTER the 2.5 s cut.  The first branch
        # structure lacked this case and mislabelled it "Mixed: 0 of 2".
        say("  NO picked arrival falls inside the census's 2017 window -- every")
        say("  arrival is later than the %.1f s cut, so that arm IS pre-event"
            % CENSUS_WINDOW_S)
        say("  noise. Both arms of the census are therefore noise and the")
        say("  cross-epoch comparison is NOT confounded by signal presence.")
        say("  Its other stated limits still apply: only ~5 s of 2017 data, and")
        say("  gauge length 10 m vs 16.335 m (a 0.1-1.1 %% effect at 3200 m/s).")
    else:
        say("  Mixed: %d of %d picks fall inside the census window. Do not rely on"
            % (sum(p < CENSUS_WINDOW_S for p in picks), len(picks)))
        say("  the cross-epoch comparison; report the per-record picks instead.")
    fig.savefig(str(STEM) + ".png", dpi=190)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
