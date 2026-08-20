#!/usr/bin/env python3
"""Interrogator blame and surface illumination, both measured correctly.

v1 (`interrogator_blame_test.py`, `ambient_directional_asymmetry.py`) produced two
striking numbers that each carried a trap. Both traps are fixed here and both v1
readings are superseded.

TRAP 1 -- v1's T1 characterised a pattern the pipeline already throws away.
v1 read channels 0-699 to include the surface lead-in, and found the leading
spatial pattern u1 stable across a year (|corr(u1,u1)| median 0.84 against a
control 95th of 0.39), concluding "instrumental". But its own T2 showed u1's
lead-in/deep power ratio was 2e4 to 1.3e5: u1 lives almost entirely in channels
0-22. The Figure 7c pipeline starts at channel 23 (CH_LO_2024, the G0 wellhead),
so that stable instrumental pattern is ALREADY EXCLUDED and says nothing about
the contaminant inside the analysed aperture. Here T1 is computed on channels
>= 23 only, which is the population the conclusion needs to be about.

TRAP 2 -- v1's illumination test was forced to zero by the contaminant.
v1 measured |A| = (E(+k)-E(-k))/(E(+k)+E(-k)) in the body-wave fan and got
0.0001-0.0009 for 2024-25, concluding no surface illumination. That value is too
exact to be physical: a genuinely balanced field with 26 fan cells scatters by
~1/sqrt(26) ~ 0.2, not 1e-4. The cause is algebraic. A separable field
x(ch,t) = a(ch) s(t) has FFT2 = A(k) S(f), so P = |A(k)|^2 |S(f)|^2, and for real
a(ch) |A(-k)| = |A(k)| EXACTLY -- perfect k-symmetry, hence |A| = 0 regardless of
illumination. That separable static pattern is precisely what
AMBIENT_LOWK_MECHANISM.md found dominating the field, so v1 measured the
pattern's symmetry, not the wavefield's balance. Here |A| is measured AFTER
projecting out the leading spatial subspace, and the rank swept so the reader can
see the answer emerge (or not) as the pattern is removed.

So v1's "no surface illumination" verdict is WITHDRAWN pending this script, and
v1's "instrumental" verdict is RESTRICTED to the lead-in, where it stands.

CONTROLS
  - T1 keeps v1's control: |corr(u1 of day i, u_r of day j)| for r >= 2 bounds
    what unrelated spatial patterns score on this data.
  - |A| keeps v1's null: fan cell powers shuffled between the +k and -k halves,
    preserving total fan energy while destroying any directional preference.
  - The 2017 pre-event noise is the positive control for illumination, since
    Lellouch et al. (2019) reported downgoing/upgoing asymmetry in this data.
    Verified genuine pre-event by lellouch2017_window_audit.py (arrivals at
    4.100 s and 4.916 s of 5.00 s, both after the 2.5 s cut).
  - Rank-0 is included in the sweep so the v1 number is reproduced in place.

Reads raw HDF5 with h5py and the 2017 binaries directly. Never
DASutils.readFile_HDF, whose median=True default would delete the quantity under
test (the cross_epoch_noise_floor.py confound).

Output: interrogator_and_illumination_v2.{png,txt}
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, detrend, resample_poly, sosfiltfilt

HERE = Path(__file__).resolve().parent
STEM = HERE / "interrogator_and_illumination_v2"
CSV = Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv")
LEL = Path("/scratch/users/nberrios/lellouch2017")

FMIN, FMAX = 5.0, 20.0
FS_COMMON = 250.0
WELLHEAD_CH = 23
APERTURE_CH = 686          # 700 m at dx = 1.0209 m, as the census uses
SECONDS = 30.0
N_DAYS = 6
WIN_S = 2.0
EDIT_LO, EDIT_HI = 0.2, 5.0
V_LO, V_HI = 2500.0, 4000.0
RANKS = (0, 1, 2, 4, 8)
NULL_COUNT = 400
SEED = 20260819


def corrected_path(p):
    return Path(str(p).replace("/data/SAFODAS1-harddrive-transfer",
                               "/data/SAFOD/SAFODAS1-harddrive-transfer"))


def edit_traces(x):
    rms = np.sqrt(np.mean(x ** 2, axis=1))
    med = np.median(rms[rms > 0]) if np.any(rms > 0) else 0.0
    bad = (rms < EDIT_LO * med) | (rms > EDIT_HI * med) | ~np.isfinite(rms)
    good = np.flatnonzero(~bad)
    if good.size < 8:
        return x, -1
    for i in np.flatnonzero(bad):
        lo, hi = good[good < i], good[good > i]
        if lo.size and hi.size:
            a, b = lo[-1], hi[0]
            w = (i - a) / (b - a)
            x[i] = (1 - w) * x[a] + w * x[b]
        else:
            x[i] = x[good[0] if not lo.size else good[-1]]
    return x, int(bad.sum())


def preprocess(x, fs):
    x, dropped = edit_traces(np.asarray(x, dtype=np.float64))
    x = np.diff(x, axis=1)
    if abs(fs - FS_COMMON) > 1e-6:
        factor = int(round(fs / FS_COMMON))
        if factor >= 1 and abs(fs / factor - FS_COMMON) < 1e-6:
            x = resample_poly(x, 1, factor, axis=1)
            fs = FS_COMMON
    x = detrend(x, axis=1)
    x = sosfiltfilt(butter(4, [FMIN, FMAX], btype="bandpass", fs=fs, output="sos"),
                    x, axis=1)
    return x, fs, dropped


def remove_rank(x, rank):
    """Project out the leading `rank` left singular vectors (spatial patterns)."""
    if rank <= 0:
        return x
    xc = x - x.mean(axis=1, keepdims=True)
    u, _, _ = np.linalg.svd(xc, full_matrices=False)
    u = u[:, :rank]
    return xc - u @ (u.T @ xc)


def spectrum(x, fs, dx):
    nw = int(WIN_S * fs)
    nch = x.shape[0]
    if x.shape[1] < nw:
        return None
    wt, wx = np.hanning(nw)[None, :], np.hanning(nch)[:, None]
    P = None
    for s in range(0, x.shape[1] - nw + 1, nw):
        F = np.fft.fftshift(np.fft.fft2(x[:, s:s + nw] * wt * wx))
        P = np.abs(F) ** 2 if P is None else P + np.abs(F) ** 2
    k = np.fft.fftshift(np.fft.fftfreq(nch, dx))
    f = np.fft.fftshift(np.fft.fftfreq(nw, 1.0 / fs))
    return P, k, f


def fan_masks(k, f):
    K, F = np.meshgrid(k, f, indexing="ij")
    pos_f = (F >= FMIN) & (F <= FMAX)        # positive frequencies ONLY
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(np.abs(K) > 0, F / np.abs(K), np.inf)
    fan = pos_f & (v >= V_LO) & (v <= V_HI)
    return fan & (K > 0), fan & (K < 0), pos_f


def asymmetry(P, k, f, rng):
    up, dn, pos_f = fan_masks(k, f)
    if up.sum() < 4 or dn.sum() < 4:
        return None
    a, b = float(P[up].sum()), float(P[dn].sum())
    if a + b <= 0:
        return None
    obs = abs(a - b) / (a + b)
    pool = np.concatenate([P[up].ravel(), P[dn].ravel()])
    n_up = int(up.sum())
    nulls = np.empty(NULL_COUNT)
    for i in range(NULL_COUNT):
        s = rng.permutation(pool)
        x_, y_ = float(s[:n_up].sum()), float(s[n_up:].sum())
        nulls[i] = abs(x_ - y_) / (x_ + y_) if (x_ + y_) > 0 else 0.0
    return dict(asym=obs, null95=float(np.percentile(nulls, 95)),
                p=float((np.sum(nulls >= obs) + 1) / (NULL_COUNT + 1)),
                fan_share=100.0 * (a + b) / float(P[pos_f].sum()),
                cells=int(up.sum() + dn.sum()))


def main():
    log = []
    def say(s):
        print(s, flush=True); log.append(s)
    rng = np.random.default_rng(SEED)

    say("Interrogator blame and surface illumination, measured correctly (v2)")
    say("  band %g-%g Hz | %.0f Hz | channels %d-%d (lead-in EXCLUDED, as the"
        % (FMIN, FMAX, FS_COMMON, WELLHEAD_CH, WELLHEAD_CH + APERTURE_CH - 1))
    say("  Figure 7c pipeline does) | %.0f s/day | fan %.0f-%.0f m/s"
        % (SECONDS, V_LO, V_HI))
    say("")
    say("  v1 is superseded on both counts. Its 'instrumental' verdict is")
    say("  RESTRICTED to the surface lead-in (channels 0-22), where u1's")
    say("  lead-in/deep power ratio was 2e4-1.3e5 -- a real and stable pattern")
    say("  that the pipeline already excludes. Its 'no illumination' verdict is")
    say("  WITHDRAWN: |A| = 1e-4 was the algebraic k-symmetry of a separable")
    say("  static pattern, not a property of the wavefield.")
    say("")

    db = pd.read_csv(CSV, sep=r"\s+")
    db = db[db.nSamples == 30000].copy()
    db["t"] = pd.to_datetime(db.startTime, utc=True, errors="coerce")
    db = db.dropna(subset=["t"]).sort_values("t").reset_index(drop=True)
    db["day"] = db.t.dt.strftime("%Y-%m-%d")
    days = sorted(db.day.unique())
    picks = sorted(set(days[int(round(i * (len(days) - 1) / (N_DAYS - 1)))]
                       for i in range(N_DAYS)))

    loaded = []
    for d in picks:
        sub = db[db.day == d]
        row = sub.iloc[len(sub) // 2]
        f = corrected_path(row.file)
        if not f.is_file():
            continue
        with h5py.File(f, "r") as h:
            g = h["Acquisition/Raw[0]"]
            fs = float(g.attrs.get("OutputDataRate", 500.0))
            dx = float(h["Acquisition"].attrs.get("SpatialSamplingInterval", 1.0))
            n = int(min(SECONDS * fs, g["RawData"].shape[0]))
            hi = WELLHEAD_CH + APERTURE_CH
            x = g["RawData"][:n, WELLHEAD_CH:hi].astype(np.float32).T
        x, fs2, dropped = preprocess(x, fs)
        xc = x - x.mean(axis=1, keepdims=True)
        u, s, _ = np.linalg.svd(xc, full_matrices=False)
        frac = (s ** 2) / float(np.sum(s ** 2))
        loaded.append(dict(day=d, x=x, fs=fs2, dx=dx, u=u[:, :6], frac=frac[:6],
                           dropped=dropped, file=f.name))
        say("  %s  u1 %5.2f %%  u1..u4 %5.2f %%  (%d channels edited)"
            % (d, 100 * frac[0], 100 * frac[:4].sum(), dropped))
    if len(loaded) < 2:
        raise SystemExit("fewer than two days loaded")
    say("")

    # ---------------- T1 on the ANALYSED aperture ----------------
    say("=== T1  is the dominant spatial pattern of channels >= %d stable? ==="
        % WELLHEAD_CH)
    same, ctrl = [], []
    n = len(loaded)
    for i in range(n):
        for j in range(i + 1, n):
            c = abs(float(np.corrcoef(loaded[i]["u"][:, 0], loaded[j]["u"][:, 0])[0, 1]))
            same.append(c)
            say("  |corr(u1)| %s vs %s = %.4f" % (loaded[i]["day"], loaded[j]["day"], c))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            for r in range(1, loaded[j]["u"].shape[1]):
                ctrl.append(abs(float(np.corrcoef(
                    loaded[i]["u"][:, 0], loaded[j]["u"][:, r])[0, 1])))
    same, ctrl = np.array(same), np.array(ctrl)
    c95 = float(np.percentile(ctrl, 95))
    say("")
    say("  cross-day |corr(u1,u1)| : median %.4f  min %.4f  max %.4f  (n=%d)"
        % (np.median(same), same.min(), same.max(), same.size))
    say("  control   |corr(u1,u2+)|: median %.4f  95th %.4f            (n=%d)"
        % (np.median(ctrl), c95, ctrl.size))
    stable = float(np.median(same)) > max(0.7, c95)
    random_like = float(np.median(same)) <= c95
    say("  -> %s" % ("STABLE inside the analysed aperture: an instrumental or "
                     "fixed-installation response is the leading explanation"
                     if stable else
                     ("NOT distinguishable from unrelated patterns: a fixed "
                      "instrumental fingerprint is NOT supported here"
                      if random_like else "partially stable: does not decide")))
    say("")

    # ---------------- illumination vs rank ----------------
    say("=== illumination: |A| in the fan as the static pattern is removed ===")
    say("  rank 0 reproduces the v1 number in place. If |A| stays at ~1e-4 as")
    say("  rank rises, the field really is balanced; if it lifts above the null,")
    say("  the v1 zero was the separable pattern's symmetry.")
    say("")
    rows = []
    rec = loaded[len(loaded) // 2]
    say("  --- 2024-25 (%s) ---" % rec["day"])
    for r in RANKS:
        sp = spectrum(remove_rank(rec["x"], r), rec["fs"], rec["dx"])
        if sp is None:
            continue
        res = asymmetry(sp[0], sp[1], sp[2], rng)
        if res is None:
            continue
        rows.append(("2024-25", r, res))
        say("    rank %-2d removed: |A| = %.4f  (null95 %.4f, p = %.4f)  fan %.2f%%"
            % (r, res["asym"], res["null95"], res["p"], res["fan_share"]))

    parts = []
    for name in ("M1p33", "M2p46"):
        p = LEL / name
        if p.is_file():
            v = np.fromfile(p, dtype="<f4").reshape(1250, 800).T.astype(np.float64)
            parts.append(v[:, :int(2.5 * 250.0)])       # pre-event, per the audit
    if parts:
        x17, fs17, _ = preprocess(np.concatenate(parts, axis=1)[:APERTURE_CH], 250.0)
        say("  --- 2017 pre-event (positive control) ---")
        for r in RANKS:
            sp = spectrum(remove_rank(x17, r), fs17, 1.0)
            if sp is None:
                continue
            res = asymmetry(sp[0], sp[1], sp[2], rng)
            if res is None:
                continue
            rows.append(("2017", r, res))
            say("    rank %-2d removed: |A| = %.4f  (null95 %.4f, p = %.4f)  fan %.2f%%"
                % (r, res["asym"], res["null95"], res["p"], res["fan_share"]))
    say("")

    say("=== reading ===")
    got24 = [(r, d) for e, r, d in rows if e == "2024-25"]
    got17 = [(r, d) for e, r, d in rows if e == "2017"]
    sig24 = [r for r, d in got24 if d["p"] < 0.05]
    sig17 = [r for r, d in got17 if d["p"] < 0.05]
    say("  2024-25: significant asymmetry at ranks %s"
        % (sig24 if sig24 else "NONE of %s" % list(RANKS)))
    say("  2017   : significant asymmetry at ranks %s"
        % (sig17 if sig17 else "NONE of %s" % list(RANKS)))
    say("")
    if not got17:
        say("  No 2017 arm: the positive control is missing, so a 2024-25 null")
        say("  cannot be interpreted. No conclusion on illumination.")
    elif not sig17:
        say("  THE POSITIVE CONTROL FAILED. Lellouch reported downgoing/upgoing")
        say("  asymmetry in this 2017 data, and this measurement does not find it")
        say("  at any rank. The measurement therefore lacks the sensitivity to")
        say("  decide -- most likely because only ~5 s of 2017 noise exists -- and")
        say("  NO conclusion about 2024-25 illumination may be drawn from it.")
    elif sig24:
        say("  Illumination IS present in 2024-25 once the static pattern is")
        say("  removed, so v1's 'no illumination' verdict was indeed an artefact of")
        say("  the pattern's k-symmetry. Recovery is not blocked by illumination:")
        say("  the blocker is the static pattern, and the remedy is spatial.")
    else:
        say("  With a working positive control (2017 significant at ranks %s) and" % sig17)
        say("  2024-25 not significant at ANY rank tested, the 2024-25 field")
        say("  carries no net directional preference in the fan even after the")
        say("  static pattern is projected out. No filter can create a propagation")
        say("  direction that is absent from the data, so body-wave retrieval in")
        say("  this band is not achievable on this recording. Consistent with Behm")
        say("  (2016), where 30 s sufficed under good illumination, and with our")
        say("  stacks degrading with more data (96 h coherent stack p = 0.9184).")
    say("")
    say("  LIMITS: one 30 s file per day; 2017 is ~5 s from two records; the fan")
    say("  is the frozen %.0f-%.0f m/s selection; |A| measures a preferred"
        % (V_LO, V_HI))
    say("  direction ALONG THE FIBRE, near-vertical but not exactly vertical;")
    say("  and rank-k removal takes out any genuinely plane-wave arrival that is")
    say("  itself low-rank, so high ranks are conservative against detection.")

    fig, ax = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
    for rec2 in loaded:
        ax[0].plot(np.arange(WELLHEAD_CH, WELLHEAD_CH + rec2["u"].shape[0]),
                   rec2["u"][:, 0], lw=0.8, label=rec2["day"])
    ax[0].set(xlabel="channel", ylabel="u1",
              title="T1: dominant spatial pattern, channels >= %d" % WELLHEAD_CH)
    ax[0].legend(fontsize=6); ax[0].grid(alpha=.3)
    ax[1].hist(ctrl, bins=30, color="0.7", density=True, alpha=.85,
               label="control u1 vs u2+")
    ax[1].hist(same, bins=15, color="tab:red", density=True, alpha=.75,
               label="cross-day u1 vs u1")
    ax[1].set(xlabel="|correlation|", ylabel="density", title="T1 control")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    for epoch, colour in (("2024-25", "tab:red"), ("2017", "tab:blue")):
        pts = [(r, d) for e, r, d in rows if e == epoch]
        if not pts:
            continue
        ax[2].plot([r for r, _ in pts], [d["asym"] for _, d in pts], "o-",
                   color=colour, label="%s |A|" % epoch)
        ax[2].plot([r for r, _ in pts], [d["null95"] for _, d in pts], "--",
                   color=colour, alpha=.6, label="%s null 95th" % epoch)
    ax[2].set(xlabel="spatial rank removed", ylabel="|A| in the fan", yscale="log",
              title="Illumination vs static-pattern removal")
    ax[2].legend(fontsize=7); ax[2].grid(alpha=.3)
    fig.savefig(str(STEM) + ".png", dpi=190)
    Path(str(STEM) + ".txt").write_text("\n".join(log) + "\n")
    say("")
    say("wrote %s.{png,txt}" % STEM.name)


if __name__ == "__main__":
    main()
