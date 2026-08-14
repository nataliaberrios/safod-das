"""Generate AWD_reproduce_analysis.ipynb, the end-to-end reproducibility notebook.

The notebook is the artifact; this script exists so it can be regenerated and so
its structure is reviewable as code rather than as JSON. Edit here and re-run,
or edit the notebook directly once it is being used interactively — but not both.

    python build_notebook.py
"""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "AWD_reproduce_analysis.ipynb"

cells: list[dict] = []


def _source(lines: tuple[str, ...]) -> list[str]:
    """nbformat wants each entry to carry its own newline, except the last."""
    flat = "\n".join(lines).split("\n")
    return [f"{line}\n" for line in flat[:-1]] + [flat[-1]]


def md(*lines: str) -> None:
    cells.append({"cell_type": "markdown", "metadata": {}, "source": _source(lines)})


def code(*lines: str) -> None:
    cells.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": _source(lines)})


# ===========================================================================
md(
 "# SAFOD repeated-source DAS — reproducible analysis",
 "",
 "**From burst stacks to a minimum detectable velocity change, one step at a time.**",
 "",
 "This notebook reproduces every analysis step behind the paper's central result:",
 "that two DAS installations in the same borehole resolve different smallest",
 "velocity changes — 1% on the cemented fiber, 0.5% on the wireline outbound",
 "branch — and that the difference is set by the trade between propagation lever",
 "arm and source timing repeatability.",
 "",
 "Every cell is meant to be read and run. Where a step is too expensive to redo",
 "inline, the notebook says which script does it and loads that script's frozen",
 "output instead of hiding the gap.",
 "",
 "---",
 "",
 "## What this notebook recomputes, and what it does not",
 "",
 "| Step | Here? | Script of record |",
 "|---|---|---|",
 "| Drop inventory from GPS times and file coverage | no — load `awd_manifest.csv` | `build_manifest.py` |",
 "| Cut 3.5 s windows from raw fiber, stack into bursts | no — load the stacks `.npz` | `paired_stack_job_deep_all.py` |",
 "| Inventory QC and the 49 → 46 burst question | **yes** | — |",
 "| Band-pass and channel selection | **yes** | — |",
 "| Record sections for both fibers | **yes** | — |",
 "| Semblance scan that picks the arrival | **yes**, reduced grid, checked against the frozen value | `deep_dvv_injection_recovery.py --stage freeze` |",
 "| Virtual source: correlation, deconvolution, wiggles | **yes** | `awd_virtual_source.py` |",
 "| The dv/v estimator, worked through for one burst | **yes** | `deep_dvv_injection_recovery.py` |",
 "| Injection–recovery, live demonstration | **yes**, a few trials | — |",
 "| Injection–recovery, full 1380-trial calibration | no — load results | `deep_dvv_injection_recovery.py --stage all` |",
 "",
 "The two 'no' rows at the top read raw fiber data off `$OAK` and take about",
 "20 minutes; the full calibration is 16 minutes at 96 GB. Everything else runs",
 "here in a couple of minutes.",
 "",
 "**Resources.** Needs ~8 GB of memory — the Deep stack array alone is 2.2 GB.",
 "Run it in `sh_dev --mem=16G` or a Slurm job, not on a login node.",
)

# --------------------------------------------------------------------------
md("---", "", "## 1. Environment and paths",
   "",
   "Nothing is installed for this project; `DASutils` is reached through",
   "`PYTHONPATH`. The analysis scripts live one directory up, and the notebook",
   "imports from them rather than copying code, so there is a single source of",
   "truth for every constant.")

code(
 "import sys, json, csv",
 "from pathlib import Path",
 "",
 "import numpy as np",
 "import matplotlib.pyplot as plt",
 "",
 "ANALYSIS = Path.cwd().parent if Path.cwd().name == 'manuscript' else Path('awd_clean')",
 "sys.path.insert(0, str(ANALYSIS))",
 "",
 "import deep_dvv_injection_recovery as D      # frozen constants + estimator",
 "from fk_dispersion import weighted_stack",
 "",
 "STACKS = ANALYSIS / 'canonical_epoch_stacks_paired_deep_all.npz'",
 "print('analysis dir :', ANALYSIS)",
 "print('stacks       :', STACKS.name, f'({STACKS.stat().st_size/1e9:.2f} GB)')",
 "print('numpy        :', np.__version__)",
)

md("### The frozen constants",
   "",
   "These were fixed before any sensitivity result was computed and are recorded",
   "in `DEEP_DVV_PREREGISTRATION.md`. Printing them here so nothing in the",
   "notebook can silently disagree with the analysis.")

code(
 "print(f'analysis range      {D.COORD_RANGE_M[0]:.0f}-{D.COORD_RANGE_M[1]:.0f} m along each leg')",
 "print(f'channel stride      {D.CHANNEL_STRIDE}')",
 "print(f'turnaround channel  {D.TURNAROUND_CH}')",
 "print(f'aperture            {D.APERTURE_M:.0f} m wide, {D.STEP_M:.0f} m step')",
 "print(f'primary band        {D.BAND_CONFIG[D.PRIMARY_BAND][\"band\"]} Hz')",
 "print(f'slowness search     {1/D.SLOWNESS_RANGE[1]:.0f}-{1/D.SLOWNESS_RANGE[0]:.0f} m/s')",
 "print(f'QC: correlation >=  {D.MIN_APERTURE_CORRELATION}, at least {D.MIN_APERTURES} apertures')",
 "print(f'injection levels    {len(D.INJECTED_DVV)}: {D.INJECTED_DVV[-3:]} ...')",
)

# --------------------------------------------------------------------------
md("---", "", "## 2. The data, and the burst inventory",
   "",
   "The stacks file holds one waveform per burst per channel, for both fibers,",
   "already cut and averaged by `paired_stack_job_deep_all.py`. Load the metadata",
   "first — the big arrays come later, one at a time.")

code(
 "meta = np.load(STACKS)",
 "fs        = float(meta['fs'])",
 "dx_nano   = float(meta['dx_nano'])",
 "dx_deep   = float(meta['dx_deep'])",
 "counts    = np.asarray(meta['n_common'], dtype=int)   # drops per burst, both fibers",
 "begtimes  = meta['begtimes_str'].astype(str)",
 "",
 "print(f'sample rate      {fs:.0f} Hz')",
 "print(f'Nano  {meta[\"nano_stacks\"].shape[1]:5d} channels at {dx_nano:.4f} m'",
 "      f'  -> {meta[\"nano_stacks\"].shape[1]*dx_nano:.0f} m of fiber')",
 "print(f'Deep  {meta[\"deep_stacks\"].shape[1]:5d} channels at {dx_deep:.4f} m'",
 "      f'  -> {meta[\"deep_stacks\"].shape[1]*dx_deep:.0f} m of fiber')",
 "print(f'slots {counts.size}, of which {(counts>0).sum()} have drops on both fibers')",
 "print(f'total common drops {counts.sum()}')",
)

md("### Why 46 bursts and not 49",
   "",
   "This is the first thing anyone asks. The survey fired 49 bursts, but three of",
   "them have no drops recorded on **both** fibers. They are not scattered",
   "dropouts and nothing was cut for quality — they are the last three bursts,",
   "because the Deep fiber stopped recording before the survey ended.")

code(
 "from datetime import datetime, timezone",
 "",
 "times = [datetime.fromisoformat(t).astimezone(timezone.utc) if t else None",
 "         for t in begtimes]",
 "t0 = next(t for t in times if t is not None)",
 "hours = np.array([(t - t0).total_seconds()/3600 if t else np.nan for t in times])",
 "",
 "fig, ax = plt.subplots(figsize=(11, 3.2), constrained_layout=True)",
 "ok = counts > 0",
 "ax.bar(hours[ok],  counts[ok],  width=0.32, color='#2a78d6', label='drops on both fibers')",
 "ax.bar(hours[~ok], np.full((~ok).sum(), 20), width=0.32, color='#eb6834',",
 "       label='Nano only — Deep had stopped')",
 "ax.set(xlabel='hours from first burst', ylabel='drops in burst',",
 "       title='Burst inventory: the last three bursts have no Deep data')",
 "ax.legend(frameon=False); ax.grid(alpha=.2)",
 "plt.show()",
 "",
 "for i in np.flatnonzero(~ok):",
 "    print(f'  burst {i}: {begtimes[i][:19]}  common drops = {counts[i]}')",
)

md("**The drop counts, and why four different numbers are all correct.**",
   "Reading them off the manifest rather than asserting them.")

code(
 "man = list(csv.DictReader(open(ANALYSIS / 'awd_manifest.csv')))",
 "nano_av = sum(int(r['nano_available']) for r in man)",
 "deep_av = sum(int(r['deep_available']) for r in man)",
 "print(f'rows in awd_manifest.csv          {len(man)}')",
 "print(f'  available on Nano               {nano_av}')",
 "print(f'  available on Deep               {deep_av}')",
 "print(f'used by the stacks (both fibers)  {counts.sum()}   <- the analysed dataset')",
 "print()",
 "print('The manifest-to-stack drop is window truncation at file boundaries:')",
 "print('paired_stack_job_deep_all.py keeps a drop only if its full 0.5/3.0 s')",
 "print('window fits inside the recording file.')",
)

# --------------------------------------------------------------------------
md("---", "", "## 3. Preprocessing — what has already been applied",
   "",
   "Most preprocessing happened inside the reader, with defaults, before these",
   "stacks were written. It is invisible unless you go looking, so it is stated",
   "here (full detail in `PREPROCESSING.md`):",
   "",
   "```",
   "detrend -> Tukey taper -> Butterworth 1-100 Hz -> across-channel median removal -> scaling",
   "```",
   "",
   "Three consequences worth knowing:",
   "",
   "1. **Median removal is a modelling choice, not cleanup.** Subtracting the",
   "   across-channel median removes interrogator noise common to all channels —",
   "   and equally removes any arrival that reaches every channel at once.",
   "2. **The taper spans the whole file** (300 s Nano, 60 s Deep), not the 3.5 s",
   "   drop window cut from it.",
   "3. **The two fibers are in different units.** Nano is strain *rate*, Deep is",
   "   strain. That is a factor of ω — a 90° phase rotation — which changes",
   "   waveform shape but *not* an F–K ridge or a moveout velocity. Since",
   "   everything compared here is a velocity or a delay, it does not affect the",
   "   comparison. It would affect any amplitude comparison between fibers, and",
   "   none is made.",
   "",
   "On top of that, this analysis applies its own band-pass. Both stages are",
   "zero-phase Butterworth, so the arrival is not shifted in time.")

code(
 "def bandpass(x, fs, band, order=4):",
 "    \"\"\"Zero-phase Butterworth. Same as D._bandpass; written out so the\"\"\"",
 "    \"\"\"filtering is visible rather than hidden behind an import.\"\"\"",
 "    from scipy.signal import butter, sosfiltfilt",
 "    sos = butter(order, band, btype='bandpass', fs=fs, output='sos')",
 "    return sosfiltfilt(sos, np.nan_to_num(x), axis=-1)",
 "",
 "NANO_BAND = (30.0, 60.0)",
 "DEEP_BAND = D.BAND_CONFIG[D.PRIMARY_BAND]['band']",
 "print(f'Nano analysis band {NANO_BAND} Hz')",
 "print(f'Deep analysis band {DEEP_BAND} Hz')",
)

# --------------------------------------------------------------------------
md("---", "", "## 4. What each fiber records",
   "",
   "Average the bursts into one section per fiber and look at them. This is the",
   "paper's organising observation: the same source, in the same hole, produces",
   "a different coherent arrival on each installation.")

code(
 "# Count-weighted average over bursts, the same reduction every script uses.",
 "# fk_dispersion.weighted_stack returns the section only; the drop total is",
 "# just the sum of the per-burst counts.",
 "nano_section = weighted_stack(meta['nano_stacks'], counts)",
 "n_drops = int(counts.sum())",
 "print(f'Nano section {nano_section.shape} from {n_drops} drops')",
)

code(
 "deep_all = meta['deep_stacks']                       # 2.2 GB — the big one",
 "deep_section = weighted_stack(deep_all, counts)",
 "del deep_all                                          # free it immediately",
 "print(f'Deep section {deep_section.shape}')",
)

code(
 "def record_section(ax, section, dx, band, tmax, title, decimate=4):",
 "    filt = bandpass(section[::decimate], fs, band)",
 "    n = int(tmax * fs)",
 "    filt = filt[:, :n]",
 "    # display scaling only: each trace by its own 99th percentile, so a loud",
 "    # shallow interval cannot set the colour scale for the whole section",
 "    scale = np.percentile(np.abs(filt), 99, axis=1, keepdims=True)",
 "    scale[scale == 0] = 1.0",
 "    ax.imshow(filt/scale, aspect='auto', cmap='RdBu_r', vmin=-1, vmax=1,",
 "              extent=[-D.PRE_S, tmax - D.PRE_S,",
 "                      section.shape[0]*dx, 0])",
 "    ax.set(xlabel='time from drop (s)', ylabel='distance along fiber (m)', title=title)",
 "",
 "fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)",
 "record_section(axes[0], nano_section, dx_nano, NANO_BAND, 1.0,",
 "               f'Nano, cemented — {NANO_BAND[0]:.0f}-{NANO_BAND[1]:.0f} Hz')",
 "record_section(axes[1], deep_section[:D.TURNAROUND_CH], dx_deep, DEEP_BAND, 3.0,",
 "               f'Deep outbound, wireline — {DEEP_BAND[0]:.0f}-{DEEP_BAND[1]:.0f} Hz')",
 "plt.show()",
)

md("Nano's arrival is steep and dies within a few hundred metres. Deep's is much",
   "shallower — a slower mode — and persists for kilometres. Those two slopes are",
   "the whole story: slope is 1/speed, and the length of the arrival is the lever",
   "arm the estimator gets to work with.")

# --------------------------------------------------------------------------
md("---", "", "## 5. Finding the arrival: the semblance scan",
   "",
   "An arrival travelling along the fiber is a straight line in distance–time, so",
   "describing it takes two numbers: a start time `t0` and a speed `v`.",
   "",
   "To find them: shift every channel by its predicted arrival time and sum. A",
   "correct guess makes the arrival add constructively; a wrong one cancels it.",
   "**Semblance** is that coherent energy divided by the total energy. Scan a grid",
   "of `(t0, v)` and take the peak.",
   "",
   "**This scan runs on half the bursts only** — the odd-indexed 23. The other 23",
   "are held back so sensitivity is never measured on the data that chose what to",
   "look at. That split is the backbone of the whole calibration.",
   "",
   "Below is a reduced-resolution version so it runs in the notebook; the frozen",
   "value came from the full grid in the freeze stage. The reduced scan should",
   "land on the same peak.")

code(
 "# Discovery half only: keep epoch % 2 == 1, matching deep_tube_validation.py",
 "valid = np.flatnonzero(counts > 0)",
 "discovery_counts = counts.copy()",
 "discovery_counts[valid[valid % 2 == 0]] = 0",
 "print(f'discovery bursts {np.count_nonzero(discovery_counts)}, '",
 "      f'held out {len(valid) - np.count_nonzero(discovery_counts)}')",
 "",
 "deep_all = np.load(STACKS)['deep_stacks']",
 "discovery = weighted_stack(deep_all, discovery_counts)",
 "del deep_all",
)

code(
 "channels = D._leg_channels(discovery.shape[0], dx_deep)",
 "absolute, coordinate = channels['outbound']",
 "cfg = D.BAND_CONFIG[D.PRIMARY_BAND]",
 "",
 "# Reduced grid: 41 speeds instead of 101, every 3rd channel. Full-resolution",
 "# values live in deep_dvv_frozen_trajectory.json.",
 "section  = D._rms_normalize(D._bandpass(discovery[absolute][::3], fs, cfg['band']))",
 "coord    = coordinate[::3]",
 "sample_t = np.arange(section.shape[-1]) / fs - D.PRE_S",
 "",
 "intercepts = np.arange(D.INTERCEPT_GRID_S[0], D.INTERCEPT_GRID_S[1] + 1e-9, 0.004)",
 "pad = cfg['semblance_half_s'] + 2/fs",
 "q   = np.arange(int(round(((intercepts[-1]+pad) - (intercepts[0]-pad))*fs))+1) / fs \\",
 "      + (intercepts[0] - pad)",
 "slowness = np.linspace(D.SLOWNESS_RANGE[0], D.SLOWNESS_RANGE[1], 41)",
 "",
 "grid = np.empty((slowness.size, intercepts.size))",
 "for i, p in enumerate(slowness):",
 "    aligned = D._align(section, sample_t, coord, 0.0, p, q)",
 "    grid[i] = D._semblance_profile(aligned, q, fs, cfg['semblance_half_s'], intercepts)",
 "",
 "peak = np.unravel_index(np.argmax(grid), grid.shape)",
 "v_found, t0_found = 1/slowness[peak[0]], intercepts[peak[1]]",
 "frozen = json.loads((ANALYSIS/'deep_dvv_frozen_trajectory.json').read_text())",
 "expect = frozen['trajectories']['outbound|15_30']",
 "print(f'reduced scan : {v_found:7.1f} m/s at t0 {t0_found:+.3f} s')",
 "print(f'frozen value : {expect[\"velocity_mps\"]:7.1f} m/s at t0 {expect[\"intercept_s\"]:+.3f} s')",
 "print(f'speed differs by {abs(v_found-expect[\"velocity_mps\"]):.1f} m/s')",
)

code(
 "fig, ax = plt.subplots(figsize=(8, 5.2), constrained_layout=True)",
 "im = ax.pcolormesh(intercepts, 1/slowness, grid, cmap='Blues', shading='auto')",
 "ax.plot(expect['intercept_s'], expect['velocity_mps'], 'o', ms=15,",
 "        mfc='none', mec='#eb6834', mew=2.5, label='frozen trajectory')",
 "ax.set(xlabel='start time $t_0$ (s)', ylabel='apparent speed (m/s)',",
 "       title='Semblance scan, Deep outbound (reduced grid)')",
 "ax.legend(frameon=False); plt.colorbar(im, ax=ax, label='semblance')",
 "plt.show()",
)

md("**Read the shape, not just the peak.** Speed is tightly constrained — a narrow",
   "horizontal band. Start time is *not*: the band runs across most of the searched",
   "range, modulated at the wavelet period (cycle skipping).",
   "",
   "That is fine, and it is worth understanding why. The estimator in §7 carries a",
   "free intercept that absorbs any constant timing offset, so the recovered",
   "velocity change depends on the **slowness** alone. A poorly determined `t0`",
   "costs nothing.")

# --------------------------------------------------------------------------
md("---", "", "## 6. Virtual source: correlation, deconvolution, wiggles",
   "",
   "A separate check on the same arrival, using none of the trajectory above.",
   "",
   "Every channel records the *same* source wobble, so cross-correlating channel A",
   "with channel B cancels the source term and leaves the response between them —",
   "as though a source sat at A. This **redatums** the surface weight drop down",
   "into the borehole, taking the near-surface path and the source timing jitter",
   "with it (Bakulin & Calvert 2006).",
   "",
   "Correlation removes the source's *timing* but leaves its *shape* smeared over",
   "the result. **Deconvolution** divides it out instead, so the arrival comes out",
   "sharp. Both are computed because the AWD never moved: without a spread of",
   "source positions the classic construction is unavailable, so the correlation",
   "gather is a redatumed section rather than a clean Green's function.")

code(
 "def _center(x, fs, max_lag):",
 "    n = x.shape[-1]",
 "    k = int(round(max_lag * fs))",
 "    x = np.roll(x, n//2, axis=-1)",
 "    mid = n//2",
 "    lags = np.arange(-k, k+1) / fs",
 "    return lags, x[:, mid-k:mid+k+1]",
 "",
 "def correlate_gather(section, source_trace, fs, max_lag):",
 "    \"\"\"B correlated with A: cancels the source term, keeps its spectrum.\"\"\"",
 "    n = section.shape[-1]",
 "    nfft = 1 << int(np.ceil(np.log2(2*n - 1)))",
 "    A = np.fft.rfft(source_trace, n=nfft)",
 "    B = np.fft.rfft(section, n=nfft, axis=-1)",
 "    return _center(np.fft.irfft(B * np.conj(A)[None, :], n=nfft, axis=-1), fs, max_lag)",
 "",
 "def deconvolve_gather(section, source_trace, fs, max_lag, water=0.01):",
 "    \"\"\"B deconvolved by A: cancels the source term *and* its spectrum.",
 "    Water-level regularisation guards the division where |A| is near zero.\"\"\"",
 "    n = section.shape[-1]",
 "    nfft = 1 << int(np.ceil(np.log2(2*n - 1)))",
 "    A = np.fft.rfft(source_trace, n=nfft)",
 "    B = np.fft.rfft(section, n=nfft, axis=-1)",
 "    power = np.abs(A)**2",
 "    eps = water * float(power.mean())",
 "    return _center(np.fft.irfft(B*np.conj(A)[None, :]/(power+eps)[None, :],",
 "                                n=nfft, axis=-1), fs, max_lag)",
)

md("The lag window matters and is easy to get wrong. Over Deep's 2.8 km aperture",
   "at ~1545 m/s the differential travel time reaches ~1.8 s, so a window sized for",
   "Nano's 600 m would silently truncate the gather — producing a plausible-looking",
   "figure with the arrival cut off. Sizing it from the geometry instead:")

code(
 "VS_APERTURE = (200.0, 3000.0)",
 "c0, c1 = int(VS_APERTURE[0]/dx_deep), int(VS_APERTURE[1]/dx_deep)",
 "vs_sec = bandpass(deep_section[c0:c1:3], fs, DEEP_BAND)",
 "vs_z   = np.arange(c0, c1, 3) * dx_deep",
 "",
 "span   = vs_z[-1] - vs_z[0]",
 "MAX_LAG = 1.3 * span / expect['velocity_mps']      # 30% headroom over the geometry",
 "print(f'aperture {vs_z[0]:.0f}-{vs_z[-1]:.0f} m ({vs_sec.shape[0]} channels)')",
 "print(f'span/speed = {span/expect[\"velocity_mps\"]:.2f} s  ->  max lag {MAX_LAG:.2f} s')",
 "",
 "src = int(np.argmin(np.abs(vs_z - 900.0)))         # virtual source near 900 m",
 "lags, corr = correlate_gather(vs_sec, vs_sec[src], fs, MAX_LAG)",
 "_,    deco = deconvolve_gather(vs_sec, vs_sec[src], fs, MAX_LAG)",
 "print(f'virtual source at {vs_z[src]:.0f} m')",
)

code(
 "def norm_rows(g):",
 "    peak = np.max(np.abs(g), axis=1, keepdims=True)",
 "    peak[peak == 0] = 1.0",
 "    return g / peak",
 "",
 "fig, ax = plt.subplots(1, 3, figsize=(16, 6), constrained_layout=True)",
 "for a, g, name in ((ax[0], corr, 'correlation'), (ax[1], deco, 'deconvolution')):",
 "    a.imshow(norm_rows(g), aspect='auto', cmap='RdBu_r', vmin=-1, vmax=1,",
 "             extent=[lags[0], lags[-1], vs_z[-1], vs_z[0]])",
 "    a.axhline(vs_z[src], color='#1baf7a', lw=2)",
 "    a.set(xlabel='lag (s)', ylabel='distance along fiber (m)', title=name)",
 "",
 "# wiggles: the same deconvolution data as traces, so waveform shape is visible",
 "step = max(1, deco.shape[0]//40)",
 "for i in range(0, deco.shape[0], step):",
 "    tr = deco[i] / (np.max(np.abs(deco[i])) or 1.0)",
 "    ax[2].plot(lags, vs_z[i] - tr*90, color='k', lw=0.6)",
 "ax[2].axhline(vs_z[src], color='#1baf7a', lw=2)",
 "ax[2].invert_yaxis()",
 "ax[2].set(xlabel='lag (s)', ylabel='distance along fiber (m)', title='deconvolution, wiggles')",
 "plt.show()",
)

md("The green line is the virtual source. The arrival radiating away from it is",
   "the redatumed wave: further from the source means later lag, hence the slant.",
   "Its slope is a speed — measured below, with nothing from §5 used as input.")

code(
 "# Slant-stack: sweep speeds, sum along each predicted moveout, keep the best.",
 "def slant_speed(gather, lags, z, z_src, speeds):",
 "    best, best_v = -np.inf, np.nan",
 "    for v in speeds:",
 "        pred = np.abs(z - z_src) / v",
 "        idx  = np.searchsorted(lags, pred)",
 "        ok   = (idx > 0) & (idx < lags.size)",
 "        if ok.sum() < 20:",
 "            continue",
 "        picked = gather[np.flatnonzero(ok), idx[ok]]",
 "        score  = picked.sum()**2 / (ok.sum() * (picked**2).sum() + 1e-30)",
 "        if score > best:",
 "            best, best_v = score, v",
 "    return best_v, best",
 "",
 "v_vs, semb = slant_speed(norm_rows(deco), lags, vs_z, vs_z[src],",
 "                         np.arange(1000, 2001, 25))",
 "print(f'virtual-source speed  {v_vs:.0f} m/s   (semblance {semb:.3f})')",
 "print(f'semblance-scan speed  {expect[\"velocity_mps\"]:.1f} m/s')",
 "print(f'difference            {abs(v_vs-expect[\"velocity_mps\"])/expect[\"velocity_mps\"]*100:.1f}%')",
)

md("Two estimators that share almost nothing agree. The semblance scan needs to",
   "know *when the drop happened* and fits a predicted straight line; the virtual",
   "source throws the source time away entirely and works only from differences",
   "between channels. They still land on the same speed.",
   "",
   "They are not fully independent — same stacks, same band, same geometry — so the",
   "honest claim is: two different estimators, one of which needs source timing and",
   "one of which cancels it, agree.")

# --------------------------------------------------------------------------
md("---", "", "## 7. The estimator: turning delays into a velocity change",
   "",
   "This is the measurement itself, and it is one line fit.",
   "",
   "Split the leg into overlapping apertures. In each, shift channels onto the",
   "frozen trajectory and average them into a beam — that stacking is what buys",
   "usable signal-to-noise 3 km down the fiber. Compare one burst's beam against",
   "a reference built from **all the other bursts**, and measure the delay.",
   "",
   "Then the key step:",
   "",
   "- if the **source** fired early, every aperture's delay shifts by the same amount",
   "- if the **medium** sped up, the delay grows with travel time",
   "",
   "so fit",
   "",
   "$$\\delta t = a - \\varepsilon\\,T_0$$",
   "",
   "The intercept $a$ absorbs source timing. The slope is the fractional velocity",
   "change $\\varepsilon$. Worked through below for a single burst.")

code(
 "leg = 'outbound'",
 "absolute, coordinate = channels[leg]",
 "epochs = valid[valid % 2 == 0]              # the 23 held-out bursts",
 "print(f'{len(epochs)} held-out bursts')",
 "",
 "deep_all = np.load(STACKS)['deep_stacks']",
 "raw = np.asarray(deep_all[np.ix_(epochs, absolute)], dtype=float)",
 "del deep_all",
 "",
 "keep  = D._channel_qc(raw)",
 "raw   = np.nan_to_num(raw[:, keep])",
 "coord = coordinate[keep]",
 "print(f'{raw.shape[1]} channels pass QC of {keep.size}')",
)

code(
 "# Align every burst onto the frozen trajectory, per channel.",
 "p_frozen  = expect['slowness_s_per_m']",
 "t0_frozen = expect['intercept_s']",
 "q0, q1 = cfg['extract_q_s']",
 "q_ax   = np.arange(int(round((q1-q0)*fs))) / fs + q0",
 "",
 "filtered = D._rms_normalize(D._bandpass(raw, fs, cfg['band']))",
 "sample_t = np.arange(filtered.shape[-1]) / fs - D.PRE_S",
 "aligned  = D._align(filtered, sample_t, coord, t0_frozen, p_frozen, q_ax)",
 "print(f'aligned gathers {aligned.shape}  (bursts, channels, lag samples)')",
 "",
 "# Reference for each burst = count-weighted mean of all the OTHER bursts.",
 "references = D._leave_one_out(aligned, counts[epochs])",
 "layout  = D._aperture_layout(coord, D.APERTURE_M)",
 "centres = np.array([c for _, c in layout])",
 "travel  = t0_frozen + p_frozen * centres",
 "print(f'{len(layout)} apertures, centres {centres[0]:.0f}-{centres[-1]:.0f} m')",
 "print(f'lever arm in travel time: {travel.max()-travel.min():.3f} s')",
)

code(
 "b = 0                                   # one burst, worked through",
 "target = D._beams(aligned[b],    q_ax, coord, layout, None)",
 "ref    = D._beams(references[b], q_ax, coord, layout, None)",
 "",
 "use = (q_ax >= cfg['recovery_q_s'][0]) & (q_ax <= cfg['recovery_q_s'][1])",
 "delay = np.full(len(layout), np.nan); corr_ = np.full(len(layout), np.nan)",
 "for i in range(len(layout)):",
 "    delay[i], corr_[i] = D._normalized_correlation_lag(",
 "        ref[i][use], target[i][use], fs, cfg['max_lag_s'])",
 "",
 "good = (np.abs(delay) <= 0.9*cfg['max_lag_s']) & (corr_ >= D.MIN_APERTURE_CORRELATION)",
 "beta, se, _ = D._robust_line(travel[good], delay[good], np.maximum(corr_[good], 0)**2)",
 "eps = -beta[1]",
 "print(f'burst {epochs[b]}: {good.sum()} apertures kept')",
 "print(f'  intercept a  = {beta[0]*1e3:+.3f} ms   (source timing)')",
 "print(f'  slope        = {eps:+.3e}             (velocity change)')",
)

code(
 "fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)",
 "ax.plot(travel[good], delay[good]*1e3, 'o', ms=10, color='#eb6834',",
 "        mec='white', mew=1.6, label='aperture delays')",
 "ax.plot(travel[~good], delay[~good]*1e3, 'x', ms=9, color='#9a9a95', label='failed QC')",
 "line = np.array([travel.min(), travel.max()])",
 "ax.plot(line, (beta[0] + beta[1]*line)*1e3, color='#2a78d6', lw=2.5,",
 "        label=f'fit: intercept {beta[0]*1e3:+.2f} ms, $\\\\varepsilon$ = {eps:+.2e}')",
 "ax.set(xlabel='reference travel time $T_0$ (s)', ylabel='delay vs reference (ms)',",
 "       title=f'One burst: delay against travel time (burst {epochs[b]})')",
 "ax.axhline(0, color='0.4', lw=.8); ax.legend(frameon=False); ax.grid(alpha=.2)",
 "plt.show()",
)

# --------------------------------------------------------------------------
md("---", "", "## 8. Minimum detectable change: injection and recovery",
   "",
   "How small a velocity change can that estimator see? Rather than derive it,",
   "**inject a known change into the real data and try to recover it.**",
   "",
   "For an injected fraction ε, shift each channel by −ε·T₀(s) — exactly as a real",
   "velocity change would — then run the whole chain and see what comes back. In",
   "the production run the recovery code never sees the injected value: injection",
   "writes perturbed gathers plus a sealed key, and the two are joined only at the",
   "end.",
   "",
   "A live demonstration on a handful of levels first.")

code(
 "def recover(gather, reference, layout, travel, q_ax, cfg):",
 "    tgt = D._beams(gather,    q_ax, coord, layout, None)",
 "    ref = D._beams(reference, q_ax, coord, layout, None)",
 "    use = (q_ax >= cfg['recovery_q_s'][0]) & (q_ax <= cfg['recovery_q_s'][1])",
 "    d = np.full(len(layout), np.nan); c = np.full(len(layout), np.nan)",
 "    for i in range(len(layout)):",
 "        d[i], c[i] = D._normalized_correlation_lag(ref[i][use], tgt[i][use],",
 "                                                   fs, cfg['max_lag_s'])",
 "    ok = (np.abs(d) <= 0.9*cfg['max_lag_s']) & (c >= D.MIN_APERTURE_CORRELATION)",
 "    if ok.sum() < D.MIN_APERTURES:",
 "        return np.nan",
 "    beta, _, _ = D._robust_line(travel[ok], d[ok], np.maximum(c[ok], 0)**2)",
 "    return -beta[1]",
 "",
 "travel_ch = t0_frozen + p_frozen * coord      # per-channel travel time",
 "levels = np.array([-1e-2, -5e-3, -1e-3, 0.0, 1e-3, 5e-3, 1e-2])",
 "",
 "demo = []",
 "for eps_in in levels:",
 "    shifted = D._shift_gather(aligned[b], q_ax, -eps_in * travel_ch)",
 "    demo.append(recover(shifted, references[b], layout, travel, q_ax, cfg))",
 "demo = np.array(demo)",
 "",
 "for a, r in zip(levels, demo):",
 "    print(f'  injected {a:+.1e}   recovered {r:+.3e}' +",
 "          ('' if a == 0 else f'   ratio {r/a:.3f}'))",
)

code(
 "fig, ax = plt.subplots(figsize=(6.5, 6.2), constrained_layout=True)",
 "lim = 1.2e-2",
 "ax.plot([-lim, lim], [-lim, lim], ls='--', color='0.6', label='perfect recovery')",
 "ax.plot(levels, demo, 'o', ms=12, color='#eb6834', mec='white', mew=1.8,",
 "        label='this burst')",
 "ax.set(xlabel='injected fractional change', ylabel='recovered fractional change',",
 "       title='Injection–recovery, single burst', xlim=(-lim, lim), ylim=(-lim, lim))",
 "ax.legend(frameon=False); ax.grid(alpha=.2); ax.set_aspect('equal')",
 "plt.show()",
)

md("### The full calibration",
   "",
   "One burst shows the machinery works. The sensitivity number needs all 23",
   "held-out bursts × 15 injection levels × 2 legs, plus timing controls — 1380",
   "trials, about 16 minutes at 96 GB. That is the production run:",
   "",
   "```bash",
   "sbatch awd_clean/deep_dvv_injection_recovery.sbatch     # --stage all",
   "```",
   "",
   "Its results are loaded below rather than recomputed.")

code(
 "rows = [r for r in csv.DictReader(open(ANALYSIS/'deep_dvv_summary.csv'))",
 "        if r['population']=='heldout' and r['band']=='15_30' and r['pass']=='primary']",
 "",
 "fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)",
 "for legname, colour in (('outbound', '#eb6834'), ('return', '#1baf7a')):",
 "    sub = sorted([r for r in rows if r['leg']==legname],",
 "                 key=lambda r: float(r['injected_dvv']))",
 "    x = np.array([float(r['injected_dvv']) for r in sub])*100",
 "    y = np.array([float(r['median_estimated_dvv']) for r in sub])*100",
 "    e = np.array([float(r['robust_scatter_1p4826_mad']) for r in sub])*100",
 "    ax.errorbar(x, y, yerr=e, fmt='o', ms=8, capsize=0, elinewidth=2,",
 "                color=colour, mec='white', mew=1.4, label=f'Deep {legname}')",
 "ax.plot([-1.2, 1.2], [-1.2, 1.2], ls='--', color='0.6')",
 "ax.set(xlabel='injected change (%)', ylabel='recovered change (%)',",
 "       title='Full calibration, 23 held-out bursts', xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))",
 "ax.legend(frameon=False); ax.grid(alpha=.2); ax.set_aspect('equal')",
 "plt.show()",
)

md("### The headline numbers",
   "",
   "**Noise floor** = the 95th percentile of what comes back when *nothing* was",
   "injected. **Reliable detection** = the smallest injected level where at least",
   "95% of trials both clear that floor and get the sign right, in both directions.",
   "The same definition is applied to Nano, so the two are comparable.")

code(
 "comp = list(csv.DictReader(open(ANALYSIS/'deep_dvv_nano_comparison.csv')))",
 "print(f\"{'observable':30s} {'leg':10s} {'noise floor':>12s} {'reliable':>10s}\")",
 "for r in comp:",
 "    if r['population'] not in ('heldout', 'all 46 bursts') or r['band_hz'] == '60-120':",
 "        continue",
 "    print(f\"{r['observable'][:29]:30s} {r['leg']:10s} \"",
 "          f\"{float(r['null_threshold_dvv'])*100:11.3f}% \"",
 "          f\"{float(r['reliable_level_detection'])*100:9.1f}%\")",
 "",
 "print()",
 "print('Note: the CSV also carries a lever_arm_s column. That is the CHANNEL SPAN,')",
 "print('not the regression lever arm — the regression runs over aperture centres,')",
 "print('which is a shorter span. The next cell computes the operative value.')",
)

md("### Why the lever arm does not convert one-for-one",
   "",
   "Deep's regression spans 11.8× more travel time than Nano's, which alone would",
   "predict roughly an order of magnitude better sensitivity. The realised gain is",
   "about 2×. The arithmetic below shows where it goes.")

code(
 "nano_scatter, nano_lever = 1.862e-3, 0.1210",
 "deep_scatter = float(next(r['robust_scatter_1p4826_mad'] for r in rows",
 "                          if r['leg']=='outbound' and float(r['injected_dvv'])==0))",
 "deep_lever   = travel.max() - travel.min()",
 "",
 "print(f'{\"\":26s} {\"Nano\":>12s} {\"Deep outbound\":>15s}')",
 "print(f'{\"lever arm (s)\":26s} {nano_lever:12.3f} {deep_lever:15.3f}')",
 "print(f'{\"null scatter\":26s} {nano_scatter:12.3e} {deep_scatter:15.3e}')",
 "print(f'{\"implied timing (ms)\":26s} {nano_scatter*nano_lever*1e3:12.3f} '",
 "      f'{deep_scatter*deep_lever*1e3:15.3f}')",
 "print()",
 "print(f'lever-arm ratio       {deep_lever/nano_lever:5.2f}x  (geometry, favours Deep)')",
 "print(f'timing ratio          {(deep_scatter*deep_lever)/(nano_scatter*nano_lever):5.2f}x  (repeatability, favours Nano)')",
 "print(f'net scatter ratio     {nano_scatter/deep_scatter:5.2f}x  = the two divided')",
)

# --------------------------------------------------------------------------
md("---", "", "## 9. What was checked, and what this does not show",
   "",
   "### Controls that passed",
   "",
   "| Control | Result |",
   "|---|---|",
   "| Synthetic injection, two independent constructions | agree to 5 significant figures |",
   "| Constant 5 ms timing shift | 4.87 ms into the intercept, 7×10⁻⁷ into ε |",
   "| Channel-order permutation, 499 scrambles | zero reach the observed beam power, all four leg/band tests |",
   "| Leave-one-aperture-out, leave-one-burst-out | nothing moves the result by as much as the noise floor |",
   "",
   "### A control that failed, reported as such",
   "",
   "The perturbed-trajectory test recovers injections even when pointed at the",
   "wrong arrival. The injection shifts *whole channel traces*, so any coherent",
   "energy in the extraction window carries the same imposed gradient and stays",
   "recoverable. It fires identically on both legs and therefore cannot",
   "discriminate. The return branch's weaker sensitivity is consequently reported",
   "as **unclassified** rather than explained after the fact.",
   "",
   "### Interpretation limits",
   "",
   "- The measured quantity is a fractional change in the **apparent along-fiber",
   "  speed of a guided mode**. Not formation V_P or V_S, not stress, pore",
   "  pressure, permeability, or strain.",
   "- **No depth resolution.** Position is distance along fiber; the depth mapping",
   "  is provisional and no result depends on it.",
   "- Attribution of the sensitivity to the guided mode specifically rests on the",
   "  permutation nulls, not on the injection experiment.",
   "",
   "### Where everything lives",
   "",
   "| | |",
   "|---|---|",
   "| Claim/evidence table, authoritative on numbers | `../DEEP_DVV_STATUS.md` |",
   "| Frozen analysis design, dated before results | `../DEEP_DVV_PREREGISTRATION.md` |",
   "| Pipeline order, script by script | `REPRODUCE.md` |",
   "| Assembled manuscript | `MANUSCRIPT.md` |",
   "| Preprocessing provenance | `../PREPROCESSING.md` |",
)

# ---------------------------------------------------------------------------
notebook = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "das", "language": "python", "name": "das"},
        "language_info": {"name": "python", "version": "3.9"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}
OUT.write_text(json.dumps(notebook, indent=1))
n_code = sum(1 for c in cells if c["cell_type"] == "code")
print(f"wrote {OUT.name}: {len(cells)} cells ({n_code} code, {len(cells)-n_code} markdown)")
