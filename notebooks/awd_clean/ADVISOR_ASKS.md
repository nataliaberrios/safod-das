# Two things the advisor asked for — captured 2026-08-05, not yet scoped

Written down so they don't get lost. Both were relayed second-hand and neither
is fully specified yet; the open questions are listed so they can be settled in
one conversation rather than guessed at.

---

## 1. "Treat a channel as a virtual source and cross-correlate it with all the others"

**This already exists, twice.** Before building anything, check whether the ask
is actually for something new.

| Where | What it does |
|---|---|
| `stack_daily.py` → `cc_tools.computeCC` | 2024–25 ambient, one virtual-source channel against all others, 30 s windows at 50% overlap, 5–20 Hz, torch FFT. The legacy pipeline. |
| `sanity/sanity_cc.py` | Strict Lellouch et al. 2019 reproduction of the same thing, continuous files only, TN window 0.1 s. |
| `awd_clean/ambient_transfer_test.py` → `normalized_corr_pairs` | Batched FFT correlation, `max_lag` 0.35 s, feeding the whole signed F-K chain. |

So the likely readings of the ask, in order of plausibility:

1. **Do it on the June 2026 AWD active-source data**, not on ambient. That would
   be new — every existing virtual-source product is ambient noise.
2. **Do it as a plain picture** — one channel against all others, correlation
   versus channel and lag, no F-K filter, no wedge, no null. A gather rather than
   a statistic. This fits the "boring look" theme and is cheap.
3. **Vary which channel is the virtual source** and show the result is not an
   artifact of picking channel 0. All the existing work fixes the source at
   channel 0 (`ambient_fk_full_pipeline_null_v2.py` keeps channel 0 at
   coordinate zero even in its surrogates).

**RESOLVED 2026-08-05.** Reading 1 + 2 + 3 together, and the advisor also said
"deconvolution", which settles it. On a known active source you do not correlate
to *find* a Green's function -- you correlate to **remove the source**.
Correlating channel A with channel B cancels the source term common to both and
leaves the response between them, as though a source sat at A. That is the
virtual source method (Bakulin & Calvert 2006), and it redatums the surface
weight drop into the borehole.

Why it is the right move here: `docs/paper1/STATUS.md` already concludes the
source is the limiter (sigma_alpha ~ 0.30 ms common-mode, 39% amplitude CV).
Source cancellation attacks the measured bottleneck head-on, and would add an
entry to the 65x -> 5.7x -> 3.1x -> 1.6x sensitivity lineage.

Caveat: classic Bakulin-Calvert sums over a range of surface source positions to
satisfy stationary phase, and the AWD never moved. Deconvolution interferometry
needs no source aperture -- dividing B by A cancels the source spectrum whatever
its shape (Snieder & Safak 2006; applied to borehole arrays by Nakata & Snieder).
That is presumably why the advisor mentioned deconvolution.

Implemented in `awd_virtual_source.py`: correlation and deconvolution gathers on
the canonical stacks, four virtual-source positions (50/150/250/350 m) so the
result cannot be a channel-0 artifact, plain lag-vs-distance gathers with no F-K
filter, no wedge, no semblance and no null. Both operations also cancel any gain
applied equally to the whole record, so the read-time taper cannot bias them.

`LITERATURE.md` section 6 is missing Bakulin & Calvert 2006 and Snieder & Safak
2006; both should be added.

---

## 2. "Do something with the borehole seismometer and the DAS"

**What we have already:** the repeaters project uses HRSN borehole stations as an
independent instrument — `faultzone/repeaters/hrsn_control.py`,
`hrsn_extend.py`, `dvv_hrsn.py`, plus ~96 cached HRSN events and
`station_geometry.py`. `METHODS_STATUS.md` records the key result: DAS CC
0.63–0.79 against HRSN CC 0.917–0.995 on identical pairs, a systematic 0.15–0.30
deficit, cause under test.

So a DAS-vs-borehole-seismometer comparison **exists for the repeaters
workstream**. What does not exist is any use of a seismometer alongside the
**AWD active-source** data, which is probably what was meant.

Things that would be genuinely new, roughly in order of how well-posed they are:

- **Instrument-response / units cross-check.** Use a co-located borehole
  seismometer to check DAS strain-rate against velocity at the same depth. This
  is the one that would actually calibrate the amplitude story, which is
  currently uncalibrated everywhere (`AWD_advisor_figure_guide.tex` line 99:
  "amplitude before instrument and coupling calibration").
- **Independent confirmation of the candidate guided arrival.** A seismometer at
  depth either sees the ~1.4–1.5 km/s arrival or it doesn't.
- **Depth registration.** A seismometer at a known depth is a fixed point for the
  fiber-coordinate-to-depth mapping, which is the open problem on the Deep fiber.

**To ask:** *which* seismometer? The candidates are not the same instrument:
- the SAFOD main-hole / Pilot Hole geophone string,
- an HRSN borehole station (already in use, but at Parkfield distances, not in
  this hole),
- the Paulsson array under `$OAK/.../Paulsson_array_data/`.

Also: does it have coverage during the June 2026 AWD window, or only during the
2024–25 ambient period? That determines which workstream it can join.

---

## The Atterholt lead

`LITERATURE.md:120` already has **Atterholt et al. 2024, JGR 129** —
`04_methods_das/atterholt_2024_jgr_garlock-fault-zone-with-fiber.pdf`, Garlock
fault-zone structure from a DAS array. That is fault-zone imaging with fiber, and
is the most likely paper meant, but it is not obviously a DAS-plus-seismometer
methods paper. Worth reading the methods section to see which of the two asks it
actually speaks to before assuming.

Not verified: whether there is a different Atterholt paper specifically on
combining DAS with a borehole seismometer. Do not cite one without checking.
