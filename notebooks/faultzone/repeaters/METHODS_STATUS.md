# Repeaters project — status and methods spec

Written 2026-08-04. Durable record of the literature methods review and where the
project stands, so none of it lives only in terminal scrollback.

---

## 1. Where the project actually stands

| Claim | Status | Evidence |
|---|---|---|
| Fiber resolves earthquake coda | **Established** | G1: 785/900 channels, coda SNR > 3 |
| Repeating earthquakes recorded during DAS period | **Established** | HRSN control: 7 pairs, CC 0.917–0.995, baselines 39–440 d |
| DAS detects those repeaters | **Established** | DAS CC 0.63–0.79 vs DAS null max 0.328 |
| DAS CC systematically below HRSN | **Established** | 0.15–0.30 deficit, identical pairs |
| Cause of the deficit | **Confirmed** — flat stacking. Aligned CC 0.956 vs HRSN 0.954 | `moveout_test.csv` (§17) |
| The +0.124 null bias has the same cause | **REFUTED.** Alignment raises the null to +0.376 | §2.5, §17 |
| Moveout correction improves *detection* | **REFUTED.** d′ falls 2.18 → 1.74 | §17 |
| dv/v at HRSN depth (251–284 m) | **NULL, quantified** | G3: 0/10 pairs above a 1.94% control floor (§12) |
| dv/v in the shallow DAS interval | **Estimator not yet trustworthy** | G5 first run: common-mode systematic (§13) |
| Depth of any velocity change | Open — the remaining question | bounded above 75 m only once §13 is fixed |

### The seven confirmed repeater pairs

| DAS CC | HRSN CC | stations | separation | events |
|---|---|---|---|---|
| 0.792 | 0.987 | 8 | 39 d | 2025-04-02 / 2025-05-11 |
| 0.790 | 0.995 | 5 | 272 d | 2024-07-08 / 2025-04-06 |
| 0.775 | 0.993 | 3 | 305 d | 2024-06-23 / 2025-04-23 |
| 0.742 | 0.962 | 8 | 39 d | 2025-04-02 / 2025-05-11 |
| 0.683 | 0.991 | 5 | 440 d | 2024-05-13 / 2025-07-27 |
| 0.677 | 0.917 | 5 | 331 d | 2024-05-10 / 2025-04-06 |
| 0.631 | 0.945 | 4 | 59 d | 2024-05-10 / 2024-07-08 |

2024-05-10, 2024-07-08 and 2025-04-06 all correlate mutually — a **three-event
sequence spanning 331 days**. Nulls: DAS median 0.124 / max 0.328; HRSN median
−0.004 / max 0.347.

These supply the CWI baselines (272, 305, 331, 440 d) that Week 2 requires.

---

## 2. Diagnosis of the DAS/HRSN deficit

### 2.1 It is additive noise, not lost similarity

Converting CC to in-band SNR via ρ = SNR/(1+SNR):

- ratio α = ρ_DAS/ρ_HRSN = **0.749 ± 0.053** across all seven pairs (CV 7%)
- a genuine decorrelation mechanism would make α vary; it does not
- **one constant instrument attenuation explains every pair** → fixable in processing

Deficit range 7.2–17.2 dB, mean ~13 dB.

### 2.2 It is not fiber ageing or coupling drift

r(CC_DAS, days) = **−0.11** over a 39–440 day range. Ruled out.
r(CC_DAS, CC_HRSN) = +0.68 — DAS tracks real similarity, just compressed.
DAS SNR spans 1.71–3.81 (2.2×); HRSN spans 11.1–199 (18×). DAS is pinned at SNR ≈ 3.

### 2.3 The organizing insight I had been missing

**Every deterministic DAS effect cancels in a repeater pair**, because it is
identical for both occurrences: gauge length, cos²θ directional response,
strain-rate differentiation, coupling, cable packaging. These explain why DAS
differs from a *seismometer*. None can explain why two DAS recordings of the same
source differ from each other. Only the *random* and *event-dependent* parts
survive. This demotes four of my six candidate causes.

### 2.4 Cause #1 — stacking across uncorrected moveout. ~27 dB. Essentially all of it.

**My code had the physics backwards.** `correlate_all.py:61` asserts the arrival is
flat across channels because incidence is near-vertical. For a **vertical** fiber,
near-vertical incidence gives the *steepest* moveout the geometry allows —
apparent velocity equals formation velocity, ~919/3000 = 0.31 s end to end.
Broadside incidence is the flat case.

Stacking N traces smeared over delay T = convolution with a boxcar of width T:
signal survives as |sinc(πfT)|, noise falls as 1/√N.

| V_app | T | sinc at 12.5 Hz | 1/√700 | net gain vs ONE channel |
|---|---|---|---|---|
| 2000 m/s | 0.460 s | 0.040 | 0.0378 | 1.06 (0.5 dB) |
| 3000 m/s | 0.306 s | 0.043 | 0.0378 | 1.12 (1.0 dB) |
| 4000 m/s | 0.230 s | 0.043 | 0.0378 | 1.15 (1.2 dB) |

**A 700-channel stack is buying ~1 dB over a single channel.** Correcting the
moveout gives √700 = 28.5 dB, i.e. **+27 dB recovery** — the whole deficit.

### 2.5 The same bug explains the +0.124 null bias — **REFUTED, see §17**

~~For independent signals CC must be zero-mean. HRSN's null sits at −0.004; DAS's at
**+0.124** under identical processing. Boxcar smearing of width T imposes nearly
the same narrowband ringing on every event, so the stack's own impulse response
dominates the waveform shape and every pair correlates positively. **Two symptoms,
one cause** — the consistency is what makes the diagnosis credible.~~

**This was wrong.** The prediction below was written to be falsifiable and it
failed: correcting the moveout drove the null the *opposite* way, +0.124 → +0.376.
The two symptoms do not share one cause. §2.4 stands on its own; the origin of the
null bias is **still undiagnosed**. Do not cite "two symptoms, one cause" as support
for the moveout diagnosis — that support is withdrawn.

~~Independent prediction: correcting moveout must pull the DAS null toward 0.000.~~

### 2.6 Cause #2 — genuinely noisier per channel. Real and published.

Lellouch, Lindsey, Ellsworth & Biondi 2020 (SRL, doi:10.1785/0220200149), on
collocated downhole DAS and geophones: *"Whereas single DAS channels are
substantially noisier than geophones at the same location, their large number and
spatial coherency allow for the application of effective array processing
techniques."* Magnitude completeness −1.4 (DAS) vs −1.7 (geophones) — a 0.3 unit
(~6 dB) penalty *after* array processing. Raw per-channel penalty ~15 dB.

Note the structure: DAS is noisier per channel, and **array processing is how you
get it back**. I was not doing array processing — I was averaging traces across an
unflattened wavefront, which is not the same thing.

### 2.7 Cause #3 — S-wave blindness. Second-order, but does NOT cancel.

For a vertical borehole, ε_axial = ∂u_z/∂z. A vertically propagating S wave has
horizontal particle motion, so ∂u_z/∂z → 0: **a vertical fiber is nearly blind to
vertically-incident S**, while seeing P well. HRSN is 3-component and captures S at
full amplitude. Corroborated by Wang et al. 2018 (GJI, doi:10.1093/gji/ggy102):
collocated DAS and geophones matched on the first several P cycles but S
comparisons were distinctly poorer.

Consequence: **weight windows toward P, not S.**

### 2.8 Ruled out — gauge length

|sinc(πfG/V)| with G = 16.34 m, worst case V_app = formation velocity:

| f | V=3000 | V=4000 | V=5000 |
|---|---|---|---|
| 10 Hz | 0.995 | 0.997 | 0.998 |
| **20 Hz** | **0.981** | **0.989** | **0.993** |
| 50 Hz | 0.882 | 0.933 | 0.957 |

−1 dB point is 48–80 Hz, above the 50 Hz Nyquist. **1.9% attenuation at 20 Hz.**
Cannot account for any part of 13 dB. And it would cancel in the pair anyway.

---

## 3. Independent confirmation from Lellouch's SAFOD papers

On **this exact fiber**, which settles the flat-vs-dipping contradiction in my code:

- velocities estimated "from near-vertically propagating earthquakes using **slant
  stacks**" → the moveout is real and *is* the signal
- they abandoned channel-by-channel picking because "individual DAS channels are
  noisier than clamped receivers," adopting local slant-stack decomposition →
  causes #1 and #2 in one sentence
- "a **loop failure** at the fiber end limited analysis to **800 m**" → my channel
  cut at 800 was right, for the reason I guessed
- **`github.com/ariellellouch/SAFODDAS/`** hosts his SAFOD earthquake records *and*
  cross-correlated waveforms publicly — a direct benchmark for my chain
- setup: OptaSense ODH3.1, fiber to 864 m, ~1 N tension, 0.9 mm steel tube cemented
  between casing strings

---

## 4. Thresholds — my approach is the convention, not a workaround

Nobody imports an absolute DAS CC threshold. Ellsworth's own group builds the null:

| Study | Combination | Null | Threshold |
|---|---|---|---|
| Lellouch et al. 2021 (JGR, doi:10.1029/2020JB020462) | per-channel CC, SNR-weighted correlogram stack | MAD ≈ 0.03 from ~50,000 pairings | cluster 6×MAD ≈ 0.20, detect 0.09 |
| Li & Zhan 2018 (GJI, doi:10.1093/gji/ggy359) | per-channel CC, network-mean stack | MAD | 9×MAD |

Lellouch et al. validated with a **time-reversed acausal template**, which can only
produce random CC — no acausal value exceeded 0.09.

**Two fixes needed in mine:**

1. **My null is contaminated** (+0.124, should be 0.000). Any threshold derived from
   it is biased. Re-measure after fixing moveout and confirm it returns to zero.
2. **Add the acausal test** — correlate time-reversed A against B for all pairs. More
   trustworthy than a random-pair null, since random pairs can share genuine
   structure (same path, same noise field) while a time-reversed waveform cannot.

Given a clean null, 0.63–0.79 against max 0.328 is roughly **8–15 MAD** — well
above the 6×MAD convention. Report as MAD multiples *and* absolute CC, always with
the null stated.

---

## 5. Trace stack or correlogram stack? Two objectives, two answers.

This is not stated cleanly anywhere in the literature, so, explicitly:

- **Stack traces, then correlate** → maximizes the CC *value* (signal adds as N,
  noise as √N). Use when comparing against a published threshold. **Requires
  moveout correction** — that is the 27 dB.
- **Correlate per channel, then stack correlograms** → maximizes *detectability*,
  not the CC value. The stacked peak converges on the *mean per-channel* ρ, which is
  low; what improves as √N is peak-to-background. This is why DAS detection papers
  quote thresholds as low as 0.09. **This is what all the DAS papers do**, because
  their objective is detecting unknown events.
- **Count channels exceeding a threshold** → the DAS analog of "≥5 stations"
  (Waldhauser & Schaff 2008), robust to a few badly-coupled channels. Report
  alongside, don't substitute.

My objective is similarity of *known* pairs, so the moveout-corrected trace stack is
the better estimator. **Compute all three and cross-check** — they fail differently,
and disagreement is itself diagnostic.

Do **not** use cross-channel coherence for identification: Schaff & Waldhauser 2005
warn it recovers only a fraction of the window length and biases toward zero. Use it
for sub-sample delay after alignment only.

---

## 6. Standard repeater windows and bands

| Source | Window | Band | Criterion |
|---|---|---|---|
| Schaff & Waldhauser 2005, doi:10.1785/0120040221 | 1 s and 2 s on P and S; S travel time = 1.732 × P; lags ±1 s | 1.5–15 Hz | CC ≥ 0.7 at ≥4 stations |
| Waldhauser & Schaff 2008, doi:10.1029/2007JB005479 | 1 s P window (2 s as check; discard if delays differ >1 sample) | 1.5–15 Hz | mean C_f ≥ 0.9 at ≥5 stations |
| **Waldhauser & Ellsworth 2002**, doi:10.1029/2000JB000084 | **tapered 2.56 s (256 samples)** on the P or S train | **2–10 Hz** | **half of squared-coherency > 0.9** |
| Lengliné & Marsan 2009, doi:10.1029/2008JB006118 | 256 samples, 100 before P; coherency on 1.28 s centred on P | — | mean coherency, ≥5 stations |
| Igarashi 2020, doi:10.1186/s40623-020-01205-2 | P onset → S + 3 s (cap 50 s) | 1–4 / 2–8 / 4–16 Hz | CC ≥ 0.95 at ≥2 stations |

**Phase windows on picks — not the full waveform, not the coda.** My 11 s window
from catalog origin is far too long.

**Ellsworth's own physical criteria** (Waldhauser & Ellsworth 2002): source-area
overlap > 50% at Δσ = 3 MPa; ΔM < 0.5 (later ±0.3); drop the smaller event of pairs
< 1 month apart. Method lineage: Poupinet, Ellsworth & Fréchet 1984
(doi:10.1029/JB089iB07p05719).

**Magnitude/band coupling.** Brune corner frequencies at Δσ = 3 MPa, Vs = 3.5 km/s:
M0 → 228 Hz, M1 → 72 Hz, M2 → 23 Hz, M3 → 7.2 Hz. My 5–20 Hz sits below f_c for
essentially every event, which makes it forgiving of magnitude differences but blind
to fine discrimination. Add 4–16 and 10–40 Hz once SNR is fixed.

---

## 7. Recommended processing chain

Step 1 is worth ~27 dB; everything else is worth a few tenths.

1. **Moveout-correct before stacking.** Either use the SAFOD velocity model from
   Lellouch et al. 2019 (doi:10.1029/2019JB017533 — same fiber, P and S published),
   or the model-light route: scan incidence angle, align on the *relative* predicted
   travel-time curve (no absolute timing needed), pick the angle maximizing
   coherency. **Predicted: ρ → 0.993–0.999, i.e. HRSN parity even with a 10× haircut
   for imperfect alignment.** If it does not move substantially above 0.79, the
   diagnosis is wrong.
2. **Re-measure the null**; confirm it returns to ~0.000.
3. **Add the acausal (time-reversed) test** as the false-positive floor.
4. **Compute all three combination statistics** (§5) and report together.
5. **Per-channel preprocessing** (Lellouch et al. 2021 step 1, verbatim: "median
   removal, band-pass filter, removal of noisy channels, and trace-by-trace L2
   normalization"): per-sample median across channels; **drop channels below ~800 m**;
   demean, detrend, 5% cosine taper; zero-phase Butterworth order 4; L2-normalize.
   **No whitening, no running-absolute-mean** — those are ambient-noise tools from
   the CC pipeline and must not cross over.
6. **Windows: 1.0 s and 2.0 s starting 0.2 s before the P pick**, plus Igarashi's
   P→S+3 s. **Weight toward P** (§2.7). Require 1 s and 2 s delays to agree within
   one sample.
7. **Lag search ±1.0 s as a detector, not a correlation function** — fix the window
   on A, slide over B's un-padded continuous data, recompute the normalization
   denominator at each lag. A zero-padded FFT correlation recovers lags only to half
   the window length and biases CC downward on partial overlap. *Check whether my
   current code does this — plausible second contributor.*
8. **Sub-sample delay** by cross-spectral phase slope (Poupinet et al. 1984). ~1 ms.
9. **Multi-band 2–8, 4–16, 10–40 Hz.** True repeaters survive all three.
10. **Physical criteria**: ΔM ≤ 0.3; separation < source radius at Δσ = 3 MPa; drop
    the smaller event of pairs < 1 month apart. Where relocation is unavailable,
    substitute median ΔS–P across channels — scale-free, no absolute timing needed,
    which suits DAS.

### The one step with no published precedent

**No paper does moveout-corrected trace-stacking then pair correlation for repeater
ID on DAS.** The alignment machinery is published (Lellouch et al. 2021) but for
*detection*. Combining it with pair correlation is an inference from the sinc
arithmetic. That is a defensible novelty, but it is the step that cannot be
supported by citation — hence validating against the HRSN pairs first. The
seven-pair benchmark is an unusually good ground-truth set for tuning.

---

## 8. Running now

`moveout_test.py`, job 37524673. Scans slowness **including p = 0**, so the data
decides flat-vs-dipping rather than my argument deciding. Four variants isolate
alignment from window choice:

- flat / origin window — must reproduce the published 0.63–0.79
- aligned / origin window — alignment alone
- flat / P window — window alone
- aligned / P window — both (literature-standard)

Prints one of three verdicts including **"NO recovery — the diagnosis is WRONG."**
Outputs `moveout_test.csv`, `moveout_test.png`.

---

## 9. What I need help accessing

Top three are all Stanford-authored; #1 and #2 are the actual recipe for the fix on
this fiber. Stanford library login should make them one click each.

1. **Lellouch, Yuan, Ellsworth & Biondi 2019**, BSSA 109(6):2491 —
   `10.1785/0120190176` — semblance aperture, moveout/velocity scan ranges, channel
   quality criteria on *this* array
2. **Lellouch, Yuan, Spica, Biondi & Ellsworth 2019**, JGR 124(7):6931 —
   `10.1029/2019JB017533` — measured P/S velocity profile along this fiber; turns an
   assumed 2–6 km/s range into the real moveout curve
3. **Lellouch, Lindsey, Ellsworth & Biondi 2020**, SRL 91(6):3256 —
   `10.1785/0220200149` — measured per-channel SNR deficit in dB and the array gain
   actually achieved; direct benchmark for the 13 dB and the +27 dB prediction

Lower priority: Ichinose et al. 2022 (`10.1029/2022JB025101`, DAS fidelity vs
frequency); Nadeau et al. 1994 BSSA 84:247 and Nadeau & Johnson 1998 BSSA 88:790
(original Parkfield HRSN thresholds and measured repeater stress drops);
Abercrombie 2014 (`10.1002/2014GL062079`); Uchida & Bürgmann 2019
(`10.1146/annurev-earth-053018-060119`, Suppl. Table 2 compiles criteria for every
global repeater study); Rubinstein & Beroza 2005 (`10.1029/2005GL023189`, their
borehole-vs-surface comparison is the closest published analog to my HRSN control).

Note: Scite is exhausted, and both top-priority DAS papers returned
`contentDenied` / purchase-only through it.

### Two questions for Ellsworth, which no paper answers

1. For a single ~900 m downhole DAS array with no azimuthal coverage, what
   substitutes for the "≥5 stations" requirement in your 2002/2008 criteria — is a
   channel-count threshold on the correlogram stack sufficient, or does repeater ID
   require ΔS–P as an independent geometric check?
2. For measuring pair similarity (as opposed to detecting events), is
   moveout-corrected trace stacking followed by correlation preferable to
   per-channel correlation with correlogram stacking? Your FORGE workflow uses the
   latter for detection; the sinc arithmetic suggests the former should recover
   ~27 dB for a fixed pair.

---

## 10. Errors in this project so far, for the record

1. **Median subtraction destroying the signal** — `Xb -= median(...)` then
   `.mean(axis=0)` gives mean−median ≈ 0. Produced a believable but meaningless null.
2. **Applying HRSN's 0.9 threshold to DAS** — hid seven genuine repeaters.
3. **Asserting the arrival is flat across channels** — backwards for a vertical
   fiber; this is what led to (1) being justified by a wrong argument.
4. **Zero-lag CC with no lag search** — halved every correlation (0.42 → 0.893).
5. **Wrong catalog** — standard locations (±0.2–0.5 km) cannot predict which events
   share a patch.
6. **Half-merged manifests** — merged coverage for selection but not extraction, so
   second members were silently skipped.
7. **Second-resolution cache tags** — two events in one second collapsed to one file,
   giving a spurious CC = 1.000.
8. **A self-check that checked nothing** — an "expected under null" column recomputed
   the observed count.
9. **`sep=r'\s+'`** forces pandas onto its Python parser; 25-min timeout on a 108 MB
   manifest.
10. **574k `os.path.exists` on Lustre** — violates this repo's CLAUDE.md guidance.
11. **Window far too long** — 11 s from catalog origin, vs the 1–2 s phase windows on
    picks that the literature uses.

---

## 11. G0 — channel-to-depth registration (2026-08-04)

**Scripts:** `channel_depth_registration.py` (206-event scan, ~22 min),
`g0_refine.py` (post-processing). Outputs `channel_depth_registration.png`,
`g0_refine.png`.

**Result: PARTIAL PASS. Wellhead ≈ channel 15–25; absolute depth uncertain ±25 m.
The shallow channels ARE usable, which is the answer the project needed.**

### Three determinations

| method | wellhead channel | quality |
|---|---|---|
| pre-event noise transition | **23** | clean — lead-in is **35× noisier** |
| velocity-structure match to 2005 VSP | **8–23** | clean — two features matched |
| fibre-length arithmetic | 50 | **outlier, unexplained** |
| travel-time match | 155 | uninformative (minimum spans ch 11–199) |

### The cross-validation that matters

| | DAS local slant stack | 2005 check shot |
|---|---|---|
| slow layer | **1493 m/s** at ch 62 | **1496 m/s** over 20–75 m |
| fast peak | 3731 m/s at ch 118 | 3439 m/s over 75–150 m |

Independent instruments 20 years apart, agreeing to **0.2%** on shallow-layer
velocity. This also shows the DAS resolves structure shallower than Lellouch et al.
attempted (they stop at 50–75 m and flag the shallow section unreliable for V_P/V_S).

### Shallow channel verdict — the gating question

Channels 0–23: 35× baseline noise → uncemented lead-in, unusable.
Channels 23–896: noise flat within ±10% → uniformly well coupled.
Semblance 0.37 at ch 30 rising to 0.80 by ch 86; the lower shallow values are
expected (slow layer ⇒ shorter wavelength ⇒ more moveout across a 61-ch window),
not a defect.

**`CH_LO = 100`, used in every script in this directory, is depth 79 m — it was
discarding the entire depth range where the seasonal signal is largest.**
Li & Ben-Zion's ~17 m peak sensitivity is **channel ≈40**, which is coherent and at
baseline noise. The fibre samples the peak of the target signal.

### Two errors in the first run, corrected in `g0_refine.py`

1. **Semblance is not a discriminator here.** Median semblance is 0.3–0.85 across
   *all* 900 channels and never approaches the incoherent floor of 1/61 = 0.016.
   The arrival is coherent along the whole array, so "where coherent moveout
   begins" has no answer; the threshold crossing at ch 30 was an arbitrary cut.
   Pre-event noise is the real discriminator.
2. **Matching velocities instead of travel times.** The reference was built by
   smoothing t(z) and differentiating, which amplified digitisation noise (the
   reference swung to 8574 m/s) and produced a broad monotonic plateau with no true
   minimum. Travel time is the integral and is smooth — but it turned out to be
   *too* smooth to constrain the offset at all. **Localised velocity features are
   what carry the information**; the fix was to match the slow layer and the fast
   peak individually, not the whole curve.
3. **Incomplete arithmetic.** The 54-channel prediction ignored channels past the
   fibre terminus. With a dead tail of D, lead-in = 900 − 846 − D.

### Unresolved

The ~27 m gap between the noise/velocity determination (~20) and the length
arithmetic (50). Not explained by well deviation (PGSI `WELL_DEP` vs `REC_DEP`
differ by only 0.9 m at 1250 m) nor by cable overstuff (0.1–0.3% ≈ 2 m over 864 m).
Per Madsen et al. 2016 this residual is the signature of a **fibre accumulation in
the wellhead area**, which breaks linearity worst in the shallow section. **The
acquisition-metadata ask stands** (`StartLocusIndex`, lead-in length, OTDR/tap test).

### Also corrected here

The near-surface velocity is **not ~300 m/s**. That figure came from dividing depth
by traveltime when the check-shot source sits **45.72 m** from the wellhead: 66.4 ms
at 20.3 m depth is along a **50.0 m slant path** → 754 m/s apparent. Vertical-
corrected interval velocities: **1496** (20–75 m), 3439 (75–150), 3361 (150–300),
3091 (300–500, a low-velocity zone), 3495 (500–800), 5111 m/s (800–1250). This
removes the "factor of 8" tension with Lellouch's regional model — there was none.
Consequence for sensitivity: 4% over the top 75 m gives **2.7 ms** (≈5σ against
0.52 ms repeatability); 1% gives 0.68 ms, marginal.

### Reference data located (was not being used)

`awd_clean/pgsi_reference/PGSIarray_rec_coords_pos1.txt` — the 2005 PGSI 80-level
array with **known** `WELL_DEP` and `REC_DEP`, 46.68 → 1250.64 m at 15.24 m (50 ft)
spacing, UTM NAD27. Two further deployments (pos2 930–2134 m, pos3 1539–2743 m).
Raw check-shot SEG2 files are in `pgsi_reference/Check shots/`. This is the depth
ground truth and is far better than the digitised PNG.

---

## 12. G3 — dv/v on HRSN: NULL at 251–284 m (2026-08-04, job 37585615)

`dvv_hrsn.py`, coda stretching on the seven confirmed pairs, three overlapping
lapse windows (1–6, 2–12, 4–14 s), bootstrapped over stations.

**0/10 pairs exceed max(2·err, control floor). Control floor 1.944%** (90th
percentile of |dv/v| on 40 random non-repeating pairs). Repeater |dv/v| medians
are 0.024–0.035% against a control median of 0.78–1.02% — the repeaters are an
order of magnitude *quieter* than the null, which is what a real null looks like
when the estimator is working.

Supporting controls:

- coda halves: median |half1 − half2| 0.062% vs median |dv/v| 0.029%;
  r(half1, half2) = −0.14. No coherent change in either half.
- r(|dv/v|, seasonal phase) = +0.16, r(|dv/v|, elapsed) = +0.18. Neither.
- selection bias: r(selection CC, |dv/v|) = −0.46. Mildly negative, so the pair
  set is if anything censored *against* signal — the null is conservative.

**Why this is a result and not a failure.** Li & Ben-Zion 2023 put the seasonal
signal in the top few tens of metres, peak sensitivity ~17 m. HRSN sits *below*
that. A null at 251–284 m localises any change to shallower depth and makes the
shallow DAS channels — which G0 showed are usable from channel 23 — the whole
point rather than a bonus.

Note the asymmetry this creates: HRSN is the better instrument everywhere except
at the one depth that matters, where it has no sensors at all.

---

## 13. G5 — shallow differential timing: first run, and why its numbers are not
usable yet (job 37588342)

`g5_shallow_dvv.py`. Per-channel differential delay of the direct arrival between
the two occurrences of a pair, then dv/v as the slope of delay against one-way
travel time (`dvv_core.slope_dvv`), interval by interval down the fibre. The
intercept absorbs the origin-time difference, so no absolute timing is needed.

**The deep half of the measurement works.** 250–550 m gives −0.29 ± 0.34% and
550–850 m gives +0.28 ± 0.52%, both consistent with zero and consistent with G3.
That is the machinery validating itself on the depth range where an independent
instrument already says the answer is zero.

**The shallow half is a common-mode systematic, not a signal.** The fine profile
returned −26% (2–25 m), −21% (25–50 m), −3.7% (50–100 m), −2.4% (100–200 m),
turning positive below 200 m. Three things say this is processing, not medium:

1. **Pair-to-pair scatter far below the mean.** 25–50 m gave −21.4% with ±6.3%
   across ten pairs whose baselines run 39 to 440 days at different seasonal
   phases. A seasonal change *must* vary between those pairs. One that is
   identical for all of them is something being done to every pair alike.
2. **The two bands disagree.** 5–30 Hz gave −1.36 ± 2.37% over 2–75 m, 5–20 Hz
   gave −6.11 ± 1.97% on the same data. A medium change shifts the arrival by the
   same *time* in every passband. A cycle skip does not.
3. **The magnitude is not physical.** −21% over 25–50 m is not a seasonal
   velocity change in any published record.

Diagnosed cause: **residual cycle skipping**. The per-channel delay residual is
~2–3 ms (printed as `sigma ms` in the new output), while the rejection guard was a
fixed |dt| < 20 ms — a 7σ outlier that a least-squares slope over a 65 ms interval
cannot absorb. The shallow intervals are hit hardest because they have both the
fewest channels (~72) and the least travel-time leverage.

Two controls were also undersampled to the point of being misleading: the acausal
test returned "CONTAMINATED" off **two** surviving pairs, and the random-pair
floor came back NaN because dissimilar waveforms fail the coherence gate and yield
no measurement at all. The gate refusing to time incoherent pairs is the control
working, but it leaves the floor unsampled, and that needs saying rather than
reporting a NaN.

### Fixes in the rerun (job 37589998)

- `slope_dvv` now rejects iteratively on the **residual about the fitted line** at
  3.5 robust sigma, not on |dt|. Delays legitimately trend with depth, so a fixed
  |dt| cut either keeps the skips or discards real signal at the interval ends.
- Per-channel MAD gate at 4σ after the coarse 20 ms guard.
- Reported error is now **max(inverse-variance, pair scatter)**, so a claim never
  rests on the more flattering of two defensible numbers.
- **Common-mode test** printed per interval: flags any interval whose pair-to-pair
  spread is below 35% of its mean as systematic rather than seasonal.
- **Band-agreement check** promoted to a verdict-level veto.
- Controls require n ≥ 5 pairs before they are allowed to clear *or* condemn.
- `sigma ms` printed per interval — the per-channel timing precision is what
  limits every number in the table and it belongs in the output.

### The sensitivity arithmetic, for reference

err(dv/v) ≈ σ / (√N · spread of tau across the interval). With σ ≈ 2 ms:

| interval | tau | N ch | floor per pair | over 10 pairs |
|---|---|---|---|---|
| 2–75 m | 65.4 ms | ~72 | ~1.2% | ~0.4% |
| 250–550 m | 97.9 ms | ~300 | ~0.4% | ~0.13% |

So the shallow interval is intrinsically the *least* sensitive part of the fibre,
which is unwelcome but not fatal: even a clean 0.4% upper bound above 75 m,
combined with G3's 1.94% at HRSN depth, is a two-depth bound worth reporting.

### Still unaddressed, and probably real

The shallowest bins sit where two effects are strongest and neither is in the
current model: **free-surface interference** (above the reflection point the
fibre records upgoing and downgoing waves together, so "delay" is not a travel
time) and **gauge-length averaging** — 16.34 m in a 754 m/s layer at 17 Hz is
0.37 wavelength, far more aggressive than the 1.9% figure in §2.8, which assumed
V ≥ 3000 m/s. Both are depth-dependent and both could survive the fixes above.
Treat any residual shallow trend as suspect until they are modelled.

---

## 14. Deliverable B — catalog closeout (queued 2026-08-04)

**Duplicate tags: resolved, and it was not what §10 item 7 assumed.** There is
exactly one collision, 2024-12-30 04:39:33. The catalog shows nc75109596 (M0.98,
18 stations, gap 60, rms 0.07) and nc75109601 (M1.02, 16 stations, gap 84, rms
0.09) at the **same origin time to the millisecond**, 0.3 km apart. That is the
NCSN double-listing one earthquake with two solutions, not two earthquakes in one
second. The cache collapse was therefore physically correct and no re-extraction
is needed; the second listing just has to be dropped so it cannot pair with
itself at CC = 1.000. Handled in `correlate_perchannel.py`, with the reasoning in
the code so it is not re-litigated.

**Per-channel correlation (job 37588824).** `correlate_perchannel.py` had been
written but never run — the explicit loop is 206 events × 21,115 pairs × 700
channels × 3 bands, which does 206 transforms' worth of work 21,115 times. Now
precomputes one rFFT per event per band and forms each pair's stacked correlogram
as a single frequency-domain sum, since correlogram stacking is linear. About an
hour per band becomes about two minutes. The explicit loop is retained and the
fast path is **asserted against it on four pairs at startup**, so a transcription
error cannot quietly produce a plausible wrong threshold.

Emits `perchannel_candidates.csv`: the top 40 pairs ranked by their **weakest**
band, so a pair that survives 2–8, 4–16 and 10–40 Hz outranks one that is
spectacular at low frequency and gone by 10–40.

**HRSN extension (job 37589048, chained afterok).** `hrsn_extend.py` takes those
40 to HRSN. The existing confirmation covers only the top 10 from the *trace-
stacked* search, and that ranking is compromised by the moveout bug of §2.4 — the
comb filter that depresses every DAS CC by 0.15–0.30. Pairs it ranked 11th–40th
are not reliably worse; they are samples from a blurred ordering. Fetches are
cache-first and fail per event, so a node without an outbound route degrades to a
cache-only run instead of dying.

The interesting outcome is not only "more pairs". If the per-channel ranking
recovers the same seven and adds nothing, that is itself informative: the moveout
bug depressed the CC **values** but did not scramble the **ordering**.

Every additional confirmed pair adds a CWI baseline, and baselines are the binding
constraint on separating the seasonal term from secular drift — that separation
works only because elapsed time and seasonal phase are nearly uncorrelated across
the pair set, which improves with more pairs.


---

## 15. Correction — there IS a seismometer in the SAFOD hole (2026-08-05)

**§11 and the plan both state that no seismometer occupies the SAFOD borehole
during the DAS period. That is wrong.** The query behind it
(`station_geometry.py`) restricted to networks `BP`, `NC`, `PB`; the SAFOD
instruments are on network **`SF`**, which was never requested. A plain search
error, found only when the claim was challenged.

`safod_hole_check.py` repeats it with no network restriction against NCEDC and
IRIS.

### What is actually there

**`SF.MH029` — SAFOD Main Hole, 2555.1 m depth, active 2022-05-21 to present.**
Three components (GP1/GP2/GP3) at **1000 Hz**, orientations az 35.5/dip 41.7,
az 81.33/dip 31.95, az 338.14/dip 20.1 — a triaxial package in a deviated hole, so
it must be rotated into a geographic frame before use.

The history is long: `SF.MH001`–`MH029` since 2004 at 1844–2765 m, plus the
`SF.PH*` Pilot Hole series including a 32-level array from 856 to 2096 m.

**Near miss:** `SF.MHVSM` had a sensor at **782 m** — inside the fibre's range —
but ran only to 2023-12-31, five months before DAS recording began.

### What survives, and what does not

**Survives:** the fibre spans 0–864 m and MH029 sits at 2555 m, about 1.7 km below
the deepest channel. *Nothing else instruments 0–864 m.* The depth-range claim was
right; the blanket "no seismometer in the hole" was not.

**Does not survive:** the assertion that no same-hole reference exists for
amplitude calibration. One does, which makes Atterholt's suggestion about
collocated sensors far more sensible than the HRSN-based reading given earlier.

### Coverage — the number that decides its usefulness

`mh029_coverage.py`, over the 12 events in pairs above HRSN CC 0.90:

| | |
|---|---|
| events with MH029 data | **5 / 12** |
| **pairs with BOTH events** | **1 / 8** |

Every 2024 event has data (2024-05-10, 05-13, 05-22, 06-23, 07-08); no 2025 event
does (04-02, 04-06, 04-23, 05-11, 07-13, 07-27). Archiving stops between 2024-07
and 2025-04. Since every confirmed pair straddles that gap, only
2024-05-10 / 2024-07-08 is usable.

**Conclusion: MH029 cannot serve as a general confirmation channel for the
repeaters.** It is usable for amplitude calibration on the 5 events it recorded
(calibration needs single events, not pairs) and as one independent check on that
one pair. It is also **not continuously archived** — 1 of 6 random quiet windows
returned data — so it cannot support ambient-noise interferometry.

Worth asking NCEDC or the SAFOD data managers whether the 2025 gap is an archiving
lapse or an instrument failure; if the former, the data may be recoverable and the
coverage would change substantially.

### Attribution correction

An earlier message credited "collocated three-component sensors allow amplitude
calibration and provide local particle motion observations" to Atterholt's Garlock
paper. That came from a search summary, not from reading the paper, and the summary
did not distinguish between Atterholt, Zhan & Yang 2022
(doi:10.1029/2022JB025052) and Atterholt, Zhan, Yang & Zhu 2024. The attribution is
unverified.


---

## 16. The published-method pipeline, end to end (2026-08-05)

Replaces the ad-hoc sequence of scripts that produced the earlier "8 confirmed
pairs". Follows Nadeau's Parkfield procedure, with each threshold traced to source.

### Ellsworth's actual guidance (email, 2026-05-20/21)

DDRT catalog, NCEDC, 2024/05/01–2026/04/01, `delta=35.982,-120.544,0,15` → 329
events. *"Spots with multiple colors are the ones I would focus on. Larger events
would be better too."* And the framing sentence: *"you will want to find events
that are close in magnitude and location, **which will then need to be verified as
either repeaters or neighbors**."*

### The pipeline

| step | script | output |
|---|---|---|
| coverage audit, both manifests, file-interval resolution | `phaseA_coverage.py` | 208/329 events, 21,528 pairs |
| similarity β on HRSN, 3C, P and S | `beta_similarity.py` | `beta_similarity.csv` |
| clusters → sequences → creep | `sequences_and_creep.py` | `sequences_creep.csv` |
| presentation | `REPEATERS_dashboard.ipynb` | validated, all cells OK |

### Result

**2 sequences pass every published gate.**

| seq | events | interval | ΔM | radius | crack | Nadeau–Johnson |
|---|---|---|---|---|---|---|
| 0 | 2024-05-13 M0.78 → 2025-07-27 M0.84 | 440 d | 0.06 | 14.4 m | 0.87 mm/yr | **31.8 mm/yr** |
| 3 | 2024-07-08 M0.65 → 2025-04-06 M0.81 | 272 d | 0.16 | 13.2 m | 1.29 mm/yr | **49.1 mm/yr** |

N&J rates bracket Parkfield's geodetic creep (~25–30 mm/yr). The crack model is 25×
low — the known small-repeater underestimate that motivated N&J.

**The agreement is partly circular**: N&J was calibrated against geodetic creep at
Parkfield, so reproducing it is a consistency check, not an independent measurement.

### β does not reach Nadeau's 0.98, and why that is not a null

Max β = 0.9695; zero pairs at 0.98. But Nadeau found 63% of 1700 events in clusters,
so a fall to zero is not credible — the absolute scale differs by ~0.02 (band,
window, aggregation, 1987-era instruments). The *population* matches: their 294
clusters from 1700 events implies ~0.1% of pairs above threshold; we have 0.06%
above 0.90. β ≥ 0.90 is used as the scale-shifted equivalent.

**Caveat, stated because it is a judgement not a measurement:** our cluster count
declines monotonically (20, 20, 17, 16, 13, 10, 7, 5, 2, 0) with no plateau, so
Nadeau's own stability criterion does not select a value in our data.

### Three estimator errors found and fixed here

1. **Inflated null.** Maximising over ±1 s of lag inside a 1.5 s window maximises
   over ~30 lags of a ~45-DOF correlation; expected max ≈ 0.4, and the measured null
   median was 0.427. Fixed by one bulk alignment per station then ±0.1 s residual.
   Null median fell to 0.270.
2. **Window wraparound.** The bulk shift was applied with `np.roll` to an
   already-cut window; a 0.5 s shift is 125 samples against a 375-sample window, so
   a third of it folded in from the opposite end. Fixed by shifting the extraction
   indices.
3. **Bursts as recurrence.** Clusters with events 0–2 days apart gave creep rates up
   to 866,000 mm/yr. Waldhauser & Ellsworth's 1-month rule was in the plan and not
   implemented; now applied, collapsing bursts to one loading cycle.

### Retraction

An earlier note reported separations of 1–8 m for the confirmed pairs and concluded
source overlap was "decisively satisfied". That came from a broken
`horizontal_separation_km` column in a prior screen (its depth column was correct).
Real DDRT separations are **167–771 m**. Nadeau 1995 predicts exactly this: routine
locations scatter genuine repeaters over ~200 m, and only waveform-based relative
relocation tightens them to 10–20 m. **Collocation therefore cannot be tested with
DDRT** and is reported, never gated on.

### Settled along the way

- **Wellhead coordinate.** The surveyed PGSI collar is 216 m from the notebook value
  (35.974204, −120.552141 ≈ MH029) and 1243 m from Ellsworth's 35.982, −120.544.
  Use the notebook value; "MH030" is a different monument.
- **2017 data unnecessary.** 23 events near SAFOD in the 10-day 2017 window, but
  0 of 26 sequence locations had one within 500 m. The 8.4 TB on scratch adds
  nothing to this analysis. Disk 1 (with the tap test) had already purged.
- **The other student's catalog** is `ettore/research/projects/SAFOD/Catalogs/
  RepeatingEq.csv`, 10 events, Oct 2025. Four are absent from DDRT, so it used a
  different catalog; it shares **zero** events with this analysis.

---

## 17. Moveout correction: the deficit is confirmed, the null prediction is refuted

`moveout_test.py`, run 2026-08-04, results in `moveout_test.csv`. The output sat
unread for a day while four other directions were pursued and abandoned. Read the
results before building anything downstream.

Four variants isolate the two factors. "Aligned" applies one rigid per-channel time
shift per event, set by a slowness scan that **includes p = 0**, so the data chooses
between flat and dipping rather than the argument choosing.

| variant | rep median | HRSN | null median | null MAD | d′ | repMIN − nullMAX |
|---|---|---|---|---|---|---|
| published (flat / origin win) | 0.680 | 0.954 | +0.124 | 0.255 | **2.18** | **+0.299** |
| flat / origin | 0.679 | | −0.017 | 0.283 | 2.45 | +0.275 |
| **aligned / origin** | **0.944** | | +0.184 | 0.455 | 1.67 | +0.217 |
| flat / P window | 0.720 | | +0.289 | 0.165 | 2.61 | +0.221 |
| **aligned / P window** | **0.956** | | +0.376 | 0.333 | 1.74 | +0.177 |

### 17.1 Confirmed — §2.4 was right about coherent signal

Median repeater CC **0.680 → 0.956** against an HRSN median of 0.954. The deficit
vs HRSN goes 0.212 → **−0.006**: parity. All 10 pairs improve; the weakest goes
0.627 → 0.880. `correlate_all.py:69` was stacking 700 channels across an
uncorrected 0.31 s moveout, and that single line accounts for the whole
DAS/HRSN gap.

**`correlate_all.py` still contains the wrong physics in a confident comment**
(lines 61–64, "the arrival is FLAT across all 900 channels"). `correlate_perchannel.py`
and `g5_shallow_dvv.py` were fixed; this one was not. Every DAS CC number published
in this document above §17 came from the unfixed path.

### 17.2 Refuted — the null moves the wrong way

The written prediction was that the null must fall toward 0.000. It rose to
**+0.376**, tripling, with its spread growing too. §2.5 is withdrawn.

The cause appears to be physical rather than a coding error: alignment turns every
event's stack into a clean coherent P wavelet, and P wavelets from comparable
distances resemble one another. Coherence rises for unrelated pairs as well as for
repeaters. What discriminates a repeater lives in the scattered coda detail, which
an aligned P-window stack suppresses. Consistent with this, the *narrowest* null
(MAD 0.165) is flat / P, and the widest is aligned / origin.

### 17.3 Consequence — this does NOT buy detection

On every discrimination metric, correction makes things **worse**: d′ 2.18 → 1.74,
and the gap between the weakest repeater and the strongest control 0.299 → 0.177.

Therefore **moveout correction does not open template matching, catalog
enhancement, or detection below completeness.** That inference was drawn and is
retracted before use.

Where it does help is *waveform fidelity*, which is a different quantity from
discrimination and does not involve the null: dv/v, spectral fitting, differential
timing. Every such measurement in this project so far was made on a 0.68-coherence
stack when 0.956 was available — including the ±3.6 % dv/v figure recorded as a
negative result.

### 17.4 Limits of this test

- The null is **40 control pairs**, against 21,528 for the published null. The d′
  ordering is indicative, not established. Rebuild the null at full size before
  treating the detection conclusion as final.
- 0 of 40 controls exceed the weakest repeater in *any* variant, so all four
  configurations separate perfectly at this sample size; the d′ differences are
  margin, not error rate.
- Everything here is at **100 Hz** (`extract_all.py:57`, `desampling=True`), one
  fifth of the 500 Hz the files actually carry. Orthogonal to the moveout question
  and still untested.

---

## 18. Re-measurement at 500 Hz with aligned stacks: two reversals and one clean death

Everything in sections 12-13 and the stress-drop closure was measured at 100 Hz with
FLAT channel stacking, i.e. through a ~27 dB array-gain loss and a 5x bandwidth loss.
Three verdicts were re-run at the corrected configuration. Predictions were written
into each script's docstring before submission.

### 18.1 Stress drop — DEAD, and now for the correct reason

`stress_drop_500.py`, job 37781283. Registered gate S1 was: fc must fall with
magnitude, or the fit is reading the instrument rather than the source.

    r(M, log fc) = +0.246        WRONG SIGN -> VOID
    two M1.86 events give 0.47 and 0.01 MPa   (factor 47 at identical magnitude)
    only 3 of 8 events land in 0.1-100 MPa
    S2 passed: fc is depth-stable, median MAD/fc = 0.37

The earlier closure blamed magnitude ("M0.65's corner is above the band"). That was
the wrong reason: the spectral fit does not recover a source corner even for M3.17,
whose corner sits an order of magnitude below the band edge with SNR > 3 on 845
channels. **Do not reopen this on "use larger events" grounds again.**

### 18.2 dv/v — REOPENS. The +/-3.6 % negative result was an artefact.

`recheck_dvv_500.py`. Same estimator (`stretch_dvv_bootstrap`), same lapse windows
as G3; the only changes are moveout-aligned subarray stacks and the wider band.

| | G3, 100 Hz flat | 500 Hz, aligned |
|---|---|---|
| coda CC | 0.22 median (0.01-0.63) | **0.958** |
| precision | +/-3.6 % | **+/-0.033 %** |

Best configuration 5-20 Hz, lapse 2-12 s: repeater median |dv/v| 0.017 % with
bootstrap error 0.033 %, random-pair floor 0.99 %. Values are consistent with zero
at the 0.03 % level, i.e. a tight null rather than a detection.

The registered prediction was "improved but still insufficient, ~+/-0.6 %". The
measurement beat that by 20x, so the prediction was wrong in the conservative
direction. Section 12's "0/10 pairs above a 1.94 % control floor" is superseded.

**PROVISIONAL.** The bootstrap resamples only 4 subarrays. An error bar from 4
samples is not trustworthy, so +/-0.033 % is not yet a quotable number. The 110x
improvement in coda CC is solid; the precision figure needs many-subarray or
channel-level resampling before use.

### 18.3 Same-patch discriminant — significant, one control still missing

`interval_dvv_gate.py`, job 37772600. G2 synthetic recovery now PASSES at all three
injected levels (0.1/0.3/1.0 % recovered as 0.112/0.339/1.296 %), so the estimator
is unbiased. G3 acausal passes. G1 fails: null sigma 54 %, per-pair error ~2.6 %.

That per-pair error is **100x worse than the Cramer-Rao estimate** of 0.027 % quoted
in the plan, so the "1.9 m resolvable offset" figure is wrong; real per-pair
resolution is hundreds of metres. The measurement fails as a ruler.

It works as a discriminant:

| group | n | median inferred offset | IQR |
|---|---|---|---|
| family | 7 | 50.2 m | 547.7 m |
| control | 33 | 1551.1 m | 3164.5 m |

ratio 0.03, Mann-Whitney **p = 0.0008**.

**NOT YET A RESULT.** The controls were random, not matched. Family pairs were
selected for high HRSN correlation, and better-correlated pairs give
better-determined delays, hence less slope scatter and a smaller inferred offset —
with no geometry involved. Until controls are matched on distance, magnitude and
SNR, and until inferred offset is regressed against CC *within* the family group,
this separation may be entirely a similarity artefact. Do not cite p = 0.0008
without that control.
