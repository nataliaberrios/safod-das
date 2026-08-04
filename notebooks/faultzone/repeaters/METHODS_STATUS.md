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
| Cause of the deficit | **Diagnosed, under test** | job 37524673, `moveout_test.py` |
| dv/v measurable | Not started | Week 2; inputs now exist |

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

### 2.5 The same bug explains the +0.124 null bias

For independent signals CC must be zero-mean. HRSN's null sits at −0.004; DAS's at
**+0.124** under identical processing. Boxcar smearing of width T imposes nearly
the same narrowband ringing on every event, so the stack's own impulse response
dominates the waveform shape and every pair correlates positively. **Two symptoms,
one cause** — the consistency is what makes the diagnosis credible. Artifact-
corrected, (ρ−0.124)/(1−0.124), the repeaters are really 0.579–0.763.

Independent prediction: correcting moveout must pull the DAS null toward 0.000.

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
