# C2 Phase 0 — results

Run 2026-08-13, jobs `38918261`, `38918803`, `38919086`. Scripts
`notebooks/c2_phase0_significance.py`, `c2_phase0_power.py`,
`c2_phase0_diagnose.py`. Reads `epoch_stacks_paired_deep_all.npz` only; the gate's
measurement is reproduced unchanged, so every difference below is statistical.

## Verdict

**The seven candidates are retired. The measurement that produced them could not
have detected a permeable fracture of any size, so C2 was never tested.** But the
reason is fixable and quantified, so C2 is *untestable as run*, not *refuted*.

I expected Phase 0 to end the thread. It ends the candidate list, and then finds
5× of usable headroom underneath it.

## 1. The gate is reproduced exactly

Same trajectory (wireline 1440 m/s at 252 ms, cemented 1330 m/s at 140 ms), same
σ = 1.20, same seven depths: 129, 131, 137, 696, 1280, 1295, 1609 m. Nothing
below is a reprocessing difference.

## 2. Why the candidate list means nothing

| Check | Wireline | Cemented |
|---|---|---|
| Channels in the fit | 1603 | 630 |
| Expected below −2σ if Gaussian | 36 | 14 |
| **Observed** | **7** | **0** |
| Binomial P(≤ observed) | 2.3×10⁻⁹ | 5.1×10⁻⁷ |
| Skew / excess kurtosis | +1.92 / +7.13 | +0.65 / −0.43 |
| Variance at wavelengths > 200 m | **71 %** | 11 % |
| Candidates at robust σ instead of std | **42** | 0 |
| Distinct features at the 51 m correlation length | **4**, not 7 | — |

Four separate ways of saying the same thing:

- **The count is a deficit, not an excess.** Seven where chance gives 36. A
  detection criterion satisfied by *fewer* events than chance is not detecting
  anything.
- **σ is not a scatter.** 71 % of the wireline residual variance sits at
  wavelengths longer than 200 m, so σ measures unmodelled long-wavelength
  structure that the linear log-amplitude trend fails to absorb. A threshold in
  units of that σ is not a threshold on channel-to-channel scatter.
- **The list is an artefact of the σ estimator.** Swapping the standard deviation
  for a robust MAD — a choice with no scientific content — turns 7 candidates
  into 42. Strongly non-Gaussian residuals (kurtosis +7.1) inflate the std, so
  −2σ lands far out in a light lower tail.
- **Seven channels are four features.** The measured residual correlation length
  is 25 channels (51 m), five times the gauge length, so 129/131/137 is one
  feature and 1280/1295 another.

## 3. Nothing is significant against a proper null

Phase-randomised surrogates preserve the residual autocorrelation and destroy
localisation — the correct null for "is this dip special?", which a Gaussian
assumption is not.

| Statistic | Wireline p |
|---|---|
| Count below −2σ | 0.368 |
| Deepest residual | 0.733 |
| Largest sustained step | 0.289 |

## 4. The candidates are not the shape the physics predicts

A permeable fracture makes amplitude *stay lower below it* — a step. The gate
detected excursions. Measured sustained offset across each candidate, positive
meaning "stays lower below":

| Depth | Offset | |
|---|---|---|
| 129, 131, 137 m | — | at the range edge |
| 696 m | +0.044 ± 0.233 (+0.2σ) | right sign, not significant |
| 1280 m | −0.304 ± 0.579 (−0.5σ) | **wrong sign** |
| 1295 m | −0.273 ± 0.493 (−0.6σ) | **wrong sign** |
| 1609 m | +0.142 ± 0.282 (+0.5σ) | right sign, not significant |

**2 of 7 even point the right way**, none above 0.6σ.

## 5. The cemented "step" was the top of the hole

Part 1's automatic verdict fired on a cemented step at p = 0.007. It is an edge
artefact: the breakpoint stays pinned to whichever channel is the first allowed
one as the guard is widened (161 m → 168 m → 209 m), and dropping the shallowest
300 m takes it to p = 0.648. Log-amplitude is not linear in depth near the
source, where geometric spreading dominates. Retired.

## 6. The measurement had no power — this is the real finding

Planting synthetic features of known size into surrogates carrying the measured
autocorrelation:

| | Wireline | Cemented |
|---|---|---|
| Detection of a planted **dip**, as measured | max **34 %** at a 95 % amplitude loss | max 13 % |
| Detection of a planted **step**, as measured | max **41 %** at a 95 % amplitude loss | max 56 % |

**Neither fiber reaches 95 % detection anywhere on the grid — not even for a
95 % amplitude loss (26 dB).** The measurement could not have found a permeable
fracture if one were there. "No significant candidates" was never going to mean
anything, and C2 PASS and C2 FAIL would have been equally uninformative.

## 7. But 96 % of the scatter is static, and static is removable

Splitting the residual with the two burst halves — `o = s + n_o`, `e = s + n_e`,
so `(o−e)/2` isolates noise and its variance equals the full stack's noise
variance:

| | Wireline | Cemented |
|---|---|---|
| Half-to-half correlation ρ | **+0.922** | +0.326 |
| Total σ | 1.202 (3.33× in amplitude) | 1.313 (3.72×) |
| **Static** per-channel response σ | **1.177** (3.25×) | 1.050 (2.86×) |
| **Random** measurement noise σ | **0.241** (1.27×) | 0.789 (2.20×) |
| Static fraction of variance | **96 %** | 64 % |
| Headroom | **5.0×** | 1.7× |

The factor-3.3 channel-to-channel scatter that destroys the test is *reproducible
between independent halves of the bursts*. Reproducible scatter is not noise — it
is static per-channel response: coupling, sensitivity, gauge-length position.
That can be calibrated out. Random noise cannot.

Note ρ = 0.922 is not a quality metric here. It is high because the static term
dominates. What matters is that the wireline's random noise, σ = 0.241, is 3×
lower than the cemented fiber's.

### What calibration would buy

With the static response perfectly removed:

| | Wireline | Cemented |
|---|---|---|
| **Step** at 95 % detection | **18 % amplitude loss (1.7 dB)** | 78 % (13.0 dB) |
| **Dip** at 95 % detection | 78 % (13.0 dB) | never reached |

An 18 % step threshold on the wireline would be a genuinely useful
permeable-fracture test. This is an upper bound on what is achievable — no real
calibration is perfect — but it is the number that makes Phase 1 worth doing.

The step is the sensitive statistic and the dip is not, which is the right way
round: the step is what the physics predicts.

## 8. One loose end, deliberately not pursued

On the wireline, dropping the shallowest 600 m gives a step of 1.609 at **1033 m,
p = 0.000**, and the location is stable against the 300 m cut (1035 m) even
though the significance is not (p = 0.162). Four analysis windows were tried and
one gave a small p, after seeing the data.

That is post-hoc window selection and this project does not count it. It is
recorded as the one thing in C2 that is not obviously dead, and it would need a
pre-registered window and a held-out split to become anything. It is **not** a
result.

## What follows

The Phase 1–4 plan in `C2_PERMEABILITY_FOLLOWUP.md` stands, with one item now
first and mandatory: **estimate and divide out the static per-channel amplitude
response.** Until that is done, no threshold on this profile can work, and the
existing Phase 1 item "separate instrument from medium" is not one option among
four — it is the whole problem.

Cheapest route: a per-channel gain from an independent window — a different
frequency band, a pre-arrival time window, or ambient RMS — so the calibration is
not derived from the same samples being tested. Then re-run this script and check
that σ falls from 1.20 toward 0.24. If it does not, the static term is not a
simple gain and C2 should be retired for good.

**Do not cite the seven depths anywhere.** They are a threshold artefact.
