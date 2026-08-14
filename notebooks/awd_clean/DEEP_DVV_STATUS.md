# Deep guided-mode sensitivity — status and claim table

Written 2026-08-06. Companion to the frozen
[`DEEP_DVV_PREREGISTRATION.md`](DEEP_DVV_PREREGISTRATION.md). Where numbers here
disagree with an older notebook cell or a figure caption, this file wins.

The question: what is the smallest fractional change in the apparent speed of the
Deep guided mode that can be recovered reliably, and does it beat the Nano
apparent-moveout observable?

---

## 1. Where this stands

| Claim | Status | Evidence |
|---|---|---|
| Deep outbound is more sensitive than Nano | **Established, restricted to the outbound branch** | reliable detection 5×10⁻³ vs Nano 1×10⁻²; null threshold 1.84×10⁻³ vs 3.25×10⁻³ |
| The Deep *installation as a whole* beats Nano | **Not claimable** | return leg reaches 1×10⁻², equal to Nano, not below it — pre-registered Rule 2 |
| The 11.8× lever arm buys 11.8× sensitivity | **Refuted** | it buys ~2×; per-burst timing repeatability is 7.6× worse and consumes the rest |
| Estimator separates source timing from propagation | **Established** | 5 ms constant shift → 4.87 ms into the intercept, 7×10⁻⁷ into ε |
| Recovery is unbiased on the primary branch | **Established** | bias ≤6×10⁻⁵ at every level, outbound |
| Result carried by one aperture or one burst | **Refuted** | no single aperture moves it >0.55× threshold, no burst >0.27× |
| Combining the two legs improves sensitivity | **Partial — improves the noise floor, not the detection level** | scatter 1.36–1.48× better; reliable detection unchanged at 5×10⁻³ |
| A more sensitive observable lowers the tidal upper limit | **Refuted** | UL 1.03×10⁻³ at 46 bursts, ~18× the Niu SAFOD expected-response scale; ~1.3× worse than Paper 1's Nano 7.98×10⁻⁴ |
| Return-leg weakness is SNR-limited | **Unresolved** | 5 of 6 pre-registered criteria met; the 6th rests on a control that turned out non-diagnostic |
| Sensitivity is attributable to the guided mode specifically | **Inherited, not shown here** | rests on the earlier permutation nulls (p = 0.002), not on this experiment |

## 2. Headline numbers

Held-out population (23 bursts, `epoch % 2 == 0`), 15–30 Hz, primary pass.
Nano row from `nano_dvv_injection_recovery.npz`, same estimator code.

| Observable | Null threshold | Reliable detection | Sign-only | Null scatter (MAD) |
|---|---|---|---|---|
| Nano apparent moveout, 46 bursts | 3.25×10⁻³ | 1×10⁻² | — | 1.86×10⁻³ |
| **Deep outbound** | **1.84×10⁻³** | **5×10⁻³** | 2×10⁻³ | 1.20×10⁻³ |
| Deep return | 3.05×10⁻³ | 1×10⁻² | 5×10⁻³ | 1.32×10⁻³ |
| Deep paired, equal weight | 1.82×10⁻³ | 5×10⁻³ | 2×10⁻³ | 8.84×10⁻⁴ |
| Deep paired, inverse-variance | 1.53×10⁻³ | 5×10⁻³ | 2×10⁻³ | 8.12×10⁻⁴ |
| Deep paired, covariance-aware | 1.67×10⁻³ | 5×10⁻³ | 2×10⁻³ | 7.87×10⁻⁴ |

"Reliable detection" is the Nano-comparable metric: smallest tested |ε| where
≥95% of trials both exceed the empirical null threshold and carry the correct
sign. Frozen trajectories: outbound 1544.6 m/s at t₀ = +0.100 s, return
1549.7 m/s at t₀ = +0.346 s — two independent halves agreeing to 0.3%.

## 3. The lever arm does not convert one-for-one

The central quantitative result of this analysis.

| | Nano | Deep outbound |
|---|---|---|
| Regression lever arm | 0.121 s | 1.428 s (11.8×) |
| Null scatter | 1.86×10⁻³ | 1.20×10⁻³ (1.55×) |
| Implied per-burst timing repeatability | 0.225 ms | 1.71 ms (7.6× **worse**) |

Deep's geometric advantage is real and large, but poorer burst-to-burst timing
repeatability consumes about four fifths of it. This is consistent with the
Deep-fiber per-shot repeatability deficit already recorded in
`docs/paper1/STATUS.md` (per-shot CC 0.054).

Do not quote the `lever_arm_s` field in the generated report as the lever arm:
that field is the channel span (1.808 s). The regression runs over aperture
centres, so the operative figure is **1.428 s**.

## 4. Paired-leg combination

Pre-registered as a comparison against outbound, never an automatic replacement.

**Between-leg null error correlation ρ = +0.121.** The two legs record the same
source impacts, so correlated errors were expected; they are almost absent. That
is why combining helps at all.

The gain is exactly what independent averaging predicts:

| | predicted (independent) | observed |
|---|---|---|
| equal weight | 8.92×10⁻⁴ | 8.84×10⁻⁴ |
| inverse-variance | 8.88×10⁻⁴ | 8.12×10⁻⁴ |

Agreement to ~1% for equal weight. **But the reliable-detection level does not
move**: all three combinations stay at 5×10⁻³. The injection grid jumps
2×10⁻³ → 5×10⁻³, so a real improvement of ~1.4× cannot resolve into a new tested
level.

At 2×10⁻³, the next lower tested change, the effect depends on the weighting.
Do not summarise this as "pairing improves detection" — only one of the three
weightings does:

| | outbound | equal | inverse-variance | covariance-aware |
|---|---|---|---|---|
| correct-sign detection at 2×10⁻³ | 0.63 | **0.54** | 0.74 | 0.63 |

None approach 0.95. Equal weighting *reduces* detection even though it reduces
scatter, because detection also requires exceeding the empirical threshold and
equal weighting's threshold barely moved (1.822 vs 1.836×10⁻³) while its MAD fell
26% — a heavier-tailed null relative to its core.

The 0.74 figure is the least robust number in this section. Inverse-variance has
the lowest threshold of the three (1.53×10⁻³), and with 23 nulls that threshold is
the ~22nd order statistic, so part of the gain reflects a favourable tail draw
rather than better precision. Rest precision claims on the MAD, not on detection
probability or threshold.

Equal weighting does not degrade *precision*: the two legs have nearly equal
ε-variance (1.20 vs 1.32×10⁻³, 10% apart) despite a 2.4× beam-power difference.

Recommended framing: outbound remains the headline observable, with the paired
inverse-variance estimator reported as a 1.4× noise-floor reduction that does not
cross to the next tested level.

## 5. Solid-Earth-tide upper limit

Better per-burst precision does **not** buy a better tidal limit. This was tested
because a scatter-ratio projection suggested the Deep observable might push the
tidal upper limit below what Paper 1 reports as out of reach at 7.98×10⁻⁴. It
does not.

**Benchmark framing.** Two literature scales are kept distinct and neither is
called "the tidal benchmark". The **Niu SAFOD expected-response scale**,
2.4×10⁻⁷ Pa⁻¹ × 240 Pa = 5.76×10⁻⁵, is the site-specific expectation and the
quantity results are scored against. The **De Fazio observed tidal-response
scale**, ~5×10⁻⁴ to 10⁻³, is historical context only: the largest reported
guided-mode tidal response, measured on an air-filled mine shaft in marble at
500 Hz. It is quoted as a range because the 1973 scan is ambiguous between the
text's "order of 10⁻³" and Figure 3's ~5×10⁻⁴, and nothing depends on which end
is taken.

The zero-injection trials are already a dv/v time series, so the tide model was
regressed directly onto them, with a constant and a linear trend, and inference by
999-realisation surrogate null (`deep_dvv_tidal_fit.py`).

| Observable | Amplitude | σ | p | 95% UL | ×Niu scale |
|---|---|---|---|---|---|
| Deep outbound | +3.64×10⁻⁴ | 4.58×10⁻⁴ | 0.576 | 1.26×10⁻³ | 21.9 |
| Deep return | +6.85×10⁻⁴ | 7.23×10⁻⁴ | 0.378 | 2.10×10⁻³ | 36.5 |
| Deep paired, inverse-variance | +4.60×10⁻⁴ | 4.02×10⁻⁴ | 0.517 | **1.25×10⁻³** | **21.7** |
| Deep paired, covariance-aware | +4.62×10⁻⁴ | 4.08×10⁻⁴ | 0.476 | 1.26×10⁻³ | 21.9 |

Every p-value is a clean null; no tidal signal is detected and none is claimed.
Fitted amplitudes sit close to Paper 1's +3.4×10⁻⁴ and the p-values close to its
0.51, so the two analyses agree.

**The projection was wrong by ~3.7×, for two reasons worth stating in the paper.**
This fit uses 23 held-out bursts against Paper 1's 46, a √2 penalty on its own; and
the tidal amplitude uncertainty runs ~1.5× larger than per-burst scatter predicts,
because over a 22-hour record the diurnal tidal shape is not orthogonal to
instrumental drift. Removing the trend term changes the limits by only 1.01–1.13×,
so the degeneracy inflates the *uncertainty* rather than biasing the amplitude.

### Matched burst count, measured rather than scaled

Because burst count is the obvious confound, the fit was repeated on all 46
bursts (`--population allbursts`). That population is trajectory-contaminated for
a *sensitivity* claim, but legitimate here: trajectory selection is
time-independent and cannot manufacture a tidal signal.

| Observable | amplitude | σ | p | 95% UL | ×Niu scale |
|---|---|---|---|---|---|
| Deep outbound | +5.79×10⁻⁴ | 3.59×10⁻⁴ | 0.701 | 1.28×10⁻³ | 2.57 |
| Deep return | +2.24×10⁻⁴ | 4.82×10⁻⁴ | 0.755 | 1.17×10⁻³ | 2.34 |
| Deep paired, equal weight | +4.01×10⁻⁴ | 3.19×10⁻⁴ | 0.682 | **1.03×10⁻³** | **2.05** |
| Deep paired, inverse-variance | +4.37×10⁻⁴ | 3.13×10⁻⁴ | 0.709 | 1.05×10⁻³ | 2.10 |

**A √N scaling of the 23-burst result was also too optimistic.** It predicted
8.8×10⁻⁴; the measured value is 1.03×10⁻³. Decomposing σ_A = σ_res / L, where L
is the tidal predictor after projecting out the constant and trend:

| | 23 bursts | 46 bursts | ratio |
|---|---|---|---|
| residual scatter σ_res | 8.11×10⁻⁴ | 8.72×10⁻⁴ | 1.076× **worse** |
| effective lever L | 2.016 | 2.788 | 1.383× (√2 = 1.414) |
| σ_A | 4.02×10⁻⁴ | 3.13×10⁻⁴ | 1.286× |

1.383 / 1.076 = 1.286, the observed ratio exactly. So 1/√N failed for two
independent reasons, both small: the added discovery-half bursts are marginally
noisier (robust per-burst scatter 7.76 → 8.03×10⁻⁴), and the effective lever grew
slightly less than √N because the record is truncated at one end (§9). 1/√N was
the wrong model because it assumes both terms are fixed.

**At matched burst count the Deep tidal limit is still ~1.3× worse than Paper 1's
Nano limit of 7.98×10⁻⁴.** Better per-burst precision does not translate into a
better tidal limit. The binding constraint on a tidal search here is the 24-hour
survey duration and its drift degeneracy, not the choice of observable.

One asymmetry worth reporting: at 46 bursts the paired estimators *do* beat
outbound for the tidal limit (1.03 against 1.28×10⁻³, 1.25× better), which they
did not do for the reliable-detection level. Pairing pays off when the metric is
a continuous amplitude and is swallowed when it is discretised onto an injection
grid.

This is the same *style* of regression as Paper 1's registered null, not a re-run:
that fit used the depth-median dv/v chain with common-mode removal over 46 epochs.
The limits are comparable in scale but not identical in construction, and should be
reported that way.

## 6. Return leg — why it is not classified

Of the six pre-registered criteria for calling the return result SNR-limited:

- **Pass** — synthetic geometry/injection validation (worst scale error 0.0040,
  marginally better than outbound's 0.0048)
- **Pass** — monotonic, correct signs at ±10⁻² (detection 1.00 both)
- **Pass** — source jitter absorbed by the intercept, ε = 7×10⁻⁷
- **Weak pass** — scatter worse than outbound by only 1.10×, threshold by 1.66×,
  both well short of the 2.4× beam-power deficit. Directionally consistent,
  smaller than predicted
- **Marginal** — a small one-signed positive bias (+6.8×10⁻⁵ at zero, +2.4%
  scale error at 10⁻²) that outbound does not carry
- **Fails as written** — wrong-observable controls, see below

**The wrong-observable control is not diagnostic.** Every perturbed-trajectory
pass — displaced intercept, both slowness offsets, both aperture lengths, the
split reference — still detects at 10⁻² with probability 1.00. The mechanism is a
flaw in the control, not in the data: the injection shifts the *entire channel
trace* by −ε·T₀(s), so any coherent energy anywhere in the extraction window
carries the same position-dependent shift and is recoverable. The control asks
"is there coherent energy here", not "is the guided mode responsible".

The decisive evidence that the control rather than the analysis is at fault: it
fires identically on outbound and return, at every aperture length and both
reference constructions. A control that behaves the same on the trusted leg and
the leg under test carries no discriminating information.

The one wrong-observable control that does work is the 60–120 Hz band, which
re-derives its own trajectory where no mode was validated: reliable detection
never reached, threshold 8.8×10⁻³ (4.8× worse than primary).

**Consequence:** under a strict reading the return leg cannot be classified,
because criterion 5 is false — and the same criterion would equally condemn
outbound, which is not a defensible scientific reading. Judging criterion 5 on
the band control instead is a post-hoc reinterpretation and must be labelled
exploratory. This is an open decision, not a settled one.

A correctly specified version would inject only inside a time gate around the
trajectory, so the perturbation attaches to the mode rather than the whole trace.
That is a re-run and is exploratory now.

## 7. Defect found and fixed during this analysis

The first complete run produced a corrupted held-out result: `n_recovered`
exceeded the 23 available bursts, and no reliable-detection level was reached at
any injection magnitude on either leg.

Cause: `_inject_population` created `np.random.default_rng(SEED)` separately for
each population. Two generators seeded identically walk the same PCG64 word
stream, so trial ids collided wherever their draw positions aligned — **1383 of
7038 ids duplicated**, about 20%, nothing like birthday chance. The truth table
was a dict keyed on id alone, so the later `allbursts` row silently overwrote the
`heldout` row and ~20% of held-out trials were scored against the wrong injected
value.

The blinded gathers and the recovery estimates were never affected; only the
join was. Re-running `summarize` on the existing data was sufficient.

Fixed three ways: per-population seeding, a group tag making ids unique by
construction, and a composite `(population, leg, band, trial_id)` join key. An
assertion now requires every injection level to contain exactly one trial per
burst — the check that would have caught this immediately.

**A stale, corrupted `deep_dvv_influence.{txt,csv}` was produced before the fix
and has been deleted.** Any figure or number matching thresholds 1.53×10⁻³
(outbound) or 9.94×10⁻³ (return) is from the corrupted run and must be discarded.

## 8. Interpretation ceiling

The recovered quantity is a fractional change in the apparent along-fiber speed
of the selected Deep guided mode. It is not formation Vp or Vs, not fault-zone
stress, pore pressure, permeability, fracture compliance, or tectonic strain, and
it carries no depth resolution.

This experiment measures the sensitivity of the observable *as constructed*.
Attribution of that sensitivity to the guided mode specifically rests on the
prior permutation nulls in `deep_tube_validation` (p = 0.002), not on this test.

## 8b. The tube-wave lineage — what became the paper, and what was dropped

Recorded because the naming hides it. The Deep guided mode **is** the tube-wave
candidate; the paper reports it phase-neutrally, so the lineage is invisible in
the manuscript and easy to think abandoned.

**The search window was chosen on tube-wave grounds.** `deep_tube_validation.py`
scans a slowness band named `TUBE_POSITIVE = (1/1800, 1/1300)` — i.e. 1300–1800
m/s, where a fluid-column tube wave should sit. Every Deep validation product
carries `deep_tube_` in its name. The window is not arbitrary and the manuscript
should say so; at present it does not, which makes a frozen choice look unmotivated.

**A gate test settled whether to build the thesis on it** (`tube_wave_gate.py`,
commit `09e5fe5`). Three claims, verdict `C1 FAIL | C2 PASS | C3 PASS`, against a
pre-set rule of "pivot only if at least C1 and C3 pass":

| | Claim | Result |
|---|---|---|
| C1 | Detectable on both fibers, and they differ | **FAIL** — semblance 0.088 cemented, 0.014 wireline |
| C2 | Tube-wave amplitude images permeable structure | ~~PASS~~ **RETIRED 2026-08-13** — the pass criterion was not a test; see below |
| C3 | It improves the dv/v floor | **PASS** — wireline floor 0.591% on the tube arrival against 2.574% on direct P, a 4.36× gain |

So the commit subject "negative result, do not pursue" refers to **not pivoting
the thesis**, not to the observable. Read claim by claim:

- **C3 became the paper.** That 4.36× floor improvement on the wireline is the
  early, rough version of this analysis's central result, and the gate's 0.591%
  estimate anticipates the final 0.5% reliable-detection level closely.
- **C1's failure is consistent with the finished result.** The two installations
  see different coherent modes rather than the same one at different strengths.
- **C2 is retired as of 2026-08-13.** Followed up in
  `C2_PERMEABILITY_FOLLOWUP.md` (plan) and `C2_PHASE0_RESULTS.md` (results).
  Banerjee & Chatterjee (2021), in the paper collection and still uncited, is the
  method reference.

### Why C2 is retired, and what survives it

The gate's pass criterion was `any(F['drops'].size >= 3 ...)` — at least three
channels below −2σ of a linear log-amplitude trend. That is not a statistical
test. With 1603 wireline channels, chance gives ~36 such channels; the gate found
**seven**, a deficit rather than an excess (binomial P(≤7) = 2.3×10⁻⁹).

Four findings, each sufficient on its own:

1. **σ is not a scatter.** 71 % of the wireline residual variance lies at
   wavelengths > 200 m, so σ measures structure the linear trend failed to
   absorb.
2. **The list is an artefact of the estimator.** Robust MAD instead of the
   standard deviation turns 7 candidates into 42.
3. **Nothing survives a proper null.** Phase-randomised surrogates preserving the
   residual autocorrelation give p = 0.37 (count), 0.73 (depth), 0.29 (step).
4. **Wrong shape, and often wrong sign.** A permeable fracture makes amplitude
   *stay lower below* it. Only 2 of the 7 point that way, none above 0.6σ. Seven
   channels are also only four features at the measured 51 m correlation length.

**The decisive point is that the measurement had no power.** Planted features of
known size are detected at most 34 % of the time even for a 95 % amplitude loss,
so C2 PASS and C2 FAIL were equally uninformative. C2 is *untestable as run*, not
refuted.

**What survives is a quantified prospect.** The factor-3.3 channel-to-channel
scatter is reproducible between burst halves (ρ = 0.922), so it is static
per-channel response, not noise: σ_static = 1.177 against σ_noise = 0.241, i.e.
**96 % of the variance is removable, 5.0× of headroom**. With the static response
divided out, a step would be detectable at an **18 % amplitude loss (1.7 dB)** —
a useful permeable-fracture threshold. That single calibration step, not the four
Phase 1 items, is what decides whether C2 can exist.

**Do not cite the seven depths.** They are a threshold artefact. The cemented
fiber's apparently significant step (p = 0.007) is separately an edge artefact —
it stays pinned to the first allowed breakpoint and reaches p = 0.648 once the
shallowest 300 m are dropped, because log-amplitude is not linear in depth near
the source.

**Unchanged consequence.** The manuscript should still state that the 1300–1800
m/s window was motivated by the tube-wave hypothesis, since that is honest
provenance for a frozen selection region.

## 9. Acquisition inventory, and why 46 and not 48

The burst count is easy to get wrong by a fencepost, and the drop count has four
different legitimate values. Full chain, verified against `awd_manifest.csv`,
`canonical_epoch_stacks_paired_deep_all.npz`, and `paired_stack_job_deep_all.py`:

| Stage | Nano | Deep | Both |
|---|---|---|---|
| GPS drops in `p26.cc9.txt` | — | — | 989 |
| Rows in `awd_manifest.csv` | 988 | 926 | 988 bursts×drops |
| Full PRE/POST window fits inside the data file | 970 | 875 | — |
| **Used by the canonical stacks** | 970 | 875 | **859 common** |

**49 bursts**, median cadence 30.0 min (range 27.6–31.2), spanning 23.96 h.
24 h at half-hour cadence gives 48 *intervals* and therefore 49 *bursts*.

**The Deep fiber stopped recording after burst 45 (22:16:23 UTC).** Bursts 46, 47
and 48 — 22:46, 23:16 and 23:44 UTC on 2026-06-17 — have full Nano coverage and
zero Deep coverage. That is the entire reason the analysis has 46 bursts and not
49. Nothing was discarded for quality.

**Requiring drops on both fibers therefore costs the Deep analysis nothing.** A
Deep-only analysis could not have used those three bursts either; the 46-burst
count is set by Deep's own coverage, not by the pairing requirement.

The 988→970 and 926→875 losses are **window truncation at file boundaries**
(`paired_stack_job_deep_all.py:239` keeps a drop only if its full PRE_S/POST_S
window fits inside the file), not quality rejection. Deep loses more because its
3.0 s post-window is long relative to the `.h5` segmentation. Checks:
859 + 111 Nano-only = 970, and 859 + 16 Deep-only = 875.

This is **not** the same count as Paper 1's "47 of 49 clean", which counts
contamination rather than coverage. Do not use them interchangeably.

Consequence for the tidal fit: because the missing bursts are all at the end, the
series is one-sidedly truncated to 21.98 h (held out) and 22.49 h (all bursts),
not the nominal 24 h. That slightly worsens the diurnal-versus-drift degeneracy.

## 10. Known caveats

- **The "De Fazio scale" of 5×10⁻⁴ — resolved by rendering the page.** The text
  layer drops superscripts, so an earlier reading of "order 10⁻⁴" was wrong. The
  rendered page states: *"The experiment has yielded a velocity change of the
  order of 10⁻³."* Figure 3's velocity-change scale bar is Δc/c = 2×10⁻⁴ with
  φ = 4°, consistent with the paper's 1° ↔ 5×10⁻⁵ phase-to-velocity conversion.
  The project's 5×10⁻⁴ is therefore a defensible peak-to-peak reading of Figure 3,
  and the authors' "order of 10⁻³" is their own rounding of that same trace.
  Beware the trap: the paper's "velocity change (Δc/c) of about 5×10⁻⁵" is the
  phase conversion, not the tidal amplitude.
  **This matters more than a factor of two.** Against 5×10⁻⁴ the Deep upper limit
  of 1.03×10⁻³ is 2.05× the benchmark; against the text's 10⁻³ it is 1.03×, i.e.
  sitting *at* it. Under a 10⁻³ benchmark Paper 1's Nano upper limit of 7.98×10⁻⁴
  becomes 0.8×, which would reverse its headline. Decide which quantity the paper
  cites — Figure 3 peak-to-peak or the authors' order-of-magnitude statement — and
  use it consistently in both papers.
  **What the benchmark actually measures.** De Fazio et al.'s monitored arrival
  was not a body wave: signal-amplitude surveys identified it as "a Rayleigh
  wave-like tube wave along the service shaft surface", guided by a 300 m
  inclined mine shaft at Ogdensburg, New Jersey, at 500 Hz with c ≈ 3000 m/s.
  The canonical tidal benchmark is therefore the tidal response of a *guided mode
  in an engineered opening*, not of bulk rock. This makes comparing it against
  our guided-mode limits like-for-like rather than a category error — a stronger
  justification than the arithmetic. At ~3000 m/s their mode is close to Nano
  (~2950 m/s), not Deep (~1547 m/s), consistent with an air-filled shaft guiding
  a wall-hosted Rayleigh-like mode while a fluid-filled borehole supports a
  Stoneley/tube wave near fluid velocity. The authors also caution that the
  strain-to-velocity relation for such a mode "is difficult to predict
  quantitatively", that expected shaft-direction compression "was not observed",
  and that local topography may matter — reinforcing that guided-mode
  sensitivity belongs to the guiding structure as much as to the rock.
- **The Niu benchmark values are now verified against the source**: stress
  sensitivity 2.4×10⁻⁷ Pa⁻¹ and tidal stress varying within 240 Pa both appear
  verbatim in Niu et al. (2008). Niu et al. further state that the resulting
  travel-time changes, of order 10⁻⁷ s, are close to their measurement error and
  therefore predicted to be undetectable — the same conclusion this work reaches
  by a different route at the same site.

- The return 3–15 Hz trajectory intercept selected at +0.396 s against a search
  ceiling of +0.400 s — a grid-edge hit. Secondary band only; the primary 15–30 Hz
  legs sit well inside. Frozen, so reported rather than re-searched.
- Held-out false-positive rates are granular at 1/23 = 4.3%; the observed 8.7% is
  two trials.
- With n = 23 nulls the 95% threshold is the ~22nd order statistic, so
  differences of 10–20% between estimators are within its sampling noise. The MAD
  is the stable statistic for ranking; the threshold is not.
- The all-46-burst population is trajectory-contaminated and is reported only as
  a check that the held-out threshold is not a small-N artefact. With 46 bursts
  the return leg also reaches 5×10⁻³.

## 11. Files

| File | Role |
|---|---|
| `deep_dvv_injection_recovery.py` | four blinded stages: freeze, inject, recover, summarize |
| `deep_dvv_synthetic_validation.py` | per-leg geometry and injection validation on noiseless synthetics |
| `deep_dvv_influence.py` | leave-one-aperture-out and leave-one-burst-out jackknife |
| `deep_dvv_paired_legs.py` | three paired-leg combinations against the outbound benchmark |
| `deep_dvv_frozen_trajectory.json` | the frozen observable definition |
| `deep_dvv_injection_recovery.sbatch` | `-p serc`, 96 GB, ~16 min |

Blinded gathers (~4.5 GB) are intermediates and live in
`$SCRATCH/deep_dvv_blind/`. The estimator and the injection grid are imported
from `nano_dvv_injection_recovery.py` so that no part of the Nano/Deep difference
can come from a difference in estimator code.
