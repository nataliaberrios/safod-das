# Deep guided-mode injection–recovery — pre-registered interpretation

**Status: frozen.** Written 2026-08-06 by the author, before any sensitivity
result existed.

## Provenance of the freeze

Recorded while job `37777984` was in its `recover` stage. State of the working
directory at the moment this file was written:

| Artefact | Present | Kind |
|---|---|---|
| `deep_dvv_frozen_trajectory.json` | yes | selection product, chosen on the discovery half (`epoch % 2 == 1`) |
| `$SCRATCH/deep_dvv_blind/*.npz` | yes | blinded gathers |
| `deep_dvv_blind_truth.csv` | yes | sealed; not read by the recovery stage |
| `deep_dvv_recovery.csv` | **no** | — |
| `deep_dvv_summary.csv` | **no** | — |
| `deep_dvv_nano_comparison.csv` | **no** | — |
| `deep_dvv_injection_recovery.txt` / `.png` | **no** | — |

No recovered fractional-velocity value, threshold, or detection probability had
been computed or viewed when the rules below were fixed.

The leg designation cites `deep_tube_validation.txt`, which predates this work:
independent even-burst validation beam power is **0.353 outbound** against
**0.147 return** at 15–30 Hz, a factor of ~2.4. That is a property of the
validated observable, not of the sensitivity test.

---

## Pre-registered interpretation of outbound and return Deep results

The outbound 15–30 Hz leg is the primary Deep sensitivity observable. This
designation was made before viewing the injection–recovery results because the
independent validation beam power is 0.353 on the outbound leg and 0.147 on the
return leg, a factor of approximately 2.4 difference.

The return 15–30 Hz leg is treated as an independent secondary replication.
Results for the two legs will be reported separately. Their trials and
thresholds will not be pooled before leg-specific performance is shown.

### Comparison with Nano

The Nano-comparable performance metric is the smallest tested absolute
fractional apparent-velocity change for which at least 95% of trials both:

1. exceed the empirical zero-injection threshold; and
2. have the correct recovered sign.

The less restrictive correct-sign probability will also be reported, but it
will not be used for the primary Nano–Deep comparison.

Interpretation will follow these rules:

- If both Deep legs have Nano-comparable reliable-detection limits below
  1×10⁻², the result may support the general statement that the selected Deep
  guided mode is more sensitive than the Nano apparent-moveout observable.

- If the outbound leg has a reliable-detection limit below 1×10⁻² but the
  return leg does not, the claim will be restricted to the best-observed
  outbound Deep branch. The manuscript will not state that the Deep
  installation as a whole outperforms Nano.

- If the primary outbound threshold is equal to or greater than 1×10⁻², the
  analysis will not claim improved Deep sensitivity, regardless of return-leg
  behavior.

### Classification of a return-leg failure

A weaker return result will be interpreted as signal-to-noise limited only if
all of the following are true:

1. the synthetic geometry and injection validation pass for the return leg;
2. recovery is monotonic and has the correct sign at the largest positive and
   negative injections;
3. zero-injection scatter, aperture-level delay uncertainty, or waveform
   similarity is worse than for the outbound leg in a manner consistent with
   the independently observed beam-power deficit;
4. common-mode jitter is absorbed by the intercept without producing a
   systematic recovered slope;
5. wrong-observable controls do not recover the injected changes;
6. no systematic nonzero bias remains after the frozen quality-control rules
   are applied.

The return result will instead be classified as a trajectory, estimator, or
processing failure if it shows any of the following:

- incorrect signs or large systematic bias at the largest injections;
- nonmonotonic recovery with increasing injection magnitude;
- leakage of common-mode source jitter into the apparent-velocity estimate;
- successful recovery using an intentionally incorrect observable;
- strong dependence on one aperture or one burst;
- disagreement that remains after accounting for measured timing uncertainty
  and usable-aperture count.

If too few return apertures pass the frozen quality controls, the return leg
will be described as insufficiently constrained or uncalibrated, not as a
successful or failed sensitivity measurement.

No frequency band, trajectory, aperture geometry, or quality threshold will be
changed after viewing these results to improve the return-leg outcome.
Alternative post-result configurations, if explored, will be labeled
exploratory and reported separately.

---

## Where each criterion is evaluated

Added at freeze time, so that no criterion depends on a diagnostic invented
after the fact. Nothing here changes a band, trajectory, aperture geometry, or
quality threshold; these are read-outs only.

| Criterion | Evaluated by | Status at freeze |
|---|---|---|
| SNR-limited 1 — return-leg synthetic validation | `deep_dvv_synthetic_validation.py`, run per leg with that leg's frozen trajectory | written before results |
| SNR-limited 2 — monotonicity, sign at extremes | `median_estimated_dvv` per level, `deep_dvv_summary.csv` | already emitted |
| SNR-limited 3 — comparative scatter / delay uncertainty / similarity | `robust_scatter_1p4826_mad`, `estimated_dvv_se`, `median_aperture_correlation` | already emitted |
| SNR-limited 4 — jitter absorbed by intercept | `deep_dvv_controls.csv`, `source_jitter` and `constant_shift` rows | already emitted |
| SNR-limited 5 — wrong-observable controls fail | summary rows for `displaced_intercept`, `slowness_offset_*`, and the 60–120 Hz control band | already emitted |
| SNR-limited 6 — no residual bias after QC | `bias` column per level | already emitted |
| Failure — dependence on one aperture or one burst | `deep_dvv_influence.py`: leave-one-aperture-out and leave-one-burst-out jackknife | **gap at freeze; script written before results were viewed** |
| Failure — usable-aperture count | `n_apertures` per trial, and the fraction of trials falling below `MIN_APERTURES = 6` | already emitted |

The only gap between the pre-registration and the shipped code was the
single-aperture / single-burst influence diagnostic. `deep_dvv_influence.py` was
written and run before any recovered value was inspected; it reads the blinded
gathers and the sealed truth table only after the primary result is fixed.

## Addendum, 2026-08-06 — paired-leg combination

**Frozen before the paired estimator was computed.** At the time of writing, the
single-leg held-out results were known (outbound null threshold 1.84×10⁻³,
reliable detection 5×10⁻³; return 3.05×10⁻³ and 1×10⁻²). No paired-leg number
existed in any form.

Motivation: a reader will ask why two legs sampling the same source were
analysed separately rather than combined. That deserves an answer backed by a
result rather than by speculation.

**Outbound remains the benchmark.** The paired estimator is evaluated *against*
the outbound branch, not as an automatic replacement for it. All three outcomes
are reportable:

- paired beats outbound → paired becomes the best Deep observable;
- paired equals outbound → combining adds no measurable benefit;
- paired is worse than outbound → the noisier return branch degrades the
  cleaner outbound measurement.

Evaluated on the same 23 held-out bursts, the same 15 injection levels, the same
frozen trajectories, and the same threshold and reliable-detection definitions.

### The three combinations, defined in advance

Combination happens at the level of the per-burst estimate, never at the level
of raw traces. Averaging traces across legs before their moveout, polarity,
timing, and waveform differences are separately corrected could smear the signal
for reasons unrelated to sensitivity.

1. **Equal weight** — `(eps_out + eps_return) / 2`. What most readers picture by
   "combine the legs". Expected to be hurt by the noisier return leg.
2. **Inverse-variance weight** — `w = 1 / sigma^2` per leg, treating the legs as
   independent.
3. **Covariance-aware joint** — `(1' C^-1 eps) / (1' C^-1 1)`, where `C` holds
   the two leg variances and their covariance. Both legs record the same source
   impacts, so their errors are expected to be correlated and the independent
   weighting in (2) is expected to be optimistic.

### Where the weights come from

`sigma` and `C` are estimated **leave-one-burst-out from the zero-injection
trials**: the weights applied to burst *b* are computed from the other 22
bursts' null pairs and never from the trial being evaluated. This keeps the
noise regime matched to the population under test, which estimating `C` from the
discovery half would not do — the discovery-half reference is built from 45 other
bursts rather than 22, so its scatter is not comparable.

The weights are therefore frozen per burst and are the same at every injection
level, since they derive only from the null.

### Explicitly not done

- The 23 outbound and 23 return estimates will **not** be concatenated and
  treated as 46 independent bursts. They are 23 paired observations of the same
  23 source bursts.
- Raw traces will not be averaged across legs.

A fourth variant — one joint regression over both legs' apertures with a shared
slope and separate per-leg intercepts — is a different definition of "combine"
and is not run here. If explored later it is exploratory.

### Reported regardless of outcome

The between-leg error correlation, since it determines how much independent
information the second leg can possibly add and therefore explains whichever
outcome occurs.

## Interpretation ceiling, unchanged

The recovered quantity is a fractional change in the apparent along-fiber speed
of the selected Deep guided mode. It is not formation Vp or Vs, not fault-zone
stress, pore pressure, permeability, fracture compliance, or tectonic strain,
and it carries no depth resolution.
