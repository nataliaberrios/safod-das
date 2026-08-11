# Manuscript methods — whole-paper draft

Growing draft following the writing sequence in the project plan (§21: Methods
first). The Deep-specific subsections live in
[`DEEP_DVV_METHODS_DRAFT.md`](DEEP_DVV_METHODS_DRAFT.md) and are cross-referenced
rather than duplicated here, so there is one source per number.

Every value below is transcribed from `awd_manifest.csv`,
`canonical_epoch_stacks_paired_deep_all.npz`, or a named script. Anything I could
not source is marked **[NEEDS SOURCE]** rather than guessed.

Status: §2 Experiment and instrumentation drafted. Remaining Methods subsections
listed at the end.

---

## 2. Experiment and instrumentation

### 2.1 Site and fibres

Two distinct borehole DAS installations at SAFOD recorded the same repeated
surface source. They are not interchangeable observations of one wavefield: the
installations differ in coupling and geometry, and each recovers a different
coherent mode.

The **Nano** fibre is cemented and shallow. It provided 732 channels at 1.266 m
channel spacing, an aperture of 927 m. Data were written as protobuf
(`timeseriesScaler_*.pb`).

The **Deep** fibre is wireline-deployed and reverses at channel 1702, so the
3200 recorded channels of 2.0419 m spacing trace an outbound limb of 3475 m and a
return limb of 3059 m along 6532 m of fibre. Acquisition parameters are recorded
in the filenames (`SAFOD-Deep-10mGL-1000HzFs-2mChDualPulse_*.h5`): **10 m gauge
length, 1000 Hz sample rate, 2 m channel spacing, dual-pulse interrogation.**

Both systems sampled at 1000 Hz.

**[NEEDS SOURCE]** interrogator make and model for each system; Nano gauge length
and pulse configuration; the physical depth of the Deep turnaround; casing and
completion details. Channel coordinate is distance along fibre throughout — the
coordinate-to-depth mapping remains provisional and no result here depends on it.

### 2.2 Source and survey

The source was an accelerated weight drop (AWD) held at a fixed position for the
survey, approximately 15 m horizontally from the Nano wellhead. Drops were fired
in bursts; drops separated by more than 60 s are assigned to different bursts
(`build_manifest.py`, `BURST_GAP_S = 60`).

The survey comprised **49 bursts at a median cadence of 30.0 min (range
27.6–31.2 min), spanning 23.96 h** from 2026-06-16 23:47 to 2026-06-17 23:44 UTC.
Note that 24 h at half-hour cadence gives 48 *intervals* and therefore 49
*bursts*; there are no gaps exceeding 1.5× the median cadence.

Each drop is windowed 0.5 s before to 3.0 s after its GPS time.

### 2.3 Drop inventory and coverage

Four different drop counts are legitimate at different stages, and they must not
be used interchangeably:

| Stage | Nano | Deep | Both |
|---|---:|---:|---:|
| GPS drops listed in `p26.cc9.txt` | — | — | 989 |
| Rows in `awd_manifest.csv` | 988 | 926 | — |
| Full 0.5/3.0 s window fits inside the data file | 970 | 875 | — |
| **Used by the burst stacks** | 970 | 875 | **859 common** |

The **analysed dataset is 859 drops common to both fibres, across 46 bursts.**

Two attritions account for the difference, neither of them quality rejection.

**The Deep fibre stopped recording after burst 45 (22:16:23 UTC).** Bursts 46, 47
and 48 — at 22:46, 23:16 and 23:44 UTC on 2026-06-17 — carry full Nano coverage
and zero Deep coverage. This alone is why the analysis uses 46 bursts and not 49.
Consequently, **restricting the analysis to drops recorded on both fibres costs
the Deep analysis nothing**: a Deep-only analysis could not have used those three
bursts either. The 46-burst count is set by Deep's own coverage, not by the
pairing requirement.

The 988→970 and 926→875 reductions are **window truncation at file boundaries**:
a drop is retained only when its full pre- and post-window falls inside the
recording file (`paired_stack_job_deep_all.py:239`). Deep loses proportionally
more because its 3.0 s post-window is long relative to the `.h5` segmentation.
The inventory closes: 859 common + 111 Nano-only = 970, and 859 + 16 Deep-only =
875.

Because the missing bursts fall at the end of the survey, the analysed series is
truncated at one end to 21.98 h (held-out bursts) or 22.49 h (all bursts) rather
than the nominal 24 h. This matters only for time-series fits over the survey
duration, where it slightly worsens the degeneracy between a diurnal signal and
instrumental drift.

### 2.4 Burst stacking

Drops within a burst are stacked to form one waveform per burst per channel.
Bursts are never silently truncated to an expected drop count; the retained count
per burst is carried through as a stacking weight. Individual drops vary
substantially in timing, amplitude, waveform similarity and signal-to-noise
ratio, whereas burst stacks are highly repeatable — the quantitative repeatability
hierarchy is given in §3.5.

---

## 3. Methods — Nano subsections

### 3.2 Nano mode identification

The Nano coherent mode was characterised by a burst-bootstrap semblance scan over
**signed** slowness and window start (`nano_mode_identification.py`). No phase
name is assigned at any point. Times are measured relative to the `UTC_Date`
column of `p26.cc9.txt`, used purely as an alignment timestamp; the fitted `t0` is
the start of the coherence window, not an onset or a phase pick.

The scan used the 80–440 m aperture at every fourth channel, a 120 ms coherence
window, window starts from −0.100 to +0.320 s in 4 ms steps, and a signed
slowness grid spanning ±1/1500 to ±1/5000 s m⁻¹ in 50 m s⁻¹ increments. Five
bands were scanned independently: 15–25, 25–35, 35–45, 45–60 and 60–80 Hz.
Uncertainty came from 999 bootstrap resamples over bursts (seed 20260802). The
population was **46 bursts and 859 common drops**.

| Band (Hz) | Apparent speed (m s⁻¹) | Bootstrap 95% interval | P(positive ordering preferred) |
|---|---:|---|---:|
| 15–25 | 3300 | [3250, 3350] | 1.000 |
| 25–35 | 2900 | [2900, 2900] | 1.000 |
| 35–45 | 2950 | [2950, 2950] | 1.000 |
| 45–60 | 2950 | [2950, 2950] | 1.000 |
| 60–80 | 2950 | [2950, 2950] | 1.000 |

Positive ordering is preferred in every band with bootstrap probability 1.000.

**Frequency-dependent slowness is resolved**, with a trend of 0.483
[0.414, 0.549] µs m⁻¹ Hz⁻¹. This bears directly on interpretation and should not
be omitted: resolved dispersion is compatible with — but not unique to — guided
propagation, and unresolved dispersion would have been necessary but not
sufficient evidence for direct body-wave propagation. **This test alone cannot
distinguish direct P energy from a weakly dispersive guided or coupled mode.**
The manuscript should therefore describe the Nano observable as a coherent fast
apparent mode and must not present its apparent speed as a formation V_P.

### 3.5 Repeatability metrics

Repeatability was quantified on a fixed, phase-neutral moveout beam
(`nano_hierarchical_repeatability.py`): 30–60 Hz, 80–440 m aperture, 2975 m s⁻¹,
intercept −0.022 s at the first aperture channel — the same frozen trajectory
later used for injection–recovery. The reduced observable is an unnormalised
equal-channel mean after trajectory alignment.

**Population note.** This analysis used **988 Nano-available drops across all 49
bursts**, not the 46-burst common set used for mode identification and
injection–recovery. The counts in this subsection therefore differ from those
elsewhere in Methods by design, because repeatability of the Nano source does not
require Deep coverage. The two populations must not be conflated.

Signal is measured over a −0.040 to +0.160 s window and compared against an
equal-duration noise window at −0.360 to −0.160 s. Delays use a ±12 ms lag
search, with positive delay meaning the observation is later than its
leave-one-out reference.

Quantiles are 16 / 50 / 84%:

| Metric | 16% | 50% | 84% |
|---|---:|---:|---:|
| Within-burst drop NCC, signal window | 0.551 | **0.889** | 0.989 |
| Within-burst drop NCC, noise window | 0.172 | **0.288** | 0.435 |
| Within-burst drop delay (ms) | −0.572 | +0.009 | +0.588 |
| Within-burst relative amplitude | 0.882 | 1.000 | 1.335 |
| Individual-drop beam SNR (dB) | 1.26 | 7.72 | 24.90 |
| Across-burst template NCC | 0.952 | **0.976** | 0.993 |
| Across-burst template delay (ms) | −0.254 | **−0.048** | +0.311 |

All 49 of 49 bursts have median signal NCC exceeding median noise NCC; the exact
one-sided burst-level sign test gives p = 1.78×10⁻¹⁵.

Independent within-burst substacks converge with drop count — median NCC 0.740,
0.782, 0.850 and 0.909 for 1, 2, 4 and 8 drops per substack — which is the
quantitative basis for the claim that stacking converts a variable impact source
into a repeatable observable.

**Two labels to keep distinct**: the median *within-burst drop* delay is +0.009 ms,
whereas the median *across-burst template* delay is −0.048 ms. The −0.048 ms
figure is the across-burst quantity.

### 3.6 Nano injection–recovery

Nano sensitivity was calibrated with the same three-stage blinded design later
applied to Deep — inject, recover, summarise — with a sealed truth table read
only at summarisation (`nano_dvv_injection_recovery.py`).

Burst stacks were bandpassed 30–60 Hz and aligned on the frozen trajectory
(2975 m s⁻¹, intercept −0.022 s) across the 80–440 m aperture at every fourth
channel. Channels were required to be finite in at least 99% of samples with
non-zero dynamic range. References were count-weighted leave-one-burst-out over
all 46 populated epochs.

Gathers were extracted over −0.080 to +0.200 s about the trajectory and delays
measured over −0.020 to +0.140 s with a 12 ms maximum lag, by normalised
cross-correlation with parabolic sub-sample refinement. Channels were retained at
correlation ≥ 0.20, with at least 20 valid channels required. The delay gradient
was fitted against propagation time by Huber IRLS with squared-correlation
weights, and ε taken as the negative slope.

Fifteen injection levels — zero and ±1, 2, 5×10⁻⁴ and ±1, 2, 5×10⁻³ and ±1×10⁻² —
across 46 epochs gave **690 trials**, of which 46 are zero-injection. Injection
applies a per-channel sub-sample shift of −ε·τ, where τ is propagation time
across the aperture.

Because Nano regresses against propagation time across an 80–440 m aperture at
2975 m s⁻¹, its lever arm is **0.121 s** — the quantity against which the Deep
lever arm of 1.398 s is compared.

---

## Remaining Methods subsections

Following the plan's §21 order. Sourcing status noted so nothing inherits
unearned confidence.

| Subsection | Source | Status |
|---|---|---|
| 3.1 Preprocessing and coordinate conventions | `PREPROCESSING.md` | not drafted |
| 3.2 Nano mode identification | `nano_mode_identification.py` | **drafted above, verified against outputs** |
| 3.3 Deep guided-mode identification and split-sample validation | `deep_tube_validation.py`, `DEEP_DVV_METHODS_DRAFT.md` | drafted in the Deep file |
| 3.4 Moveout correction and beam construction | `DEEP_DVV_METHODS_DRAFT.md` | drafted in the Deep file |
| 3.5 Repeatability metrics | `nano_hierarchical_repeatability.py` | **drafted above, verified against outputs** |
| 3.6 Nano injection–recovery | `nano_dvv_injection_recovery.py` | **drafted above** |
| 3.7 Deep local-delay-gradient estimator | `DEEP_DVV_METHODS_DRAFT.md` | drafted in the Deep file |
| 3.8 Blinding, preregistration, software audit | `DEEP_DVV_PREREGISTRATION.md`, `DEEP_DVV_STATUS.md` §7 | drafted in the Deep file |
| 3.9 Reliability definitions | `DEEP_DVV_STATUS.md` §2 | drafted in the Deep file |
| 3.10 Influence diagnostics | `deep_dvv_influence.py` | drafted in the Deep file |
| 3.11 Paired-leg estimator | `deep_dvv_paired_legs.py` | drafted in the Deep file |
| 3.12 Predefined Nano–Deep comparison | `DEEP_DVV_STATUS.md` §2, §3 | drafted in the Deep file |
| Small-signal tidal benchmark | `deep_dvv_tidal_fit.py`, `safod_tides.ipynb` | **placement undecided** — see `DEEP_DVV_STATUS.md` §5 |

Methods is now complete apart from **3.1 preprocessing** and a decision on where
the tidal benchmark lives. Both the Nano and Deep halves are drafted and every
number has been checked against a generated output.
