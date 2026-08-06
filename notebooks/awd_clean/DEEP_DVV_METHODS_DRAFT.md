# Deep guided-mode sensitivity — methods draft

Draft prose for the manuscript. Every parameter below is transcribed from
`deep_dvv_frozen_trajectory.json` and the module constants in
`deep_dvv_injection_recovery.py`; results are in
[`DEEP_DVV_STATUS.md`](DEEP_DVV_STATUS.md), the freeze in
[`DEEP_DVV_PREREGISTRATION.md`](DEEP_DVV_PREREGISTRATION.md).

Not yet edited for journal style or length. Numbers are authoritative.

---

## Observable and its selection

We measured the sensitivity of the repeatable slow guided mode on the Deep
fibre, treating it as a phase-neutral time-lapse observable. The recovered
quantity is a fractional change in the mode's apparent along-fibre speed; it is
not interpreted as a formation velocity change.

Data are the 49 count-weighted burst stacks of the June 2026 accelerated
weight-drop survey (2026-06-16 23:47 to 2026-06-17 23:47 UTC), of which 46
contain drops common to both fibres, at 1000 Hz sample rate and 2.0419 m channel
spacing over 3200 channels. The fibre reverses at channel 1702, giving an
outbound and a return leg; the return leg is indexed in reverse so that the leg
coordinate measures distance from the surface end on both. Analysis is confined
to 200–3000 m along each leg, the range spanned by the candidate apertures
validated previously, decimated by a factor of six to 12.25 m spacing (229
channels per leg).

Bursts were split by parity, inheriting the convention of the earlier
split-sample validation: the 23 odd-indexed bursts form the discovery set and
the 23 even-indexed bursts are held out. Every selection below was made on the
discovery set alone and then frozen.

For each leg and frequency band we selected one trajectory by grid search on the
discovery-set stack, maximising the semblance of a linear moveout over the full
leg. Slowness was confined to the independently validated 1300–1800 m/s window
and intercept searched from −0.10 to +0.40 s in 2 ms steps, with semblance
evaluated over a ±35 ms window at 15–30 Hz. The two legs, which share no bursts,
selected 1544.6 m/s at an intercept of +0.100 s (outbound) and 1549.7 m/s at
+0.346 s (return) — agreement in speed to 0.3%.

The primary band is 15–30 Hz. A 3–15 Hz band is reported as a secondary
robustness test, and a 60–120 Hz band, in which no coherent mode was validated,
serves as a wrong-observable control.

## Estimator

A change in propagation speed produces a delay that accumulates with travel
time, whereas an error in source timing shifts the whole wavefield equally. For
each burst we therefore measured the delay of the mode at a series of positions
along the frozen trajectory and fitted

    dt(s) = a - eps * T0(s),

where T0(s) is the reference travel time at position s, a absorbs source and
trigger timing error, and eps is the fractional change in apparent velocity.

Each leg was divided into 400 m apertures at 200 m spacing, giving twelve per
leg. Within an aperture, channels were moveout-corrected onto the frozen
trajectory and averaged into a local beam trace after normalising each channel to
unit RMS, a time-invariant scaling that cannot bias a delay but prevents a
shallow high-amplitude interval from dominating the beam. Delays were estimated
by normalised cross-correlation with parabolic sub-sample refinement against a
reference beam, over a −80 to +160 ms window with a maximum lag of 30 ms.

Apertures were retained when the correlation reached 0.30 and the delay stayed
within 90% of the maximum lag; at least six retained apertures were required. The
gradient was fitted by Huber iteratively-reweighted least squares with weights
equal to the squared correlation, and eps taken as the negative of the slope. All
thresholds were fixed before evaluation.

References were constructed leave-one-burst-out within the held-out population
and count-weighted, so no burst contributed to the reference against which it was
evaluated.

## Injection and recovery

Sensitivity was measured by injecting known fractional velocity changes into real
held-out data and recovering them blind. For an injected change eps_inj, each
channel trace was shifted by −eps_inj·T0(s) before beamforming, so the
perturbation passed through the entire pipeline.

Fifteen levels were used — 0 and ±1×10⁻⁴, ±2×10⁻⁴, ±5×10⁻⁴, ±1×10⁻³, ±2×10⁻³,
±5×10⁻³, ±1×10⁻² — identical to those used for the Nano observable, giving 345
trials per leg and band. Two further sets of 23 trials per leg tested timing
nuisance: a 5 ms shift applied equally to all channels, and a per-burst random
shift drawn with a 0.30 ms standard deviation, the common-mode source timing
scatter measured independently for this survey.

Injection and recovery were separated into stages. The injection stage wrote
randomised perturbed gathers with identifiers carrying no information about the
injected value, and a sealed truth table written separately; the recovery stage
read only the blinded gathers. Truth was joined to estimates only at the
summary stage.

Recovery performance was summarised against the zero-injection trials. Estimates
were centred on the null median, and the empirical two-sided threshold taken as
the 95th percentile of the absolute centred null. The reported detection limit is
the smallest tested magnitude at which at least 95% of trials both exceeded that
threshold and carried the correct sign, in both directions — the same definition
applied to the Nano observable, so that the two are directly comparable. The less
restrictive correct-sign probability is also reported but is not used for the
comparison.

To eliminate estimator differences as a source of the Nano/Deep contrast, the
delay estimator, the robust regression, and the injection grid were imported
directly from the Nano analysis rather than reimplemented.

## Controls

Six controls accompany the primary result.

*Zero injection* supplies the empirical null and false-positive rate. *Constant
and random timing shifts* test whether source timing error leaks into the
velocity estimate. *Outbound and return* were analysed independently and never
pooled before leg-specific performance was established. *Reference sensitivity*
compares the leave-one-burst-out reference against a split-half reference.
*Aperture length* was varied over 200, 400 and 800 m. A *wrong-observable* test
repeated the analysis with a displaced intercept, with slowness offset by
±5×10⁻⁵ s/m, and in the 60–120 Hz band where no mode was validated.

Two further diagnostics were run. A synthetic guided mode on the real geometry,
perturbed both by re-synthesis at a genuinely altered speed and by the pipeline's
own injection routine, confirms that the injection is a faithful proxy for a
medium change rather than a shift the estimator is guaranteed to find. A
leave-one-aperture-out and leave-one-burst-out jackknife bounds the influence of
any single measurement.

## Paired-leg combination

Because the two legs sample the same source bursts, we tested whether combining
them improves on the better leg alone. Combination was performed on the
per-burst estimates, never on raw traces, and evaluated on the same held-out
bursts and injection levels. Three weightings were compared: equal, inverse
variance, and a covariance-aware weighting `(1'C⁻¹eps)/(1'C⁻¹1)` accounting for
correlation between the legs. Variances and covariance were estimated
leave-one-burst-out from the zero-injection trials, so the weights applied to a
burst never derive from the trial being evaluated. The outbound leg was
designated the benchmark in advance; the combination was assessed against it
rather than substituted for it.

## Limits of interpretation

The recovered quantity is a fractional change in the apparent along-fibre speed
of the selected guided mode. It carries no depth resolution and is not a
formation Vp or Vs change, nor a measure of stress, pore pressure, permeability,
fracture compliance, or tectonic strain. Converting it would require guided-wave
physics and forward modelling not attempted here.

The experiment measures the sensitivity of the observable as constructed.
Attribution of that sensitivity to the guided mode specifically rests on the
channel-order permutation nulls of the earlier validation (p = 0.002), not on
this test.

Two caveats are carried explicitly. The perturbed-trajectory wrong-observable
control proved non-diagnostic: because the injection shifts entire channel
traces, any coherent energy within the extraction window carries the same
position-dependent shift and remains recoverable, and the control accordingly
behaves identically on both legs. The return leg is therefore reported as
unclassified rather than as signal-to-noise limited, notwithstanding that it
satisfies the other pre-registered criteria.
