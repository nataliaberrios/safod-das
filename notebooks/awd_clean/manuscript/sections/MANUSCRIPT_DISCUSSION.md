# Manuscript discussion

Discussion §5.1–5.7, following the plan's §21 step 4. Numbers are drawn from
[`MANUSCRIPT_RESULTS.md`](MANUSCRIPT_RESULTS.md), [`DEEP_DVV_STATUS.md`](DEEP_DVV_STATUS.md)
and [`DEEP_DVV_METHODS_DRAFT.md`](DEEP_DVV_METHODS_DRAFT.md); none are introduced here.

Where a mechanism is proposed rather than measured, it is marked as such. This
paper measures observability and sensitivity; it does not model why the two
installations behave differently.

---

### 5.1 Why the two installations recover different modes

The cemented Nano fibre records a coherent mode near 2950 m s⁻¹, strongest at
30–60 Hz. The wireline Deep fibre records a mode near 1547 m s⁻¹, strongest at
15–30 Hz. Same hole, same source, same nominal 24 h survey — though Deep coverage ends after burst 45, giving ~22.5 h of usable Deep record against Nano's full span.

The most natural reading is a coupling contrast. A cemented fibre is mechanically
tied to the formation and senses strain transmitted through it; a wireline fibre
hangs in borehole fluid and couples preferentially to energy guided by the fluid
column and its walls. The Deep apparent speed of ~1547 m s⁻¹ sits just above
typical borehole-fluid velocity, which is consistent with a tube-wave-like guided
mode — and the validated slowness window (1300–1800 m s⁻¹) was chosen on that
expectation. **This is a proposed mechanism, not a measured one.** No forward
model of either installation was computed here, and the mode's generation point
is not established.

The Nano mode is *also* dispersive — frequency-dependent slowness is resolved at
0.483 [0.414, 0.549] µs m⁻¹ Hz⁻¹, apparent speed falling ~12% from 15 to 80 Hz —
but this does **not** establish guided character, and it should not be presented
as if it did. Resolved dispersion rules out a direct arrival through a
homogeneous medium. It is equally consistent with a refracted or diving wave in a
vertical velocity gradient, where lower frequencies turn deeper in faster
material, and the shallow section at this site has exactly such a gradient: Hole
et al. (2006) map 1.5–3.5 km s⁻¹ in the late Cenozoic sedimentary rocks above the
Salinian basement top at 0.77 km. The Nano aperture lies wholly within that
interval, and 2950 m s⁻¹ falls squarely inside the mapped range.

The honest position is therefore narrower than "both modes are guided". The Deep
observable has positive evidence of guided character — a speed near borehole-fluid
velocity, coherent on both limbs, surviving channel-order permutation. For Nano,
a refracted P interpretation and a guided or coupled interpretation both remain
viable, and this experiment does not separate them. Reporting the Nano
observable by its apparent speed and frequency content, without a phase name, is
a deliberate choice to avoid deciding a question the data leave open.

#### Guided modes are the norm in this measurement class, not an anomaly

It is worth noting how much company these observables keep. Ellsworth & Malin
document fault-zone-guided waves crossing the SAFOD borehole at 2.7 km. And the
canonical tidal-velocity experiment, De Fazio et al. (1973), monitored what its
authors identify as "a Rayleigh wave-like tube wave along the service shaft
surface" — a mode guided by a 300 m inclined mine shaft at 500 Hz, with
c ≈ 3000 m s⁻¹.

A numerical coincidence here is worth defusing rather than exploiting. De Fazio's
mode also travels at ~3000 m s⁻¹, close to the Nano apparent speed, but the two
numbers mean different things. De Fazio et al. do not quote a P velocity, but
they do give an elastic constant "of the order of 5 × 10⁵ bar" (50 GPa) for the
Franklin marble hosting their shaft. With a marble density of ~2700 kg m⁻³ that
implies Vp ≈ 5000–5800 m s⁻¹ depending on whether the quoted modulus is Young's
or bulk — either way, far above 3000 m s⁻¹. Read as a bulk modulus with μ ≈ 30 GPa
it also gives a Rayleigh velocity of ~3070 m s⁻¹, essentially their observed
value. Their 3000 m s⁻¹ is therefore roughly half the P velocity of their own
host rock and sits at its Rayleigh velocity, which is why their Rayleigh-like
identification is self-consistent. At SAFOD the same figure lies inside the
mapped shallow Vp range. Equal apparent speeds in different rock, at 500 Hz
against 30–60 Hz, in an air-filled shaft against a cemented borehole, imply no
kinship between the modes.

What does carry over is the interpretive category. **The comparison against the De Fazio benchmark is
closer than it appears.** We are not comparing a guided-mode sensitivity against
a formation-velocity prediction; we are comparing a guided-mode tidal response
against another guided-mode tidal response, measured with a different source in a
different opening. That makes it **more directly comparable in observable type
than a bulk-rock velocity benchmark would be** — not a like-for-like match, since
the rock, frequency, geometry, coupling and guiding structure all differ, but a
comparison of the same kind of quantity.

The same lineage carries a caution the original authors stated plainly. De Fazio
et al. note that the effect of strain on such a mode "is difficult to predict
quantitatively", that the compression they expected in the shaft direction "was
not observed", and that local topography "may have to be taken into account". A
guided mode's sensitivity to strain is a property of the guiding structure as
much as of the rock. That is precisely why this paper reports apparent
guided-mode velocity changes and does not convert them into formation properties.

### 5.2 A repeatable source is not a precision source

The accelerated weight drop is unambiguously repeatable at the level of waveform
similarity. Median within-burst drop correlation is 0.889 in the signal window
against 0.288 in an equal-duration noise window, and all 49 of 49 bursts have
median signal correlation exceeding median noise correlation (p = 1.78×10⁻¹⁵).
The coherent mode is present in individual impacts, not manufactured by stacking.

It is not, however, precise. Individual-drop beam SNR spans 1.26 to 24.90 dB
between the 16th and 84th percentiles — a factor of ~20 in power — and relative
amplitude varies by roughly ±25%. Per-burst timing repeatability, inferred from
the injection–recovery scatter, is 0.225 ms for Nano and 1.68 ms for Deep.

Stacking is what converts the first property into the second. Independent
within-burst substacks converge from median correlation 0.740 at one drop to
0.909 at eight, and full burst stacks reach 0.976 against a leave-one-burst-out
reference. This is why the sensitivity results are burst-level results, and why
per-shot Deep statistics (correlation 0.054) and per-burst Deep statistics are
not in conflict — they describe different objects.

The distinction matters for comparison with earlier work. De Fazio et al. (1973)
used a continuous, phase-locked source, so their precision was set by oscillator
stability rather than by source variability. Niu et al. (2008) fired four times
per second and stacked to one high-SNR record every ~45 min. A repeated impact
source occupies a different regime: it is cheap and portable, its individual
realisations are noisy, and its usable precision is bought entirely by stacking.

### 5.3 A longer lever arm does not convert one-for-one

This is the paper's central quantitative result and the one that generalises.

The delay-gradient estimator recovers a fractional speed change from the slope of
delay against reference travel time, so the precision of ε scales as

> ε precision ≈ (per-burst timing repeatability) / (lever arm in travel time).

Deep's geometry is dramatically better on the denominator. Its regression runs
over aperture centres spanning 1.398 s of travel time, against 0.121 s for Nano —
**11.5× longer**. Taken alone, that predicts an order-of-magnitude sensitivity
gain.

The realised gain is about a factor of two: reliable detection improves from
1.0% to 0.5%, and the empirical null threshold from 0.325% to 0.184%. The
numerator explains the shortfall. Deep's per-burst timing repeatability is
1.68 ms against Nano's 0.225 ms, **7.4× worse**, and 11.5 / 7.4 = 1.55, which is
exactly the observed ratio of null scatter (1.86×10⁻³ to 1.20×10⁻³).

The design lesson is not that the Deep installation is disappointing. It is that
**fractional velocity sensitivity is a ratio, and extending the geometry buys
nothing if the extension degrades timing repeatability proportionally.** A
wireline deployment reaching deeper into the hole delivers a long propagation
path and, through poorer coupling, a noisier arrival time. Roughly four fifths of
the geometric advantage is consumed by that trade.

Stated as design guidance: for a target fractional sensitivity, an experiment
must specify both terms. Reporting aperture or depth alone does not constrain
monitoring performance.

### 5.4 Why the return branch does not support an installation-wide claim

Sensitivity is not uniform across the wireline installation. The outbound branch
reaches a reliable-detection level of 0.5%; the return branch reaches 1.0%,
equal to Nano rather than better. Under the pre-registered decision rule this
restricts the improved-sensitivity claim to the best-observed outbound branch,
and the manuscript does not state that the Deep installation as a whole
outperforms Nano.

The deficit is smaller than the headline levels suggest, and the reason is worth
reporting. The return branch's *core* precision is nearly as good as outbound's —
robust null scatter 1.32×10⁻³ against 1.20×10⁻³, only 1.10× worse — despite an
independent validation beam power 2.4× lower. Its threshold penalty (0.305%
against 0.184%, 1.66×) is therefore driven mostly by the tail of a 23-sample null
distribution, where the 95th percentile is the ~22nd order statistic. The two
branches differ far less in precision than in the discrete level they happen to
cross.

**The cause of the return deficit is left unclassified, deliberately.** Five of
the six pre-registered criteria for attributing it to signal-to-noise are
satisfied: the return geometry passes noiseless synthetic validation (marginally
better than outbound), recovery is monotonic with correct signs at the largest
injections, common-mode timing is absorbed by the intercept, aperture yield is
full, and no single aperture or burst dominates. The sixth criterion — that
wrong-observable controls should fail — could not be evaluated, because those
controls proved non-diagnostic (§5.5). Rather than reinterpret a frozen criterion
after seeing the result, the branch is reported as unclassified.

### 5.5 What these measurements do and do not represent

The recovered quantity is a fractional change in the apparent along-fibre speed
of a selected guided mode. It is not a formation V_P or V_S change, and it is not
a measure of stress, pore pressure, permeability, fracture compliance, or
tectonic strain. Converting it into any of those would require guided-wave
physics and forward modelling not attempted here.

It carries no depth resolution. Channel coordinate is distance along fibre; the
coordinate-to-depth mapping remains provisional, and no result depends on it.

The injection–recovery experiment measures the sensitivity of the observable *as
constructed*. Attribution of that sensitivity to the guided mode specifically
rests on the earlier split-sample and channel-order permutation validation — zero
of 499 permutations reached the observed beam power in any of four leg/band
tests — and not on the injection experiment itself.

One control limitation should be reported rather than quietly dropped. The
perturbed-trajectory wrong-observable tests — displaced intercept and offset
slowness — recovered the injected changes rather than failing. The mechanism is a
flaw in the control, not in the data: the injection shifts each entire channel
trace by a position-dependent delay, so any coherent energy anywhere in the
extraction window carries the same imposed gradient and remains recoverable.
Those controls therefore test whether coherent energy exists in the window, not
whether the selected mode is uniquely responsible, and they behave identically on
both legs. The 60–120 Hz band, where no mode was validated and no reliable
detection level is reached, is the wrong-observable control that does work.

### 5.6 Implications for repeated-source borehole DAS monitoring

Four points follow for anyone designing a comparable experiment.

**Installation is an experimental variable, not a deployment detail.** Two fibres
in one hole, recording one source, returned different coherent modes with
different frequency content, different apparent speeds and a factor of two
difference in monitoring sensitivity. Installation should be specified and
justified alongside source and geometry.

**Geometry and timing repeatability must be reported together.** A long
propagation path is only half the ratio that sets fractional sensitivity.

**Stacking design determines what is measurable.** With a repeated impact source,
the number of drops per burst sets the usable precision; the convergence curve,
not a single-shot statistic, is the relevant specification.

**Spatially distributed branches improve precision but not necessarily detection.**
Combining the outbound and return branches reduced robust null scatter by ~1.5×,
close to what independent averaging predicts for their measured error correlation
of ρ = 0.121. It did not lower the smallest tested reliable-detection level,
which remained 0.5% for every weighting. Precision gains show up continuously;
detection levels move in discrete steps.

### 5.7 Small-signal benchmark and the limits of this survey

Two distinct literature scales appear below and must not be conflated. The
**Niu SAFOD expected-response scale** is what a tidal response *should* be at this
site; the **De Fazio observed tidal-response scale** is the largest such response
anyone has reported, measured elsewhere on a different guided mode. Neither is
"the tidal benchmark".

The Niu SAFOD expected-response scale follows from Niu et al. (2008):
a barometric stress sensitivity of 2.4×10⁻⁷ Pa⁻¹ against a calculated tidal
stress variation of ~240 Pa gives ~5.76×10⁻⁵, or 0.00576%. This is an
order-of-magnitude benchmark for a small signal, not a prediction for the modes
measured here.

That scale lies roughly 87× below the best demonstrated reliable-detection level
of 0.5%, and ~32× below the Deep outbound null threshold. The measured Deep tidal
upper limit of 1.03×10⁻³ sits ~18× above it; Paper 1's Nano upper limit of
7.98×10⁻⁴ sits ~14× above it. Under idealised
independent stacking, reaching it from the Deep single-burst null would require
of order 1000 independent burst stacks; the survey collected 46.

Fitting the tide model directly to the per-burst velocity estimates gives a 95%
upper limit of 1.03×10⁻³ on 46 bursts, with surrogate p-values between 0.66 and
0.76 — a clean null, and no tidal response is claimed. Notably, **the more
sensitive observable did not produce a better tidal limit**: 1.03×10⁻³ is ~1.3×
worse than the Nano upper limit obtained from the same survey, despite Deep's
better per-burst precision. Two effects account for this. Amplitude uncertainty
is inflated ~1.5× beyond what per-burst scatter predicts, because over a 22 h
record the diurnal tidal shape is not orthogonal to instrumental drift; and the
record is truncated at one end, the Deep fibre having stopped after burst 45.

The binding constraint on a tidal search with this design is therefore **survey
duration and its drift degeneracy, not the choice of observable**. A longer
record that separates diurnal signal from drift would do more than a more
sensitive mode.
