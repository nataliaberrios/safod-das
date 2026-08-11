# Manuscript introduction

Drafted from papers read in `/home/groups/ettore88/nberrios/planpapers/`. Every
claim attributed to prior work below was checked against the PDF text, not cited
from memory. Extraction: `ml system poppler/0.47.0` then `pdftotext`.

---

## 1. Introduction

Repeated active-source seismic monitoring can resolve small changes in crustal
velocity, and borehole deployment removes much of the near-surface noise and
weather sensitivity that limits surface experiments. Distributed acoustic sensing
(DAS) makes such deployments dense and cheap: a single fibre returns thousands of
channels along the borehole, where a conventional array returns a handful. What
DAS does not do is remove the question of how the fibre is coupled to the rock,
and that question turns out to control what the instrument can measure.

At Parkfield, both halves of this problem have prior art at the same site. Niu
et al. (2008) monitored velocity at SAFOD with a repeated active source and
measured a stress sensitivity of 2.4 × 10⁻⁷ Pa⁻¹, calibrated against ~1.3 kPa of
barometric loading. They also computed the solid-Earth tidal stress at the site,
found it varies within 240 Pa, and concluded that the resulting travel-time
changes — of order 10⁻⁷ s — were close to their measurement error and therefore
predicted to be undetectable. Earlier, De Fazio et al. (1973) had monitored
seismic phase velocity continuously with a phase-locked shaker to a precision
finer than one part in 10⁴, and detected periodic velocity variations they
attributed to solid-Earth tides, of an order substantially larger than laboratory
stress sensitivities on small rock samples would predict. Together these frame
the small-signal regime: a tidal velocity response is real but sits near or below
what a repeated-source experiment can resolve.

The wavefield at SAFOD is also known not to be a simple body-wave field.
Ellsworth & Malin showed that the SAFOD borehole crosses a channel of damaged,
low-velocity fault-zone rock at 2.7 km depth which supports Love-type and
Rayleigh-type fault-zone-guided waves along with a leaky mode travelling in the
fault core. Guided propagation is therefore expected here, and an apparent
velocity measured along a borehole should not be assumed to be a formation
velocity. Lellouch et al. (2019) subsequently recorded the same site with
downhole DAS and estimated P- and S-wave velocities from passive records,
demonstrating that fibre in this hole resolves structure but also that the
relationship between a DAS observable and a formation property requires care.

What has not been established is how the *installation* — cemented against
wireline, shallow against deep — changes the answer, and what that costs or buys
for time-lapse monitoring. Two fibres in the same borehole, recording the same
source, provide a controlled test of exactly this. If they return the same
observable with the same sensitivity, installation is a logistical choice. If
they do not, installation belongs in the experimental design alongside source and
geometry.

A second question follows from the geometry. The natural estimator for a
fractional velocity change measures the rate at which delay accumulates with
propagation time, so its precision improves as the observed propagation path
lengthens. A fibre reaching further into the borehole, or a slower mode occupying
more travel time over the same distance, should therefore be more sensitive. That
argument treats timing precision as fixed. It need not be: the same deployment
choices that extend the path may also degrade coupling, waveform repeatability
and arrival-time stability. Whether a longer lever arm converts into better
sensitivity is an empirical question, and it is one that a repeated source
answers directly.

Answering it requires separating two things that are easily conflated. Whether a
coherent mode is *observable* — repeatable, directional, distinguishable from
noise — is a different question from whether its apparent velocity has a known
physical interpretation. This paper addresses the first and is deliberately
conservative about the second: modes are described by their apparent speed and
frequency content, no phase names are assigned, and no apparent velocity is
converted into a formation property. Sensitivity is calibrated empirically, by
injecting known velocity changes into real data and recovering them blind,
because that approach carries the actual source variability, noise and estimator
behaviour of the experiment rather than an idealised error model.

The source itself sets a further constraint. De Fazio et al. used a continuous,
phase-locked shaker, so their precision was governed by oscillator stability
rather than by source variability; Niu et al. fired four times per second and
stacked to one high-signal record roughly every 45 minutes. A repeated impact
source is portable and inexpensive, but its individual realisations vary in
timing, amplitude and waveform, and its usable precision must be recovered by
stacking. Quantifying that hierarchy — from single impacts, through substacks, to
burst stacks — is a prerequisite for interpreting any sensitivity number obtained
from it.

We therefore ask four questions of a 24-hour repeated accelerated weight-drop
survey recorded simultaneously on two borehole DAS installations at SAFOD:

1. What coherent modes does each installation recover?
2. How repeatable is the source, from individual impacts to burst stacks?
3. What fractional apparent-velocity changes are empirically resolvable through
   each observable?
4. Does the longer propagation time of the deeper installation's mode deliver
   correspondingly better fractional sensitivity?

---

## Citation notes and verification status

| Claim | Source | Status |
|---|---|---|
| Stress sensitivity 2.4 × 10⁻⁷ Pa⁻¹ | Niu et al. (2008) | **verified in text** |
| Tidal stress varies within 240 Pa | Niu et al. (2008) | **verified in text** |
| Tidal travel-time changes ~10⁻⁷ s, predicted undetectable | Niu et al. (2008) | **verified in text** |
| Precision finer than 1 part in 10⁴ | De Fazio et al. (1973) | **verified in text** |
| Observed tidal velocity change order 10⁻⁴, exceeding lab prediction of 10⁻⁵–10⁻⁶ | De Fazio et al. (1973) | order verified; **exact value not recoverable from this scan** |
| Guided-wave channel at 2.7 km, FL/FR plus a leaky fault-core mode | Ellsworth & Malin | **verified in abstract** |
| Downhole DAS P/S velocity estimation at SAFOD | Lellouch et al. (2019, JGR) | **verified in abstract** |

**Warning on the De Fazio benchmark.** The project quotes a "De Fazio scale" of
5 × 10⁻⁴. The 1973 scan's text layer drops superscripts, so the observed
amplitude reads as "of the order of 10⁻[illegible]". What is recoverable is that
the observed change exceeds the laboratory prediction of 10⁻⁵–10⁻⁶ and that the
authors emphasise in-situ sensitivity is "much greater" than laboratory values,
which places it at order 10⁻⁴.

The specific value 5 × 10⁻⁴ is **not confirmed**. Note also the trap: the paper
does contain "a velocity change (Δc/c) of about 5 × 10⁻⁵", but that is the
instrument's phase-to-velocity conversion — 1° of phase corresponds to 5 × 10⁻⁵ —
and not the tidal amplitude. Anyone re-checking the benchmark could easily take
the wrong number.

This affects only the multiplier, not any conclusion. The measured Deep upper
limit of 1.03 × 10⁻³ is 2.05× a 5 × 10⁻⁴ benchmark and ~10× a 1 × 10⁻⁴ benchmark;
"cannot resolve a De Fazio-scale tidal change" holds either way, and holds more
strongly under the smaller value. Before submission the figure should be settled
from a clean copy of the paper or from its Figure 3.
