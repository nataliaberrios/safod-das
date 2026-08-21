# Literature review: robust ambient-noise CC, and robust F-K filtering in particular

Written 2026-08-19 for the Figure 7c recovery goal. Every paper below was
retrieved and read through the Scite full-text index or the publisher; nothing
here is cited from memory. Where a paper's method has already been tried in this
tree, the script is named. Where it has not, that is stated as a gap.

Companion documents: `AMBIENT_LOWK_MECHANISM.md` (what the contaminant is),
`Ambient_FK_QC_workflow.ipynb` (why the fixed 2.5-4.5 km/s fan was rejected),
`FIG7C_MULTIDAY_RESULT.md` (the multi-day negative).

> **Reconciled later the same day (2026-08-19).** Two of the four gaps in §4 were
> closed within hours of this review being written: adaptive f-k was run (§1, five
> configurations, all negative) and the illumination prerequisite was tested (§4.3,
> and it **fails**, which makes it the binding constraint). Sections 1, 3, 4, 4a
> and 5 carry inline updates. Read §5 first — it is the current summary. Where this
> file and `AMBIENT_LOWK_MECHANISM.md` describe the same two findings, they agree
> on precedence: illumination first, static fixed-k pattern second.

---

## 1. The finding that reframes our F-K work

**Isken, Vasyura-Bathke & Dahm (2022).** *De-noising distributed acoustic sensing
data using an adaptive frequency-wavenumber filter.* Geophysical Journal
International 231(2), 944-949. https://doi.org/10.1093/gji/ggac229 (open access)

This is the most directly relevant paper found, and it is relevant because it
does something we never did. Their filter is **data-adaptive and
non-destructive**: the mask is built from the observed amplitude spectrum raised
to a power alpha, so the filter "enhances the energy of coherent waves
(wavenumbers) within individual filter windows and no energy will be suppressed
to the stop band." It is applied in **sliding t-x windows**, with overlaps
"tapered in t-x domain with a Bartlett (triangular) taper in order to fade
equally between overlapping windows", so "filtering in sliding windows makes the
operation adapt to locally present coherent signals." They demonstrate it on
"a noisy data set that was recorded in a **vertical borehole observatory**
showing active and passive seismic phase arrivals" — our geometry.

They also report a design choice worth copying and one worth knowing about:

- They **deliberately do not smooth the amplitude spectrum**, unlike the InSAR
  practice of Goldstein & Werner (1998) and Baran et al. (2003), because
  "smoothing the amplitude spectrum (particularly in a small window) distorts
  the signal amplitudes and narrows the filter band."
- A **normalized** variant (NAFK) normalises the spectrum coefficients, which
  "will retain the amplitude of the most dominant f-k component of the input."

**How this bears on our failure.** Everything in our F-K family imposes a
*fixed velocity fan* — `ambient_signed_fk_v2.VELOCITIES`, the frozen 2.5-4.5 km/s
selection region, and the brick-wall/raised-cosine masks in
`ambient_fk_taper_test.py`. An adaptive filter makes no velocity assumption at
all. So the honest answer to "have you tried hard enough at F-K filtering" is:
we tried one *family* hard (fixed-fan, two mask shapes, mask-sensitivity sweeps,
kernel-residual checks) and never tried the adaptive family.

**But it predicts failure here, for a specific reason.** An adaptive filter
enhances the *dominant coherent* f-k component, and the NAFK variant explicitly
"retain[s] the amplitude of the most dominant f-k component." In our data the
dominant coherent component is the static fixed-k pattern holding ~39 % of the
in-band energy (`AMBIENT_LOWK_MECHANISM.md` §1). AFK/NAFK applied directly would
therefore **amplify the contaminant**, not the arrival. This is a testable
prediction, not a reason to skip the experiment: the correct ordering is

    static per-channel response removal (spatial)  ->  AFK/NAFK (adaptive)

and only in that order. When this section was written this was the single
highest-value untried method in the review.

### It has now been run, and the prediction registered above was confirmed

**Updated 2026-08-19.** Adaptive f-k was implemented and run in **five
configurations** (`afk_recovery.sbatch`, `afk_aggregate.sh`; the results table is
generated from the aggregates by `fig7c_negative/append_afk_results.py` so it
cannot drift from them). **None passed the five predeclared gates.**

| configuration | peak | at (m/s) | p | pedestal |
|---|---:|---:|---:|---:|
| AFK α=1, **no** static removal | 12.53 | 5900 | **0.0060** | **+0.987** |
| median common mode, no AFK | 1.91 | 5900 | 0.8796 | +0.820 |
| median common mode + AFK α=1 | 1.21 | 5900 | 0.5199 | +0.766 |
| median common mode + AFK α=2 | 1.12 | 5900 | 0.3465 | +0.737 |
| median common mode + rank-2 + AFK α=1 | 1.08 | 3300 | 0.9392 | **−0.036** |

Two rows carry the lesson.

**Row 1 — the filter manufactured significance, exactly as predicted here.** With
no prior static removal, AFK reached **p = 0.0060**, which read in isolation is a
detection. It is not one. It carries the *worst* pedestal diagnostic of the five
(+0.987) and peaks at 5,900 m/s — the scan ceiling — not 3,200. The filter raised
the score 6.6-fold (1.91 → 12.53) by amplifying the fixed-k pedestal, and the
amplified pedestal then cleared a null that structurally could not see the
amplification.

**The generalisable mechanism, which applies to any filtered ambient result in
this tree:** the receiver-order null permutes the **finished gather**, so an
operator applied *before* the gather is formed sits outside its own null. Only the
predeclared gates caught this. Any future pre-gather operator owes an
**input-level** null.

**Row 5 — the cleanest configuration.** Common mode + rank-2 + AFK is the only row
that drives the pedestal to essentially zero (−0.036) and the only one whose peak
lands in the physical range (3,300 m/s). It gives **p = 0.9392**. Repair the
statistic and the peak moves where it should, and it is still worth nothing —
the same behaviour the stack-convergence experiment found independently.

So the honest verdict on the criticism that the F-K work was never pushed to the
adaptive family: the criticism **was** correct as to breadth, the gap **is** now
closed, and closing it did not recover the arrival.

---

## 2. What successful borehole ambient-noise body-wave retrieval required

**Behm (2016).** *Feasibility of borehole ambient noise interferometry for
permanent reservoir monitoring.* Geophysical Prospecting 65(2), 563-580.
https://doi.org/10.1111/1365-2478.12424

Reconstructed "low-frequency (2 Hz-50 Hz) P- and S-waves propagating through the
vertical borehole arrays ... by the virtual source method", in a producing oil
field in Romania during water-flooding. Two statements matter for us.

1. **The illumination was from the surface, and they could tell.** "The obtained
   interferograms clearly indicate an origin of the ambient seismic energy from
   above the arrays, thus suggesting surface activities as sources."
2. **Very little data was needed.** "Ambient noise from time periods as short as
   **30 seconds** is sufficient to obtain robust interferograms," and the
   interferograms were of "consistency and high quality ... throughout the entire
   observation period."

Point 2 is diagnostic for us and I had the inference backwards for most of this
thread. We treated non-detection as an SNR problem and answered it with more
data: 6 days, then a coherent 96-hour stack. The results went the *wrong* way
(min p = 0.1345 per day; 96 h coherent stack p = 0.9184). Under adequate
illumination 30 s suffices, and stacking that fails to converge is the signature
of an **absent** signal rather than a weak one.

**Lellouch, Yuan, Spica, Biondi & Ellsworth (2019).** *Seismic velocity
estimation using passive downhole distributed acoustic sensing records: examples
from the San Andreas Fault Observatory at Depth.* JGR: Solid Earth.
https://doi.org/10.1029/2019JB017533

The target paper. Two details recovered from the text that we had not been using:

- **His inference of surface sources is a measurement, not an assumption.** After
  5-20 Hz bandpass, "the strong signals of these downgoing waves compared to the
  upgoing ones suggest that the dominant ambient field sources originate at the
  surface", citing Zhou & Paulssen (2017). The **downgoing/upgoing asymmetry is
  therefore a prerequisite that can be tested independently of any velocity
  scan** — implemented as `ambient_directional_asymmetry.py`.
- **Acquisition:** June 2017, **OptaSense ODH3.1** on the longer SAFOD main-hole
  fibre to 864 m, analysis limited to 800 m by a loop failure at the fibre end;
  fibre under ~1 N tension in a 0.9 mm steel tube cemented between casing
  strings. Our 2024-25 files are `SAFOD-Zumberge-16mGL-500HzFs`, i.e. a
  **different interrogator** at 16 m gauge length — confirming the acquisition
  difference independently of any of our own measurements.

---

## 3. Methods already tried in this tree, with provenance

All of these separate energy **by velocity**, which is why
`AMBIENT_LOWK_MECHANISM.md` predicts they fail against a fixed-k contaminant.
Listed so the review is not mistaken for a list of untried options.

| family | provenance | our script | outcome |
|---|---|---|---|
| fixed-fan F-K, brick wall | standard exploration practice; Duncan & Beresford (1994) as cited by Isken et al. | `ambient_fk_transfer_test.fk_filter`, `ambient_signed_fk_v2.py` | apparent ridge, fails the pre-filter channel-scramble gate |
| fixed-fan F-K, raised-cosine taper | Gibbs/ringing mitigation | `ambient_fk_taper_test.py` | 229× less ringing; **no recovery**. (The stronger wording "no discriminating power gained" is **withdrawn** — the real and synthetic paths were not matched, and the corrected ratio is 1.51–2.01, not 0.63. See the note at `ambient_fk_taper_test.py:17`.) |
| linear Radon / tau-p slowness mute | Radon-constrained VSP separation | `ambient_radon_slant_stack.py` | pedestal survives every mute to 8 km/s |
| rank-k coherent subspace projection | randomised range finder | `--svd-rank` in `ambient_lellouch2019_exact_stack.py` | pedestal not suppressed |
| phase cross-correlation | Schimmel (1999), nu=2 | `ambient_pcc_pws.py` | amplitude-blind, does not remove a phase-coherent pedestal |
| phase-weighted stacking | PWS / tf-PWS / ts-PWS | `ambient_pcc_pws.py` | no recovery |
| offset-axis median flat-event removal | Rao & Wang (2016)-style f-k + median combination for VSP | `ambient_flat_event_removal.py` | no recovery |
| running-absolute-mean temporal normalisation | Bensen et al. (2007) | `ambient_transfer_test.preprocess` | standard, retained |

---

## 4. Gaps: what the literature offers that we have not tried

Ordered by expected value **as of when this section was written**. Items 1 and 3
have since been **closed** — both are marked below, and neither recovered the
arrival. Do not cite this section as a list of open options.

1. ~~**Adaptive F-K (AFK/NAFK), applied AFTER a spatial static-response
   removal.**~~ **CLOSED 2026-08-19 — run, negative.** Five configurations, none
   passing the gates; cleanest form (common mode + rank-2 + AFK) reached pedestal
   −0.036 and p = 0.9392, and the form with no prior static removal produced a
   *spurious* p = 0.0060 by amplifying the pedestal 6.6-fold. Full table and the
   generalisable null-design lesson in §1.
2. **Spatial calibration of the static per-channel response.** Not from this
   review — it is already this project's own proposal in
   `C2_PERMEABILITY_FOLLOWUP.md`, which records 5x headroom from calibrating the
   static per-channel amplitude response. §1 of `AMBIENT_LOWK_MECHANISM.md` is
   independent motivation: the contaminant is a fixed *spatial* pattern, so the
   removal must be spatial. This is a prerequisite for item 1, not an
   alternative to it.
3. ~~**Test the illumination prerequisite before any more filtering.**~~
   **CLOSED 2026-08-19 — run, and it is the BINDING CONSTRAINT.** The
   prerequisite was tested (`ambient_directional_asymmetry.py`,
   `interrogator_and_illumination_v2.py`, `illumination_window_scan.py`) and it
   **fails**, with a working positive control:

   | arm | \|A\| in the 2500-4000 m/s fan | p |
   |---|---:|---:|
   | Lellouch 2017 pre-event (positive control) | **0.348** | **0.0050** |
   | 2024-25, matched spatial rank | 0.040 | 0.7307 |

   The 2017 arm is significant at spatial ranks 0, 1 and 2; the 2024-25 arm at
   **no** rank tested (0, 1, 2, 4, 8). An archive-wide scan then sampled **240
   windows** from 2024-05-21 to 2025-05-06 under a decision rule fixed before the
   run: **11 windows reach p < 0.05 against 12.0 expected by chance** (Binomial
   95th percentile 18). No illuminated window exists anywhere in the sampled
   archive.

   Lellouch's downgoing P is an inference *from* that asymmetry, so with the
   asymmetry absent there is no arrival of that kind present to be recovered, and
   no filter can create a propagation direction that is not in the data. Items 1
   and 2 remain worth doing for the *pedestal*; they cannot produce the arrival.
   Note the stated limits: 30 s per window, ~5 s of 2017 data from two records,
   the frozen 2500-4000 m/s fan, |A| is along-fibre rather than exactly vertical,
   and rank-k removal is conservative against a genuinely low-rank plane wave.
4. **Short-window interferograms.** Behm (2016) got robust interferograms from
   30 s. We have only ever stacked upward. A sweep over window length from 30 s
   to 24 h would show whether anything appears and then *degrades* with
   stacking, which would point at non-stationary illumination rather than
   absence. **Still untried, but largely answered from another direction:**
   `ambient_stack_convergence.py` fits detectability against the number of
   stacked hourly chunks and finds **N^+0.042** (baseline) and **N^+0.019**
   (common-mode removed), where a real arrival requires **N^+0.50**. Stacking 24×
   more data should have improved detectability 4.9-fold; it improved it by about
   **13 %** (0.858 at 1 chunk to 0.971 at 24; 24^0.042 = 1.14). *Corrected
   2026-08-20 — this read "about 4 %", which was the exponent 0.042 misread as a
   percentage.* Detectability crosses 1.00 once, at 16 chunks (**1.014,
   p = 0.0170**), and falls back to 0.971 at 24; that excursion peaks at
   5,950 m/s with a pedestal diagnostic of +0.983, i.e. it is the scan-ceiling
   artefact, not a detection. *Corrected 2026-08-20 — this read "detectability
   never reaches the 1.00 a detection needs", which its own product
   (`ambient_stack_convergence.txt`) contradicts.* The
   illumination scan in item 3 also sampled 30 s windows directly and found
   nothing, so the "appears, then degrades with stacking" scenario has no
   support. Low priority now.

---

## 4a. What controls illumination, and why the two datasets differ

This is the question the review exists to answer, so it is stated explicitly.

**The geometry is unforgiving.** The virtual-source method recovers the Green's
function between two receivers only if noise sources occupy the *stationary-phase
zone* for that path. For two channels on a **vertical** line the direct P path
between them is vertical, so the stationary-phase sources lie on the extension of
that line: a narrow cone **directly above the wellhead**, radiating steeply down.
Noise arriving at oblique incidence contributes nothing to that path and simply
adds to the background. This is why the DAS ambient-noise literature is dominated
by *surface-wave* dispersion studies (Shao et al. 2023; the Ebao Basin and
Sanriku studies in §1's search results) — for surface waves, sources anywhere
along the fibre azimuth are usable. A downgoing body wave in a borehole is the
hard case.

**What fills that cone is cultural and local:** vehicles and machinery on the
wellhead pad, drilling, pumps, work on the access road, wind coupling into the
casing. Not ocean microseism, not distant earthquakes — those do not arrive
vertically.

**Consequently the medium being identical is not sufficient.** Same fibre, same
borehole, same rock guarantees the same *Green's function*; it does not supply the
*sources* needed to reconstruct it. The two recordings differ in source
environment, not geology:

- **June 2017** (Lellouch et al.): an OptaSense ODH3.1 was attached to the
  main-hole fibre that month.
- **2024-25** (this archive): unattended continuous recording.
- **Behm (2016)**, for contrast, succeeded in a **producing oil field during
  water-flooding** — about the noisiest wellhead environment obtainable. This one
  IS documented in the paper and is the only source-environment datum here that
  is not inference.

**Status: SPECULATION, and weaker than the rest of this document. Do not cite it.**
An earlier draft of this section asserted that 2017 was recorded "during a manned
field deployment" whose crew and vehicles supplied the illumination. There is no
evidence for that beyond the fact that an interrogator was installed sometime in
June 2017. Specifically, there is **no** site log, **no** operational record, and
**no** statement in Lellouch et al. (2019) attributing the noise to any source —
the paper infers surface origin from the downgoing/upgoing asymmetry alone, which
is the same inference we make, not independent corroboration of a cause.

**A gap in the comparison, stated plainly.** Our 2017 asymmetry measurement uses
the pre-event windows of the two earthquake records `M1p33` and `M2p46` (~5 s
total), because that is the only raw 2017 noise in the release. Lellouch's Figure
7c correlograms are, per the release README, "a stack of 7 different one-day
correlations" — raw data we do not have. So the arm we measure is **not the same
acquisition** as the one that produced his figure. This does not affect the
measured contrast (asymmetry significant at p = 0.005 in the 2017 records,
p = 0.73 in 2024-25, matched operations and ranks), but it does mean we cannot
attribute the difference to any specific cause, temporal or instrumental.

What would settle it: SAFOD/USGS operational records for June 2017 and for the
2024-25 window, or the raw ambient data behind his Figure 7c. Neither is in hand.

It is testable on **our** side, and `illumination_window_scan.py` has now done so.
The reasoning was that cultural noise has a **diurnal signature**, so any
illuminated windows in the 2024-25 archive should cluster in working hours.

**Result (2026-08-19): there are no illuminated windows to cluster.** Of 240
windows sampled across the full archive, 11 reach p < 0.05 against 12.0 expected
by chance (Binomial 95th percentile 18); |A| has median 0.0296, 90th percentile
0.0753, maximum 0.2385. The hour-of-day test was declared *supporting, not
gating*, before the run, and it is moot — with the count at chance level there is
no positive signal whose timing could be examined.

So the cultural-source *explanation* is neither confirmed nor refuted by this: it
remains speculation, exactly as the paragraph above says. What is now measured is
the **absence** on our side, which does not require knowing the cause. Keep those
two separate. The absence is a result; the cause is not.

## 5. Honest summary of the goal

The goal was to recover the signal, review the literature for how others do
ambient-noise CC, and focus on robust F-K implementations.

Updated 2026-08-19, after the two gaps this section left open were closed.

- **Literature reviewed:** §1-§3. The material finding was that a whole F-K family
  — adaptive, data-driven, sliding-window — existed and had never been tried here,
  so the criticism that the F-K work was not pushed hard enough was **correct as
  to breadth**. That gap is now **closed**: adaptive f-k was implemented and run
  in five configurations (§1) and none passed the gates.
- **Signal not recovered.** No method in or out of this review has produced an
  arrival that clears the controls. Three quantitative statements of that
  negative, none depending on a multiplicity correction:
  - six days give min p = 0.1345; Fisher combination p = 0.524; the coherent
    four-day 96 h stack gives p = 0.9184;
  - **0 of 181 velocities clear the per-velocity null** — at 3,200 m/s, 2.752
    against a threshold of 2.811;
  - detectability grows as N^+0.042 / N^+0.019 rather than the N^+0.50 that
    stacking a real arrival requires.
  The one selected day-pair at p = 0.039 is a **selection artefact**, not a
  result: it was chosen because it scored lowest, and it peaks at the scan ceiling
  with causal/acausal 0.97.
- **What changed, and in what order.** The reason for the failure is no longer
  unknown, and it is **two** things with a definite precedence:
  1. **Illumination is the binding constraint** (§4.3). The 2024-25 field carries
     no downgoing/upgoing asymmetry in the fan (|A| = 0.040, p = 0.7307) where
     Lellouch's 2017 records do (|A| = 0.348, p = 0.0050), and no illuminated
     window was found in 240 sampled across the archive. This is about the
     **recording**, it is upstream of all processing, and it rules the route out.
  2. **A static fixed-k spatial pattern is the contaminant**
     (`AMBIENT_LOWK_MECHANISM.md`), which explains why eight velocity-domain
     methods failed identically and why the pedestal could not be filtered out.
     This is about the **processing** failures, not the absence of the arrival.
  Do not cite (2) alone as the reason Figure 7c does not reproduce.
- **What did not change:** the F-K fan verdict in `Ambient_FK_QC_workflow.ipynb`
  **still stands** — the 2.5-4.5 km/s fan produces an apparent ridge but fails the
  pre-filter channel-scramble gate, so it is not an independently recovered
  physical arrival. Nothing in this review or in the adaptive-f-k work reopens it;
  §1 row 1 is a second, independent instance of the same failure mode.
- **The transferable methodological result.** A pre-gather operator sits *outside*
  a receiver-order null, because that null permutes the finished gather. AFK with
  no prior static removal exploited this to reach p = 0.0060 while carrying the
  worst pedestal of the five configurations. Any pre-correlation operator in this
  tree owes an **input-level** null, not just a gather-level one.
- **Ruled out as explanations:** gauge length (0.1-1.1 % at 3,200 m/s), and the
  interrogator inside the analysed aperture — over channels 23-708 the dominant
  spatial pattern is **not** stable across days (|corr| 0.0188 against a control
  95th of 0.1282), so no instrumental fingerprint is supported there. *Corrected
  2026-08-20 — this previously attributed |corr| 0.8426 vs 0.3889 to "channels
  0-22". That measurement is v1's run over channels **0-699** with the lead-in
  included (`interrogator_blame_test.txt`); no channels-0-22-only correlation was
  ever computed. See `AMBIENT_LOWK_MECHANISM.md` §5.*

## References

Baran, I., Stewart, M. P., Kampes, B. M., Perski, Z., & Lilly, P. (2003). A
modification to the Goldstein radar interferogram filter. *IEEE Transactions on
Geoscience and Remote Sensing*. https://doi.org/10.1109/tgrs.2003.817212

Behm, M. (2016). Feasibility of borehole ambient noise interferometry for
permanent reservoir monitoring. *Geophysical Prospecting*, 65(2), 563-580.
https://doi.org/10.1111/1365-2478.12424

Isken, M. P., Vasyura-Bathke, H., & Dahm, T. (2022). De-noising distributed
acoustic sensing data using an adaptive frequency-wavenumber filter.
*Geophysical Journal International*, 231(2), 944-949.
https://doi.org/10.1093/gji/ggac229

Lellouch, A., Yuan, S., Spica, Z., Biondi, B., & Ellsworth, W. L. (2019).
Seismic velocity estimation using passive downhole distributed acoustic sensing
records: examples from the San Andreas Fault Observatory at Depth. *Journal of
Geophysical Research: Solid Earth*. https://doi.org/10.1029/2019JB017533

Note on citation hygiene: Duncan & Beresford (1994), Goldstein & Werner (1998),
Zhou & Paulssen (2017), Schimmel (1999), Bensen et al. (2007) and Rao & Wang
(2016) appear above **as cited by** the papers retrieved here or by this
project's own existing documentation. They were not independently retrieved in
this review and should be checked against
`notebooks/LITERATURE.md` and `/home/groups/ettore88/nberrios/planpapers/` before
being cited in the manuscript.
