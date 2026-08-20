# Literature review: robust ambient-noise CC, and robust F-K filtering in particular

Written 2026-08-19 for the Figure 7c recovery goal. Every paper below was
retrieved and read through the Scite full-text index or the publisher; nothing
here is cited from memory. Where a paper's method has already been tried in this
tree, the script is named. Where it has not, that is stated as a gap.

Companion documents: `AMBIENT_LOWK_MECHANISM.md` (what the contaminant is),
`Ambient_FK_QC_workflow.ipynb` (why the fixed 2.5-4.5 km/s fan was rejected).

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

and only in that order. This is the single highest-value untried method in this
review. Their open-source implementation exists, which removes most of the cost.

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
| fixed-fan F-K, raised-cosine taper | Gibbs/ringing mitigation | `ambient_fk_taper_test.py` | no discriminating power gained |
| linear Radon / tau-p slowness mute | Radon-constrained VSP separation | `ambient_radon_slant_stack.py` | pedestal survives every mute to 8 km/s |
| rank-k coherent subspace projection | randomised range finder | `--svd-rank` in `ambient_lellouch2019_exact_stack.py` | pedestal not suppressed |
| phase cross-correlation | Schimmel (1999), nu=2 | `ambient_pcc_pws.py` | amplitude-blind, does not remove a phase-coherent pedestal |
| phase-weighted stacking | PWS / tf-PWS / ts-PWS | `ambient_pcc_pws.py` | no recovery |
| offset-axis median flat-event removal | Rao & Wang (2016)-style f-k + median combination for VSP | `ambient_flat_event_removal.py` | no recovery |
| running-absolute-mean temporal normalisation | Bensen et al. (2007) | `ambient_transfer_test.preprocess` | standard, retained |

---

## 4. Gaps: what the literature offers that we have not tried

Ordered by expected value.

1. **Adaptive F-K (AFK/NAFK), applied AFTER a spatial static-response removal.**
   Isken et al. (2022). See §1. Highest value, open-source implementation exists,
   and it carries a clear falsifiable prediction in both orderings.
2. **Spatial calibration of the static per-channel response.** Not from this
   review — it is already this project's own proposal in
   `C2_PERMEABILITY_FOLLOWUP.md`, which records 5x headroom from calibrating the
   static per-channel amplitude response. §1 of `AMBIENT_LOWK_MECHANISM.md` is
   independent motivation: the contaminant is a fixed *spatial* pattern, so the
   removal must be spatial. This is a prerequisite for item 1, not an
   alternative to it.
3. **Test the illumination prerequisite before any more filtering.**
   `ambient_directional_asymmetry.py`, from Lellouch's own diagnostic (§2). If
   the field carries no net directional preference in the body-wave fan, no
   filter can create one and items 1-2 cannot succeed either. This is cheap and
   should gate the rest.
4. **Short-window interferograms.** Behm (2016) got robust interferograms from
   30 s. We have only ever stacked upward. A sweep over window length from 30 s
   to 24 h would show whether anything appears and then *degrades* with
   stacking, which would point at non-stationary illumination rather than
   absence. Untried, cheap.

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

It is testable on **our** side, and `illumination_window_scan.py` does so:
cultural noise has a **diurnal signature**, so if any illuminated windows exist in
the 2024-25 archive they should cluster in working hours. A diurnal pattern in
|A| would support the cultural-source explanation even if no window reaches
significance.

## 5. Honest summary of the goal

The goal was to recover the signal, review the literature for how others do
ambient-noise CC, and focus on robust F-K implementations.

- **Literature reviewed:** §1-§3. The material finding is that a whole F-K family
  — adaptive, data-driven, sliding-window — exists and was never tried here, so
  the criticism that the F-K work was not pushed hard enough is **correct as to
  breadth**, even though the family we did try was tested thoroughly.
- **Signal not recovered.** No method in or out of this review has produced an
  arrival that clears the controls.
- **What changed:** the reason for the failure is no longer unknown. The
  contaminant is a static fixed-k spatial pattern (`AMBIENT_LOWK_MECHANISM.md`),
  which explains why eight velocity-domain methods failed identically, and it
  redirects the work to a spatial remedy plus an adaptive rather than fixed-fan
  filter — with an illumination prerequisite (§4.3) that should be checked first
  because it can rule the whole route out.

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
