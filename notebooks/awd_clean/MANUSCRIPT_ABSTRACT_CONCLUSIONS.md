# Manuscript abstract and conclusions

Following the plan's §21 sequence. Every number below is verified against a
generated output; nothing new is introduced here. The Introduction is drafted in
[`MANUSCRIPT_INTRODUCTION.md`](MANUSCRIPT_INTRODUCTION.md).

---

## Abstract

Borehole distributed acoustic sensing is increasingly used for time-lapse
velocity monitoring, but the consequences of installation and coupling for
monitoring sensitivity have not been quantified against a controlled repeated
source. We recorded a 24-hour repeated accelerated weight-drop survey on two
borehole DAS installations in the same hole at SAFOD: a shallow cemented fibre
and a deeper wireline fibre that reverses, giving outbound and return branches.

The two installations recover different coherent modes. The cemented fibre
records a fast apparent mode near 2950 m s⁻¹, strongest at 30–60 Hz, with
resolved frequency-dependent slowness. The wireline fibre records a slow
borehole-guided mode near 1547 m s⁻¹, strongest at 15–30 Hz, whose trajectory is
selected independently on each branch to within 0.3% and which survives
channel-order permutation testing at the resolution floor of 499 permutations.
Individual impacts vary substantially — beam signal-to-noise spans 1.26 to
24.90 dB between the 16th and 84th percentiles — but stacking drops into bursts
yields highly repeatable observables, with burst-stack correlation of 0.976
against an independent reference.

Blind injection–recovery on held-out bursts, using an identical estimator for
both fibres, establishes smallest tested reliable-detection levels of 1.0% for
the cemented observable, 0.5% for the outbound wireline branch, and 1.0% for the
return branch. Although the outbound regression lever arm is 11.5 times longer
than the cemented one, poorer burst-to-burst timing repeatability — 1.68 ms
against 0.225 ms — limits the realised sensitivity gain to a factor of two.
Combining the two branches lowers the empirical noise floor by roughly 1.5 times
but does not lower the smallest tested detection level.

Installation therefore controls both which coherent mode is observable and the
sensitivity of the resulting time-lapse measurement, and the advantage of a
longer propagation path is branch-specific rather than installation-wide.

**Candidate result sentence for the body**

> Blind injection–recovery established smallest tested two-direction
> reliable-detection levels of 1.0% for the Nano apparent-moveout observable,
> 0.5% for the outbound Deep guided branch, and 1.0% for the return branch.
> Although the Deep outbound regression lever arm was 11.5 times longer than
> Nano's, poorer burst-to-burst timing repeatability limited the demonstrated
> sensitivity gain to a factor of two.

---

## 6. Conclusions

**1. Contrasting installations in one borehole recover distinct coherent modes.**
A cemented fibre and a wireline fibre, recording the same repeated source over
the same nominal 24 h survey, returned a fast apparent mode near 2950 m s⁻¹ at 30–60 Hz and
a slow guided mode near 1547 m s⁻¹ at 15–30 Hz respectively. Both are dispersive;
neither apparent speed should be read as a formation velocity. Installation is an
experimental variable, not a deployment detail.

**2. Burst stacking converts a variable impact source into a repeatable
observable.** The coherent mode is present in individual drops — median
signal-window correlation 0.889 against 0.288 in an equal-duration noise window,
with all 49 of 49 bursts favouring signal over noise (p = 1.78×10⁻¹⁵) — but
individual impacts vary by a factor of ~20 in beam power. Independent substacks
converge from median correlation 0.740 at one drop to 0.909 at eight, and full
burst stacks reach 0.976. Usable precision is bought by stacking, so stacking
design determines what is measurable.

**3. Blind injection–recovery establishes mode- and branch-specific monitoring
sensitivity.** At the tested levels: 1.0% for the cemented observable, 0.5% for
the outbound wireline branch, and 1.0% for the return branch, with empirical null
thresholds of 0.325%, 0.184% and 0.305%. The improved-sensitivity claim is
restricted to the best-observed outbound branch; the wireline installation as a
whole is not shown to outperform the cemented one.

**4. A longer propagation-time lever arm improves fractional sensitivity only in
proportion to timing repeatability.** Fractional velocity sensitivity is set by
the ratio of per-burst timing repeatability to lever arm in travel time. The
outbound branch has 11.5 times the lever arm but 7.4 times worse timing
repeatability, and 11.5/7.4 = 1.55 reproduces the observed ratio of null scatter
exactly. The realised improvement in reliable detection is therefore a factor of
two, not an order of magnitude. Experiments should specify both terms; aperture
or depth alone does not constrain monitoring performance.

**5. Combining spatially distinct branches improves precision without
necessarily improving detection.** The two branches' errors are nearly
uncorrelated (ρ = 0.121), and combining them reduces robust null scatter by ~1.5
times, close to the independent-averaging prediction. The smallest tested
reliable-detection level nonetheless remains 0.5% for every weighting. Precision
gains accrue continuously; detection levels move in discrete steps.

---

## Interpretation limits carried into the conclusions

The measured quantity throughout is a fractional change in the apparent
along-fibre speed of a selected guided mode. It is not a formation V_P or V_S
change, not a measure of stress, pore pressure, permeability, fracture
compliance, or tectonic strain, and it carries no depth resolution. The Niu SAFOD expected-response scale for a tidal velocity change, 5.76×10⁻⁵,
lies roughly two orders of magnitude below the best demonstrated detection level;
a direct fit to the per-burst estimates returns a clean null with a 95% upper
limit of 1.03×10⁻³, some 18× above that scale.

## Outstanding before submission

| Item | Status |
|---|---|
| Introduction | **drafted**, from papers read; citation status tabulated inline |
| Niu benchmark values (2.4×10⁻⁷ Pa⁻¹, 240 Pa) | **verified verbatim** in Niu et al. (2008) |
| De Fazio as a threshold | **retired**; results are scored against the Niu SAFOD expected-response scale, De Fazio kept as historical context at ~5×10⁻⁴–10⁻³ |
| Methods §3.1 preprocessing | **drafted** from `PREPROCESSING.md` |
| Nano gauge length / channel spacing / raw rate | **resolved** from the protobuf header: 16.459 m, 1.2661 m, 10 kHz |
| Tidal placement | supplement to carry both the Niu arithmetic and the measured upper limit, labelled a limit rather than a search |
| Interrogator **model numbers** | open — vendors established (Sintela, OptaSense); models need Ettore |
| Physical depth of the Deep turnaround (channel 1702) | open — needs Ettore |
| Which physical cable each dataset is | open — the project's own fibre-identity question |
| Two-model tide notebook (Model A / Model B) | not yet on Sherlock; Model A intended for main text, Model B for supplement |
| Prose vs internal-record separation | pending — these files still carry filenames, defect history and reviewer-facing asides that belong in reproducibility material, not the manuscript |
