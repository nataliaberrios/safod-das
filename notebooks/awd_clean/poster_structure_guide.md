# Poster structure guide

## Recommended poster title

**Installation-dependent coherent modes and sensitivity limits of repeated-source DAS at SAFOD**

This title is intentionally accurate. It does not imply that the poster has already resolved a formation (V_P) profile, permeability, or natural Parkfield (dv/v).

## Overall layout

Use a landscape poster with three vertical columns and a strong left-to-right argument:

**Experiment → Nano result → Installation comparison and limits**

Keep the PGSI calibration in the third column as an interpretation constraint, not as the headline result.

---

# Header band

## Title, authors, affiliations

Under the title, add a one-sentence central result:

> Repeated AWD excitation produces a statistically significant, repeatable approximately 3 km s⁻¹ apparent mode on the cemented Nano fiber; the wireline Deep installation emphasizes a distinct slower mode.

Add a small QR code or link to the v8 notebook and repository if appropriate.

---

# Column 1 — Experiment and question

## 1. Motivation

Explain why repeated-source DAS at SAFOD is useful:

- DAS installation and coupling can select different observable modes.
- Repeated AWD drops allow repeatability to be measured before stacking.
- A calibrated sensitivity test is needed before making monitoring claims.

Do not lead with Earth tides or permeability. Those are context and future applications, not established results here.

## 2. Scientific questions

Use three or four short questions:

1. Is there a coherent, spatially ordered Nano mode?
2. Does that mode repeat at individual-drop and burst scales?
3. Do Nano and Deep emphasize different modes because of installation?
4. What sensitivity and propagation aperture does the current experiment support?

## 3. Experiment and geometry

Use Core Figure 1.

Include:

- AWD source location and approximate 15 m lateral offset.
- Cemented Nano installation.
- Wireline Deep installation.
- 49 bursts and 988 drops.
- A clear statement that the schematic is not a surveyed depth map.

## 4. Processing overview

Use a short flow diagram:

**catalog/manifest → burst grouping → common-drop stacks → signed slowness → repeatability → injection–recovery**

Keep this section methodological and brief.

---

# Column 2 — Main Nano result

## 5. Nano coherent mode

Use Core Figure 2 as the largest figure on the poster.

## Question

Does the Nano data contain a spatially ordered coherent mode, and what is its apparent speed and direction?

## Main result

State:

- Positive signed slowness.
- Approximately 2.95–2.99 km s⁻¹ apparent speed.
- Coherent frequency-dependent ridge.
- The mode is phase-neutral: direct P versus guided/coupled energy remains unresolved.

Avoid calling this “formation (V_P)” in the main result box.

## 6. Repeatability

Use Core Figure 3.

## Question

Does the mode recur before full-experiment stacking?

## Main result

State:

- 988 individual drops analyzed.
- 49 qualifying bursts.
- Signal NCC exceeds the noise comparison across all bursts.
- Across-burst and within-burst convergence are both strong.

Add a small caveat:

> These metrics quantify repeatability of a fixed apparent-moveout observable; they are not 988 independent velocity estimates.

## 7. Statistical validation

Use a compact box or one small supporting panel showing the spatial-permutation result and pre-source control.

The purpose is to answer:

> Is the Nano ridge associated with the actual channel ordering rather than arbitrary spatial structure?

Keep detailed null distributions in supplementary material if space is limited.

---

# Column 3 — Installation contrast, sensitivity, and calibration

## 8. Installation-dependent mode selection

Use Core Figure 4.

## Main result

State:

- Nano: approximately 3 km s⁻¹ mode at higher frequencies.
- Deep: approximately 1.4–1.56 km s⁻¹ low-frequency mode.
- The installations emphasize different coherent modes.

Use careful language:

> The Deep mode is consistent with a tube wave, but tube-wave identity and conversion depth remain preliminary.

Do not compare raw Nano and Deep amplitudes unless instrument response and coupling have been calibrated.

## 9. Frequency-dependent aperture and sensitivity

Use Core Figure 5.

Split the interpretation into two statements:

- The Nano mode remains observable through approximately the 460 m analyzed aperture in the strongest band.
- Blind injection–recovery reaches 1% as the smallest tested two-direction 95%-correct-sign level for the selected apparent-moveout estimator.

Add a bold limitation:

> This is not a demonstrated detection threshold for natural Parkfield-scale (dv/v).

## 10. PGSI timing and moveout calibration

Use Core Figure 6 as a smaller supporting figure.

## Main result

State:

- One canonical timing convention was applied to all six explosion records.
- Five of six shots pass the fixed moveout QC.
- Shot 1040 is the only PGSI fit compatible with the Nano 2.95–2.99 km s⁻¹ interval.
- The comparison with the Nano depth axis is conditional because source, sensor, coupling, and absolute registration differ.

This section should be labeled **Calibration constraint**, not **(V_P) result**.

---

# Bottom band — Conclusions and next steps

## 11. Conclusions

Use three boxes:

### Established

- A statistically significant positive-slowness Nano mode exists.
- The mode is repeatable within drops, bursts, and the full experiment.
- Nano and Deep installations emphasize different modes.
- The current estimator has a measured observable-level sensitivity floor.

### Unresolved

- Direct formation P versus guided/coupled Nano mode.
- Physical identity of the Deep slow mode.
- Absolute PGSI-to-Nano depth registration.
- Whether the faster PGSI arrivals represent a different phase or coupling regime.

### Not claimed

- A formation (V_P(z)) inversion.
- A permeability or fracture-flow inference.
- Detection of Earth tides or natural Parkfield-scale (dv/v).

## 12. Next steps

Choose only two or three:

1. Investigate why four usable PGSI shots are faster than shot 1040.
2. Build a formal forward comparison of Nano and PGSI moveout under candidate phase models.
3. Test structural or tube-wave amplitude/reflection observables without requiring a full (V_P) inversion.

## 13. Data and code availability

List:

- SAFOD June 2026 AWD-DAS data.
- The v8 dashboard notebook.
- The figure-generation scripts.
- The PGSI reference geometry and timing files.

Include a QR code if allowed.

---

# What the poster should not do

- Do not make the PGSI figure the central result.
- Do not call the Nano ridge a formation P wave in the title.
- Do not call the 1% injection level an Earth-tide detection capability.
- Do not compare raw Nano and Deep amplitudes as if the installations were calibrated.
- Do not include every diagnostic figure.
- Do not leave the audience to infer the conclusion from plots; put the established, unresolved, and not-claimed statements in text.

## Suggested visual hierarchy

- Largest: Core Figure 2, Nano coherent mode.
- Second largest: Core Figure 3, repeatability.
- Medium: Core Figures 4 and 5.
- Smaller: Core Figures 1 and 6.
- Text boxes: conclusions and limitations.

The poster should be understandable by reading only the title, central-result sentence, figure titles, and three conclusion boxes. Everything else should support that path.
