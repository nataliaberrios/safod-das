# SAFOD DAS — literature organised by project direction

The collection lives in `/home/groups/ettore88/nberrios/planpapers/`, sorted into
folders that mirror the sections below and renamed to `author_year_journal_topic.pdf`.
`MANIFEST.csv` there records every rename, so the whole reorganisation reverses
with one command (see `notebooks/organise_papers.py`).

| folder | papers |
|---|---|
| `01_repeaters` | 12 |
| `02_tides` | 7 |
| `03_ambient_noise` | 4 |
| `04_methods_das` | 7 |
| `05_site_safod_parkfield` | 14 |
| `06_methods_borehole` | 1 |
| `_duplicates` | 1 |

Organised 2026-08-05.

**On the organisation.** Three of the four directions you named — tides, repeating
earthquakes, ambient-noise correlation — are *scientific targets*: each could carry
a claim on its own. Active source is, as you said, a method rather than a result,
so it sits with the other methods. I've split the methods into two groups because
they answer different questions: **§4 is "can this instrument measure it"** and
**§5 is "what is already known about this hole."** Several papers serve more than
one direction; those are cross-referenced rather than duplicated.

---

## 1. Repeating earthquakes → fault slip at depth

*The current active direction. Target: aseismic creep rate from recurrence.*

### The method lineage — read in this order

| paper | file | why |
|---|---|---|
| **Poupinet, Ellsworth & Fréchet 1984**, JGR 89:5719 | `01_repeaters/`<br>`poupinet_ellsworth_frechet_1984_jgr_earthquake-doublets-calaveras.pdf` | Founding paper. Earthquake doublets to monitor velocity, Calaveras Fault. Cross-spectral phase for sub-sample delay — still the best delay estimator. |
| **Nadeau, Foxall & McEvilly 1995**, Science 267:503 | `01_repeaters/`<br>`nadeau_1995_science_clustering-periodic-recurrence-parkfield.pdf` | Clustering and *periodic* recurrence at Parkfield. Establishes that these things repeat on a clock. |
| **Nadeau & Johnson 1998**, BSSA 88(3):790 | `01_repeaters/`<br>`nadeau_johnson_1998_bssa_parkfield-VI-moment-release-scaling.pdf` | "Parkfield VI: Moment Release Rates and Source Parameters." **The slip–moment scaling relation.** Note: calibrated *against* geodetic creep, so using it to measure creep at Parkfield is partly circular. |
| **Nadeau & McEvilly 1999**, Science 285:718 | `01_repeaters/`<br>`nadeau_mcevilly_1999_science_fault-slip-rates-at-depth.pdf` | **"Fault Slip Rates at Depth from Recurrence Intervals of Repeating Microearthquakes."** The paper this project's current direction reproduces. Four pages. Start here. |

### Why the scaling is contested — read after the above

| paper | file | why |
|---|---|---|
| **Beeler & Hickman 2001**, JGR 106:30701 | `01_repeaters/`<br>`beeler_hickman_2001_jgr_stress-drop-strength-recovery.pdf` | Stress drop vs lab-inferred interseismic strength recovery. Magnitude-dependent aseismic slip as the explanation for N&J's scaling. |
| **Chen & Lapusta 2009**, JGR 114:B01311 | `01_repeaters/`<br>`chen_lapusta_2009_jgr_scaling-repeaters-rate-and-state.pdf` | Rate-and-state simulation reproducing the scaling: a velocity-weakening patch inside a velocity-strengthening region. |
| **Abercrombie 2014**, GRL | `01_repeaters/`<br>`abercrombie_2014_grl_stress-drops-repeating-earthquakes.pdf` | Measured stress drops of Parkfield repeaters. **Contradicts** N&J's implied stress-drop-rises-as-moment-falls trend. |
| **Gao, Kao & Wang 2021**, GRL | `01_repeaters/`<br>`gao_2021_grl_misconception-waveform-similarity.pdf` | "Misconception of Waveform Similarity in the Identification of Repeating Earthquakes." Cautionary — high CC does not imply co-located rupture. |

### Reviews and applications

- **Uchida & Bürgmann 2019**, Annu. Rev. 47:305 — `01_repeaters/`<br>`uchida_burgmann_2019_annurev_repeating-earthquakes-review.pdf`. The synthesis. §"Source characteristics and recurrence" contains every caveat above in one place.
- **Uchida 2019**, PEPS 6:40 — `01_repeaters/`<br>`uchida_2019_peps_detection-of-repeating-earthquakes.pdf`. Detection methods; Table 1 compiles windows and bands per study.
- **Lengliné & Marsan 2009**, JGR 114:B10303 — coseismic/postseismic stress change at Parkfield from repeaters.
- **Ide, Beroza, Shelly & Uchide 2007**, Nature 447 — "A scaling law for slow earthquakes." Adjacent: the aseismic end of the spectrum.

*Also relevant:* Rubinstein & Beroza 2005 (§5), which uses the SAFOD target repeaters specifically.

---

## 2. Tides → velocity modulation and poroelastic response

*A genuinely distinct target, and the one with the cleanest forcing: tides are
predictable to high precision, so the null is exactly computable.*

| paper | file | why |
|---|---|---|
| **De Fazio, Aki & Alba 1973**, JGR 78:1319 | `02_tides/`<br>`defazio_aki_alba_1973_jgr_solid-earth-tide-velocity-change.pdf` | **The founding observation** — solid Earth tide changing in-situ seismic velocity. Aki is a co-author. |
| **Takano et al. 2014**, GRL | `02_tides/`<br>`takano_2014_grl_velocity-changes-earth-tide-ambient-noise.pdf` | **The modern template.** Velocity change from the Earth tide, detected by ambient-noise correlation. This is the paper your tides direction would follow. |
| **Berger 1975**, JGR 80:274 | `02_tides/`<br>`berger_1975_jgr_thermoelastic-strains-and-tilts.pdf` | Thermoelastic strains and tilts. Short. The other periodic forcing that must be separated from tides. |
| **Ben-Zion & Leary 1986**, BSSA 76:1447 | `02_tides/`<br>`benzion_leary_1986_bssa_thermoelastic-strain-halfspace.pdf` | Thermoelastic strain in a half-space under unconsolidated cover. Gives the *depth dependence* — essential for telling a tidal signal from a thermal one. |
| **Métivier & Conrad 2008**, JGR 113:B11405 | `02_tides/`<br>`metivier_conrad_2008_jgr_body-tides-heterogeneous-earth.pdf` | Body tides of a heterogeneous, aspherical Earth. The forward model for what strain to expect. |
| **Rojstaczer 1989**, JGR 94:12403 | `02_tides/`<br>`rojstaczer_1989_jgr_well-water-levels-response-earth-tides.pdf` | Water levels in wells responding to Earth tides. The poroelastic transfer function. |
| **HESS 2022**, 26:4301 | `02_tides/`<br>`hess_2022_groundwater-response-earth-atmospheric-tides.pdf` | In-situ hydro-geomechanical properties from groundwater response to Earth *and atmospheric* tides. Atmospheric loading is the confound. |

**Why this direction is attractive:** the forcing is known exactly, so unlike a
seasonal signal you can predict the phase and amplitude in advance and test against
it. The M2 tide is ~12.42 h, far from any diurnal instrumental artifact.

**Why it's hard:** tidal strain is ~10⁻⁸, and the implied dv/v is correspondingly
tiny. Takano et al. is the place to check what amplitude is actually detectable.

---

## 3. Ambient-noise cross-correlation → continuous monitoring

*The direction with the most existing infrastructure here — `stack_daily.py` and
`cc_tools.py` already run over 413 covered days.*

| paper | file | why |
|---|---|---|
| **Wapenaar et al. 2010**, Geophysics 75:75A195 | `03_ambient_noise/`<br>`wapenaar_2010_geophysics_interferometry-tutorial-part1.pdf` | Tutorial on seismic interferometry, Part 1. Read first if the Green's-function-from-noise idea isn't yet intuitive. |
| **Shapiro & Campillo 2004**, GRL 31:L07614 | `03_ambient_noise/`<br>`shapiro_campillo_2004_grl_broadband-rayleigh-from-noise.pdf` | Foundational: broadband Rayleigh waves emerging from noise correlation. |
| **Brenguier et al. 2008**, Science 321:1478 | `03_ambient_noise/`<br>`brenguier_2008_science_postseismic-relaxation-parkfield.pdf` | **Postseismic relaxation at Parkfield from continuous seismological observations.** dv/v monitoring *at your site*, on the 2004 M6. The direct precedent. |
| **Li & Ben-Zion 2023**, JGR 128 | `03_ambient_noise/`<br>`li_benzion_2023_jgr_daily-seasonal-shallow-velocity.pdf` | Daily and seasonal shallow velocity variations, southern California. 4% seasonal, 10% daily, peak sensitivity ~17 m, driven by thermoelastic strain and soil moisture. |

*Also here:* Takano et al. 2014 (§2) uses ambient noise as its tool.

**Missing and important — see §7:** Sens-Schönfelder & Wegler 2006, which
introduced the stretching estimator and the seasonal application. It is the single
most-cited method paper for this direction and you don't have it.

---

## 4. DAS as an instrument — cross-cutting method

*Applies to all three targets. The question these answer is "can fibre measure
this," not "what does it mean."*

| paper | file | why |
|---|---|---|
| **Zhan 2020**, SRL | `04_methods_das/`<br>`zhan_2020_srl_das-fiber-optic-seismic-antennas.pdf` | "DAS Turns Fiber-Optic Cables into Sensitive Seismic Antennas." The accessible overview — good to hand to anyone unfamiliar. |
| **Lindsey, Rademacher & Ajo-Franklin 2020**, JGR 125 | `04_methods_das/`<br>`lindsey_2020_jgr_broadband-instrument-response-das.pdf` | Broadband instrument response, calibrated against a broadband seismometer. Amplitude and **phase** response — phase flatness is what lets you trust delay measurements. |
| **Ichinose et al. 2022**, JGR 127 | `04_methods_das/`<br>`ichinose_2022_jgr_das-vs-array-derived-strain-rate.pdf` | DAS strain rate vs array-derived strain rate. Where DAS agrees with conventional sensors and where it stops (~2–5 Hz). |
| **Lellouch, Lindsey, Ellsworth & Biondi 2020**, SRL 91:3256 | `04_methods_das/`<br>`lellouch_2020_srl_das-vs-geophones-forge.pdf` | DAS vs **collocated** geophones at FORGE. Magnitude of completeness −1.4 (DAS) vs −1.7 (geophones): DAS loses by 0.3 units *after* array processing. The honest benchmark. |
| **Madsen, Tøndel & Kvam 2016**, TLE 35:610 | `04_methods_das/`<br>`madsen_2016_tle_data-driven-depth-calibration-das.pdf` | Data-driven depth calibration for DAS. Directly relevant to the ±25 m registration uncertainty in this project. |
| **Martin 2018**, Stanford dissertation | `04_methods_das/`<br>`martin_2018_thesis_passive-imaging-das.pdf` | Passive imaging and characterisation with DAS. Dissertations carry the acquisition detail that papers cut. |
| **Atterholt et al. 2024**, JGR 129 | `04_methods_das/`<br>`atterholt_2024_jgr_garlock-fault-zone-with-fiber.pdf` | Garlock fault zone imaged with fibre. Fault-zone structure from a DAS array. |

---

## 5. SAFOD / Parkfield — what is already known about this hole

*Site knowledge. Mostly 2004–2006, from the SAFOD drilling campaign.*

### Your fibre specifically

- **Lellouch, Yuan, Spica, Biondi & Ellsworth 2019**, JGR 124:6931 — `05_site_safod_parkfield/`<br>`lellouch_2019_jgr_velocity-estimation-downhole-das-SAFOD.pdf`. Velocity estimation from downhole DAS at SAFOD. **Same cable.** Source of the 864 m length, the 800 m analysis limit, the loop failure, and the 50–750 m velocity model.
- **Lellouch, Yuan, Ellsworth & Biondi 2019**, BSSA 109:2491 — `05_site_safod_parkfield/`<br>`lellouch_2019_bssa_velocity-based-detection-SAFOD.pdf`. Velocity-based earthquake detection on the same fibre.
- **Chavarria, Malin, Catchings & Shalev 2003**, Science — `05_site_safod_parkfield/`<br>`chavarria_2003_science_vsp-inside-san-andreas-parkfield.pdf` (two copies). VSP inside the San Andreas at Parkfield. Relevant to the check-shot registration work.

### Structure, stress and damage

- **Hickman & Zoback 2004**, GRL — stress orientations and magnitudes in the pilot hole
- **Boness & Zoback 2004**, GRL — stress-induced velocity anisotropy in the pilot hole
- **Townend & Zoback 2004**, GRL — regional tectonic stress near the SAF
- **Li et al. 2004**, GRL — low-velocity damaged structure from fault-zone trapped waves
- **Hole et al. 2006**, GRL — fault-zone structure from seismic refraction
- **Unsworth 2004**, GRL — magnetotelluric resistivity structure
- **Schleicher et al. 2006**, GRL — clay-coated fractures in mudrock fragments
- **Ellsworth & Malin**, Sibson volume — deep rock damage from P- and S-type guided waves
- **Rubinstein & Beroza 2005**, GRL 32:L14313 — depth constraints on nonlinear ground motion, **using the SAFOD target repeaters**. Concludes the change is confined above ~100 m. Bridges §1 and §3.
- **Cochran et al. 2025**, Seismica — continental scientific drilling for fault mechanics. Current framing of why holes like this exist; useful for an intro.
- *Drillbit seismic* (Oil & Gas Journal 2005) — trade article, background only.

---

## 6. Borehole and active-source methods

*Your fourth item. A method rather than a target, so it supports the others.*

- **Banerjee & Chatterjee 2021**, Near Surface Geophysics 20:710 — Stoneley/tube waves for fracture analysis. Relevant to the AWD tube-wave observations.
- *Also:* Chavarria et al. 2003 (VSP, §5) and Madsen et al. 2016 (depth calibration, §4).

---

## 7. Worth finding — gaps in the collection

Ordered by how much each would unblock.

### Method papers you're missing

| paper | why it matters |
|---|---|
| **Sens-Schönfelder & Wegler 2006**, GRL 33:L21302 | **The stretching estimator itself**, plus seasonal velocity variation at Merapi. Every dv/v measurement in this project uses this method; you should have the paper. |
| **Waldhauser & Ellsworth 2002**, JGR 107:2054 | The double-difference method and the repeater criteria (ΔM ≤ 0.3, source overlap) currently applied here second-hand. |
| **Waldhauser & Schaff 2008**, JGR 113:B08311 | The relocated NCSN catalog — the historical cross-check for your five sequences. |
| **Nadeau et al. 1994**, BSSA 84:247 | "Parkfield III" — the original HRSN repeater criteria that everything later cites without restating. |
| **Li & Zhan 2018**, GJI 215:1583 | DAS template matching at Brady. The detection method Phase 3 would use. |
| **Lellouch et al. 2021**, JGR 126 | FORGE low-magnitude seismicity with downhole DAS — per-channel correlogram stacking, 6×MAD thresholds, acausal validation. |

### For the tides direction specifically

- **Hillers et al. 2015**, JGR — in-situ velocity changes in response to tidal deformation. The most direct modern precedent after Takano.
- **Yamamura et al. 2003**, JGR — long-term in-situ velocity and attenuation, including tidal and seasonal terms.
- **Mao et al. 2019/2020** — time-frequency domain traveltime change measurement; matters when separating a 12.42 h tidal period from longer trends.

### For the ambient-noise direction

- **Brenguier et al. 2008**, Nature Geoscience — forecasting volcanic eruptions from noise (the companion to the Science paper you have).
- **Clements & Denolle 2018**, GRL — groundwater storage from ambient noise. Directly on the hydrologic interpretation.
- **Wang et al. 2017**, JGR — seasonal velocity changes throughout Japan; the large-N version of Li & Ben-Zion.
- **Shi et al. 2026**, Science — agroseismology, DAS on buried fibre resolving soil-moisture velocity change. Very recent and very close to the shallow-DAS question.
- **Kidiwela et al. 2026**, Sci. Adv. — Cascadia noise interferometry; uses dv/v as a strain proxy with ε = −(1/β)·dv/v, β ~10³–10⁴.

### For repeaters

- **Turner, Nadeau & Bürgmann 2013**, GRL — repeaters and creep rate, the applied version of Nadeau & McEvilly.
- **Schaff & Waldhauser 2005**, BSSA 95:2446 — CC-based differential times at NCSN; the windows and bands convention.

---

## Suggested reading order if starting cold

1. **Zhan 2020** (§4) — what DAS is
2. **Nadeau & McEvilly 1999** (§1) — what repeaters can measure
3. **Lellouch et al. 2019 JGR** (§5) — what this specific fibre has already done
4. **Brenguier et al. 2008** (§3) — dv/v monitoring at this site
5. **Takano et al. 2014** (§2) — the tidal version
6. **Uchida & Bürgmann 2019** (§1) — the caveats, once the rest makes sense
