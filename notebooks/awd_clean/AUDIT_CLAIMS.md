# Claim audit — ambient-noise CC documents and the public repository

**Audited 2026-08-20.** Every quantitative claim in the documents below was traced
to the `.npz` / `.json` / `.txt` product that produced it, the number re-read from
that product, and the *strength* of the wording checked against what the number
supports.

Scope: `awd_clean/AMBIENT_FIG7C_STATUS.md`, `AMBIENT_FIG7C_STATUS_ADDENDUM.md`,
`AMBIENT_FIG7C_MULTIDAY.md`, `FIG7C_MULTIDAY_RESULT.md`, `AMBIENT_LOWK_MECHANISM.md`,
`AMBIENT_CC_LITERATURE_REVIEW.md`, `AMBIENT_FK_EXPLAINER.md`; and in
`safod-ambient-cc`: `README.md`, `docs/REPRODUCE.md`, `docs/fig7c_negative/*`,
`figures/README.md`, `figures/FIGURE_NAMING.md`, and the nine committed figures.

**`AMBIENT_FIG7C_STATUS.md` was not edited.** Its section 4 was inspected and is
flagged only. In fact nothing in that file needed correction: every number in it
verifies (see §6 below).

Verdict key: **SUPPORTED** · **OVERSTATED** (number right, wording outruns it) ·
**UNSUPPORTED** (no product) · **STALE** (superseded by later work) ·
**CONTRADICTED** (two places in the tree disagree, or the product says otherwise).

---

## 1. OVERSTATED / UNSUPPORTED / CONTRADICTED — read these first

| # | CLAIM | DOCUMENT:LINE | SUPPORTING PRODUCT | VERDICT | WHAT I ACTUALLY FOUND |
|---|---|---|---|---|---|
| 1 | "5 velocities clear at 2,900–3,000 m/s; **p = 0.0474** at the a-priori 2,950" | `safod-ambient-cc/README.md:122` (§4) | `products/nano_ambient_cc.txt`, block `--- source channel 10 ---` | **CONTRADICTED + STALE** | The number is real (0.0474, 5 of 229, band 2900–3000) but **the virtual source is Nano channel 10, which is in the air.** `nano_find_wellhead.txt` puts the borehole entry at channel 73 from three agreeing diagnostics. The *same repository* already says so: `figures/README.md:12` — "This is why the Nano 2,950 m/s result was withdrawn." The front page contradicted its own figure index. **FIXED.** |
| 2 | "A 25 h stack **is running** with a convergence test at 2,950 m/s" | `safod-ambient-cc/README.md:129` | `awd_clean/nano_long_stack.txt` (completed 17:26) | **STALE** | It finished and it is null: 25.00 h, 5,950 windows, **p = 0.1621** at 2,950 (was 0.0474 at 3.0 h), p = 0.3965 at 3,200, peak 1.211 at 4,325 m/s (p = 0.6359), **0 of 229** velocities clear, convergence **N^+0.023**. Also run from air channel 10. **FIXED.** |
| 3 | "the same picker … **returns the 2,416–4,357 m/s profile his Figure 9 reports** (depth–velocity r = 0.948)" | `safod-ambient-cc/README.md:47-49` | `awd_clean/ambient_lellouch_fig7d_profile.txt`, "fit r = 0.948" | **OVERSTATED (live error)** | r = 0.948 is `np.corrcoef(depth, velocity)` on **our** picks from **his** released 7d traces — our picks' monotonicity. His Figure 9 curve is not digitised anywhere in either repository. `docs/FIG7C_MULTIDAY_RESULT.md:295-298` and `docs/AMBIENT_LOWK_MECHANISM.md:103-104` both record this wording as *withdrawn*; the front page still carried it. Note the `.npz` key `r` is a **different quantity: −0.4453**, the velocity–depth trend of our own 2024-25 picks. **FIXED.** |
| 4 | Same Figure 9 wording, uncaveated | `docs/fig7c_negative/sections/00_abstract.md:15-16`; `01_introduction.md:38-39`; `04_discussion.md:115-116`; `make_figures.py:13` | same | **OVERSTATED (live error)** | Four uncaveated instances, including the abstract — the highest-exposure text in the set. A further four instances carry a correction that disclaims only *point-by-point* matching while leaving "the range his Figure 9 reports" / "the published velocity range" standing; those assert range agreement against the same undigitised curve. **ALL FIXED** (sections + `assemble.py` caption + `make_figures.py` docstring). |
| 5 | "the interrogator … stable cross-day spatial pattern is confined to the surface lead-in, **channels 0-22 (\|corr\| 0.8426 vs control 0.3889)**" | `AMBIENT_LOWK_MECHANISM.md:165`; `AMBIENT_CC_LITERATURE_REVIEW.md:359`; `AMBIENT_FK_EXPLAINER.md:201` | `awd_clean/interrogator_blame_test.txt` | **UNSUPPORTED as attributed** | 0.8426 / 0.3889 is v1's run over **channels 0–699, "lead-in INCLUDED"** (that file's own line 4). **No channels-0–22-only correlation was ever computed.** The table presented it as one half of a matched pair against the 23–708 row, which *is* real. The "it's the lead-in" inference rests instead on the T2 power-ratio test, whose own output says "a high ratio alone is NOT proof of an instrumental origin … does not decide on its own", and whose quoted 2e4–1.3e5 band excludes 2 of its 6 days (2024-05-21 = 323, 2025-03-13 = 0.01, i.e. inverted). **FIXED** (row relabelled; the load-bearing 23–708 result is untouched and stands). |
| 6 | "k = 0 holds **~66 %**, against **0.13 %** in Lellouch's 2017 records" | `AMBIENT_FK_EXPLAINER.md:150, 171`; `AMBIENT_LOWK_MECHANISM.md:154` | none | **UNSUPPORTED** | **66 % appears in no product in the tree.** The matched measurement — both arms through `ambient_apparent_velocity_census.py`, same band, same k = 0 handling — gives **52.12 % vs 0.04 %** (`ambient_apparent_velocity_census.txt:26`). The 0.13 % comparand came from `lellouch2017_window_audit.txt`, a different script with different processing: unmatched arms, the exact error this project says caused five of its ten withdrawals. **FIXED.** |
| 7 | "detectability **never reaches the 1.00** a detection needs" | `AMBIENT_CC_LITERATURE_REVIEW.md:232-233` | `ambient_stack_convergence.txt:13` | **CONTRADICTED by its own product** | The baseline branch reaches **detect 1.014, p = 0.0170 at 16 chunks**. The interpretation survives — that excursion peaks at 5,950 m/s with pedestal +0.983, i.e. the scan-ceiling artefact — but the sentence as written is false. **FIXED.** |
| 8 | "Stacking 24× more data … improved it by about **4 %**" | `AMBIENT_CC_LITERATURE_REVIEW.md:231-232` | `ambient_stack_convergence.txt:9-14` | **CONTRADICTED (arithmetic)** | 24^0.042 = 1.14, and the measured detectability goes 0.858 → 0.971 = **+13 %**. The "4 %" is the exponent 0.042 misread as a percentage. **FIXED.** |
| 9 | Post-filter null row "**0.222 \| 0.098 \| 0.002**" | `AMBIENT_FK_EXPLAINER.md:118` | none | **UNSUPPORTED** | None of the three numbers occurs anywhere in the tree. The matched product for the same date and same 300 files (`ambient_transfer/interim_2024-12-20_n300.json`) gives **0.2789 / 0.1022 / 0.0005**. The dependent "factor of **3.7**" (:125) is 3.6 on the real numbers. **FIXED.** |
| 10 | "day selection (**the richest day scored *worst***)" | `AMBIENT_FK_EXPLAINER.md:202` | `FIG7C_MULTIDAY_RESULT.md:37-53, 92-103` | **STALE + CONTRADICTED** | Doubly wrong. The richest day by census (2024-11-30, 4.9 %) scored **best** (p = 0.1345) and was excluded from the four-day set; and the "census fan %" ranking that defines "richest" is **itself withdrawn** for having no geometric baseline. Day selection is still closed, but by the six-day sweep and the 96 h stack, not by this comparison. **FIXED.** |
| 11 | QC row "Stack convergence, **5 h → 8 days** … N^+0.042" | `AMBIENT_FK_EXPLAINER.md:62` | `build_ambient_fk_qc_notebook.py` §8 vs §13.2 | **CONTRADICTED (two experiments merged)** | QC §8 is a day-combination *waveform-correlation* test over 5 h → 8 days and reports **no exponent**. N^+0.042 is QC §13.2 / `ambient_stack_convergence.txt`, which is **one day in 1 → 24 hourly chunks**. **FIXED** (split into two rows). |
| 12 | "The Deep arrival, **correct geometry**: fixed-top source at the wellhead …" | `figures/README.md:8` (fig 01) | `products/deep_section_depth_top211.npz` | **OVERSTATED** | The npz confirms source channel 211, TVD 2.05 m, 1.6 h, 383 windows, receivers TVD 4–700 m — the geometry claim is right. But "**The** Deep arrival" imports the statistics (p ≤ 0.0005, c/a 2.34, 61 of 229) that exist **only for source channel 400**, and **no per-velocity null has been run in the src-211 geometry.** Figure 01 is a visual claim. The same file's line 21 removes other figures precisely because src-400 "is not a Figure 7c geometry". **FIXED** (both caveats added). |
| 13 | "the causal/acausal ratio flips sign at channel **1702**, **exactly** where the fibre reverses" | `README.md:109-111` | `ambient_transfer/deep_exact_stack/aggregate_deepA_src*_…_r0.npz` | **OVERSTATED** | The ratio at peak really does cross 1: **2.34 / 1.80 / 1.47 / 1.45** at outbound ch 400 / 800 / 1200 / 1600 against **0.656 / 0.522** at return-limb ch 2000 / 2400. But the scan has only those positions, so the crossing is bracketed between **1600 and 2000** — about 800 m of fibre — not localised to a channel. Inclination also co-varies with limb across that gap. Separately the depth registration puts the reversal at **1700** (`deep_channel_depth.npz`, `FIGURE_NAMING.md`), while the older scripts say 1702. **FIXED.** |
| 14 | Deep table row "**1200** \| 1800 m/s \| 0.0005 \| 1.47 \| 44 of 229" | `README.md:106` | `products/deep_pervel.txt`; `scripts/safod_geometry.py:22-27` | **OVERSTATED** | Numbers verify exactly. But channel 1200 sits at **49.6° inclination**, and `safod_geometry.py` states in its own docstring that at that inclination "their along-fibre apparent velocity is **NOT** a vertical velocity" and that "only sources 400 … and 800 … were geometrically valid". Presented in one table with two valid positions and no marker. **FIXED** (row flagged, TVD/inclination column added). |
| 15 | "per-velocity p = **0.0005**" (three times) | `README.md:104-106`; `docs/REPRODUCE.md:72` | `deep_pervel.py:41` `NULLS = 2000` | **OVERSTATED (precision)** | 1/2001 = 0.00049975 is the **resolution floor**: it means no permutation of 2,000 exceeded the observation, i.e. p ≤ 0.0005. Stating it as an equality implies a measured value. Same for the f-k null's "p = 0.005" (`NULLS = 200`, floor 1/201). **FIXED** (both stated as floors). |
| 16 | "It **matches** instead the **1,547 m/s** mode" | `README.md:114` | `deep_pervel.txt`; `manuscript/MANUSCRIPT.md:30` | **OVERSTATED (mild)** | The measured Deep peak is **1,675 m/s** against the AWD manuscript's ~1,547 m/s — **8.3 %** apart. Both are dispersive and the AWD mode is "strongest at 15–30 Hz" while this work is 5–20 Hz. "Consistent with" is the supportable word. **FIXED**, and the "this is NOT Lellouch's 3,200 m/s P wave" statement was promoted to bold leading text (see §3). |
| 17 | "**Cause: illumination.**" | `README.md:51` | `interrogator_and_illumination_v2.txt`; `illumination_window_scan.txt` | **OVERSTATED (mild)** | The measurements support "the asymmetry Lellouch used as his own evidence is present in 2017 and absent in 2024-25". Calling that *the cause* is a step past what an absence measurement gives; the parent doc says "the **binding constraint**". Note also the 2024-25 |A| = 0.040 headline is **one 30 s window on 2025-01-17** against **~5 s** of 2017 from two records. **PARTLY FIXED** ("Most likely cause"); the single-window disclosure is a **flagged item**, see §2. |
| 18 | "The **archive-wide** scan is **complete** … no illuminated window was found **anywhere**" | `AMBIENT_LOWK_MECHANISM.md:23-28` | `illumination_window_scan.txt` | **OVERSTATED** | 240 × 30 s = **2.0 h of a 377,680-file (6,294.7 h) archive = 0.032 % coverage**, and the product's own LIMITS say "a burst shorter than 30 s could be missed". `AMBIENT_CC_LITERATURE_REVIEW.md:213` gets this right ("anywhere **in the sampled archive**"). **FLAGGED**, not edited — the wording is the PI's own framing and the fix is a scope clause, not a number. |
| 19 | "the R±10 sum hurts: moveout visible to ~450 m → **~250 m**" | `README.md:206-212` (§5) | none — the figures were culled | **UNSUPPORTED as evidence** | This is an eyeball reading of figures that `figures/README.md:22` deliberately removed. The quantitative products point partly the other way: R±10 lowers moveout energy at 1,675 m/s (3.366 → 2.057) but **raises** the downgoing/upgoing ratio (1.64 → 1.97). Also "four variants" preceded a three-row table. **FIXED** (products added alongside, tension stated, fourth row restored). |
| 20 | `p = 0.0002` presented as positive evidence | `awd_clean/README.md:74`; `awd_clean/NIGHTLY_STATUS.md:193` | `ambient_transfer/fk_seasonal_aggregate.json` | **STALE** | "permutation p=0.0002 … a stable 70-minute selected observable" is the 2.5–4.5 km/s fan result that `Ambient_FK_QC_workflow.ipynb` **rejected**, and it is the very max-over-grid statistic that moved to 0.1596 on a wider scan. **FLAGGED** — these two files are outside the audited set and are the PI's own older status headers. |

---

## 2. Flagged, not edited

1. **`AMBIENT_FIG7C_STATUS.md` section 4 — inspected, no change needed, not touched.** The three verdicts (common-mode removal is not a published Figure 7c step; Equation 6's estimator is underspecified; the R±10 sum is necessary but insufficient) all stand, and the sqrt(21) ≈ 4.6 arithmetic is correct. This matches `AMBIENT_FIG7C_STATUS_ADDENDUM.md` §1.
2. **The single-window basis of the illumination headline.** |A| = 0.040, p = 0.7307 for 2024-25 is one 30 s window on 2025-01-17; the 2017 positive control is ~5 s from two earthquake records. `AMBIENT_CC_LITERATURE_REVIEW.md:220-222` discloses this; `AMBIENT_LOWK_MECHANISM.md:12-13` and `AMBIENT_FK_EXPLAINER.md:183-185` quote the numbers without it. The PI should decide whether the headline carries the disclosure.
3. **The 2017 positive control fails in the lower half of its own band.** `ambient_directional_asymmetry.txt:8` puts the 2017 arm at |A| = 0.0007, p = 0.9975 in **5–12 Hz** (significant only in 12–20 Hz). `AMBIENT_LOWK_MECHANISM.md:20-22` cites that script alongside v2 without noting it.
4. **`p = 0.0002` in `awd_clean/README.md:74` and `NIGHTLY_STATUS.md:193`** — see row 20. Both are outside the audited set.
5. **`AMBIENT_FIG7C_STATUS_ADDENDUM.md` §2.6** says the QC notebook is now v11; `build_ambient_fk_qc_notebook.py` confirms `# v11`. The `v8` marker in `AMBIENT_FIG7C_STATUS.md`'s Provenance table is stale — the addendum already asks the PI to decide, so it is left alone.
6. **`AMBIENT_FIG7C_MULTIDAY.md`** is correctly banner-marked SUPERSEDED and its census column withdrawn. No further action; do not cite.
7. **Per-velocity nulls are missing exactly where positives appear.** `docs/fig7c_negative` declares the per-velocity null as its fix for max-over-grid p-values but reports it for **one** run (the four-day stack, "0 of 181"). The three positive-looking numbers — day-pairs p = 0.0390 / 0.0165, adaptive f-k p = 0.0060, 16-chunk p = 0.0170 — carry only the deprecated statistic. Each is repudiated in prose on other grounds. A limitation bullet recording this was added; running the per-velocity null on those three is the real fix.
8. **The 2017 illumination positive control is the best of a five-value rank sweep** (ranks 0/1/2 significant, 4/8 not) with no familywise statement. p = 0.0050 would survive Bonferroni over five, but it is never said. Note p = 0.0050 = 1/200 exactly, i.e. the null's resolution floor, while the 2024-25 arm's p-values are not on a 1/200 grid — worth confirming the two arms used the same null construction, given how often this document asserts "identical processing".

---

## 3. The two distinctions that must stay unmissable

**A — the Deep ~1,675 m/s arrival is not Lellouch's 3,200 m/s P wave.** Different
mode (fluid-guided, not body), different fibre (2026 **wireline** Deep, not the
cemented fibre Lellouch used), different band optimum. Per-velocity p at
3,200 m/s is 0.20–0.36, i.e. nothing. The disclaimer existed in `README.md:97`
but sat *after* a section heading reading "An arrival IS recovered on the Deep
fibre" and a sentence saying ambient noise "independently **reproduces**" an
active-source measurement, in a README whose lede is "an attempted reproduction
of Lellouch et al. (2019) Figure 7c". **FIXED** — the disclaimer is now the
paragraph's bold opening sentence.

One numeric collision to be aware of: `docs/fig7c_negative` uses "1,675 m/s" for
something else entirely — the scan-**floor** artefact on 2024-05-11 in the
2024-25 Nano archive (`sections/03_results.md:70`, `REPRODUCE.md:68`). Same
number, different fibre, different meaning, no cross-reference either way.
**FLAGGED.**

**B — the central negative is scoped.** One band (5–20 Hz), one geometry (fixed
source at channel 23, receivers 50–700 m in 50 m steps), one archive (2024-05 →
2025-05). Before this audit: the archive scope was stated everywhere; the **band**
was in `docs/fig7c_negative` §5.4 / §6.1 but missing from Conclusion 1, the
abstract's operative sentences, and the public README; the **geometry** was stated
once in methods and never as a scope limit, with no bullet in the limitations
list. **FIXED** in `README.md` §1, `sections/05_limitations.md` (Conclusion 1 and a
new §6.1 bullet), and `FIG7C_MULTIDAY_RESULT.md`'s conclusion.

The document is careful never to say Lellouch is wrong, and it independently
recovers his asymmetry from his own data — that framing is correct and unchanged.
Its title and Conclusion 7 do make an unqualified general claim ("borehole
ambient-noise body-wave interferometry should be treated as illumination-limited")
but argue it from a stationary-phase mechanism rather than from the null alone,
which is a defensible basis. **FLAGGED, not changed** — that is the PI's thesis.

---

## 4. Figures — labelling, provenance, and dead references

Convention (`scripts/safod_geometry.figure_label`): every figure states **hours
stacked** and the date range in **both local and UTC**, local first.

| figure | source | hours | local | UTC | verdict |
|---|---|---|---|---|---|
| `01_deep_gather_depth` | `deep_section_depth.py --source 211`, 1.6 h, 383 win | yes | yes | in-title | SUPPORTED; **stale render** (byte-differs from the current script's output; the current script moves UTC to a footnote and spells the source channel out) |
| `02_deep_gather_wiggle` | same | yes | yes | in-title | same |
| `03_deep_fk_downgoing` | `fk_downgoing.py` on `deep_record_section.npz` (**source ch 400**) | yes | yes | in-title | numbers SUPPORTED (1.644 / 1.564 / 1.972; p = 0.0050 floor of 200). **Title does not state the source channel**, so it reads as the same geometry as 01. Also a stale render. |
| `04_nano_no_moveout` | `deep_record_section.py --fibre nano --source 10` | **NO** | **NO** | **NO** | **VOID as a Nano borehole result** — source channel 10 is in the air. Also **has no suptitle at all**; panel (a) reads "592 traces, 719 windows", i.e. counts not duration, from a superseded script version. Carries a "2950 m/s (recovered)" moveout line — the withdrawn velocity. |
| `05_nano_wellhead` | `nano_find_wellhead.py`, 0.50 h | yes | **NO** | yes | UTC-only title. **FIXED in the script**; the committed PNG still needs regenerating. Its `.txt` is not in the public repo. |
| `06_illumination_2017_vs_2024` | `interrogator_and_illumination_v2.py` | **NO** | **NO** | **NO** | numbers SUPPORTED; no duration/date header |
| `07_illumination_archive_scan` | `illumination_window_scan.py` | **NO** | **NO** | **NO** | numbers SUPPORTED (11 / 240 / 12.0 / ceiling 18); no header |
| `08_stack_convergence` | `ambient_stack_convergence.py` | **NO** | **NO** | **NO** | numbers SUPPORTED; no header |
| `09_deep_per_velocity_nulls` | `deep_pervel.py` (**sources 400 / 800 / 1200**) | **NO** | **NO** | **NO** | numbers SUPPORTED; no header, and no indication the sources are down the hole or that 1200 is deviated |

No "nan h" title survives — `fk_downgoing.py` reprints duration from `n_windows`
when a section product predates the `hours` key, which is the case for
`deep_record_section.npz`.

All nine figures exist and all nine are referenced by `figures/README.md`. In
`docs/fig7c_negative/figures/` all eight `.png` are referenced; the eight `.pdf`
are orphans (print copies, harmless).

**Dead references found:** `figures/FIGURE_NAMING.md:28` pointed at
`figures/corrected/README.md`, deleted by the figure cull `fff585b` — **FIXED**,
replaced by an air-channel table in place. `docs/fig7c_negative/REPRODUCE.md`
referenced `FIGURE_NUMBERS.txt` ("the machine-checked half of this file") and
`_afk_unit_check.txt`, neither of which had been copied into the public repo —
**FIXED**, both copied in from `awd_clean/`.

`FIGURE_NAMING.md`'s tag table describes tags that **no figure in the directory
uses any more** after the rename to `01_`…`09_`, and its header claimed titles
"now spell the configuration out in full", which the committed PNGs do not.
**FIXED.**

---

## 5. Reproducibility breaks in the public repo

1. **`ambient_lellouch2019_exact_stack.py` is not in the repository.** Eight
   scripts import it by name — `ambient_stack_convergence.py`, `deep_cc_arms.py`,
   `deep_cc_steps.py`, `deep_pervel.py`, `deep_record_section.py`,
   `deep_section_depth.py`, `deep_zerolag_vs_stack.py`, `nano_long_stack.py` — so
   every step of `docs/REPRODUCE.md`'s build order raises `ModuleNotFoundError`.
   Every `scripts/jobs/*` launcher sidesteps this by `cd`-ing to
   `notebooks/awd_clean` first, which means the jobs do **not** exercise this
   repository's copies. **DOCUMENTED**, not fixed (the engine is large and is the
   parent tree's).
2. **`scripts/safod_geometry.py` looked for its spreadsheet and depth cache beside
   itself**, while both are committed under `products/`, so `geo.load()` raised
   `SystemExit` and figures 01/02/03 could not be regenerated. **FIXED** — it now
   searches `products/` first with a fallback to its own directory.
3. **`docs/fig7c_negative/make_figures.py` cannot run from its current location**
   (`AWD = HERE.parent` resolved to `notebooks/awd_clean/`, now resolves to
   `docs/`; the aggregate directory it needs is not in the repo), and
   `docs/fig7c_negative/REPRODUCE.md:8` still says "all paths are relative to
   `notebooks/awd_clean/`". **DOCUMENTED**, not fixed.
4. **Air-channel guards are not uniform.** `deep_section_depth.py` refuses a
   non-downhole source. `nano_ambient_cc.py`, `nano_long_stack.py` and
   `deep_record_section.py` do not, and `deep_cc_arms.py` still defaults
   `--source-channel` to 98 (Deep surface lead-in). **DOCUMENTED**; the guard
   should be lifted into the shared reader.

---

## 6. What verified exactly — recorded so it is not re-audited

- **`AMBIENT_FIG7C_STATUS.md`, whole file.** Re-read from
  `aggregate_2024-12-20_src23_…_ordered_r0.npz`: grid 181 velocities, 1500–6000
  at 25 m/s; peak **6.1314 at 5850**; at 3200 causal **2.7523**, acausal
  **2.8314**; scan-max null 95th **6.3241**; familywise **p = 0.1470** over 10,000
  permutations. Matches the section 3 table and the section 1 prose exactly.
- **`FIG7C_MULTIDAY_RESULT.md`.** Coherent stack: 4 days, 5,760 files, **23,036
  windows, 96.0 h**, peak **1.9121 at 5925**, at 3200 causal 1.1434 / acausal
  1.1556 = **0.99**, null95 **2.6493**, **p = 0.9184**, detectability **0.72** —
  all exact against `ambient_lellouch2019_multiday_stack.txt`. The ten-pair survey
  table is exact against `fig7c_pair_survey.txt`, including 0.0390 and 0.0165 and
  "2 of 10 against 0.5 expected".
- **Fixed-k mechanism.** Band centres 8.50 / 16.00 Hz ratio **1.882**; |k|
  centroids 0.00173 / 0.00175 ratio **1.011**; 1.21 / 1.23 cells; **39.33 % /
  39.08 %** of in-band energy — all exact against `ambient_fixed_k_test.txt` and
  `FIGURE_NUMBERS.txt`.
- **Adaptive f-k, all five configurations and twenty cells** — re-derived from the
  aggregate `.npz`, exact, including p = 0.0060 with pedestal +0.987 at 5900 m/s
  and the clean config at 1.078 / 3300 / p = 0.9392 / pedestal −0.036, and the
  6.57× (quoted 6.6×) pedestal amplification.
- **Illumination.** 2017 |A| = 0.3478, null95 0.2484, p = 0.0050; 2024-25
  |A| = 0.0399, null95 0.1930, p = 0.7307; archive scan 11 hits of 240 against
  12.0 expected, Binomial 95th = 18 — exact.
- **Convergence.** N^+0.042 / N^+0.019 baseline / common-mode, raw N^+0.158 /
  N^-0.004 — exact.
- **Deep per-velocity.** src 400 peak 3.3327 at 1675 (61 of 229, range
  1100–2750); src 800 2.2585 at 1525 (39 of 229); src 1200 2.2520 at 1800 (44 of
  229); curve ends 0.786 / 0.784 for src 400; p at 3200 = 0.2049 / 0.3608 /
  0.2344. Causal/acausal at peak re-derived from the aggregates: **2.338 / 1.800 /
  1.466** → the README's 2.34 / 1.80 / 1.47.
- **The "5 hours"** of Deep ARM A: `deep_srcscan.sbatch` uses `--nfiles 300` ×
  60 s = 5.0 h, first record 2026-06-15 23:09:48 UTC = 16:09 local. Exact.
- **Familywise arithmetic.** 1 − 0.95^17 = 0.5819, matching the "58 % over ~17
  aperture positions" in `docs/REPRODUCE.md:91-94`.
- **The `.npz` key `r` = −0.4453** with p = 5e-5 against a depth-label permutation
  null — the velocity–depth trend of our own 2024-25 picks, and a genuinely
  different quantity from the 0.948 in the `.txt`. Both real; only the conflation
  was wrong.

---

## 7. Prioritised remaining work

**Do before presenting.**

1. **Regenerate figures 01, 02, 03, 05 and re-commit.** All four are stale renders;
   05's title fix is in the script but not in the PNG. The `safod_geometry.py` path
   fix in §5.2 is what unblocks this.
2. **Regenerate or retire figure 04.** It is the only figure whose *claim* is void:
   air source, no title, counts instead of duration, withdrawn velocity line. Run
   `deep_record_section.py --fibre nano --source 73` (or below) to get the contrast
   with 01 that the caption wants.
3. **Add a duration/date header to figures 06–09** via `geo.figure_label`. Four of
   nine figures currently carry no indication of how much data they represent.
4. **Rerun the Nano pre-drop CC from a downhole source.** Nano has never been
   tested with a valid one. This is the single highest-value missing measurement:
   it is the only 2026 test of the *cemented* fibre and it is currently absent.

**Do before submitting anything.**

5. **Run a per-velocity null in the src-211 fixed-top Deep geometry.** The headline
   figure has no statistic and the statistics belong to a different geometry.
6. **Run per-velocity nulls on the three positive-looking max-over-grid numbers**
   (day-pairs, adaptive f-k, 16-chunk) so the document's declared statistic covers
   the places where it matters.
7. **Reconcile 0.1470 vs 0.1614 and −0.381 vs −0.219** in `docs/fig7c_negative` — a
   §6.4 note now records both, but one canonical null should be chosen.
8. **Regenerate or drop the unsourced pedestal `+0.951`** at `sections/04_discussion.md:47-48`.
9. **Lift the air-channel guard into the shared reader** so no script can source
   from the lead-in again.

**Housekeeping.**

10. Fix `docs/fig7c_negative/{REPRODUCE.md, make_figures.py}` paths, or state that
    the figure build runs only from the parent tree.
11. Decide the four wording questions `AMBIENT_FIG7C_STATUS_ADDENDUM.md` §5 puts to
    the PI (line-3 supersession scope; re-running the 900-channel common-mode branch
    over 23–708; section 6's "next test" paragraph; the v8 → v11 marker).
12. Retire `p = 0.0002` from `awd_clean/README.md:74` and `NIGHTLY_STATUS.md:193`.
