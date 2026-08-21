# Geophysics audit — SAFOD ambient-noise cross-correlation

**Date:** 2026-08-20
**Scope:** claims A–F as stated by the PI, plus the seven specific physics
questions attached to them.
**Method:** every number below was recomputed from the stored product or from
first principles. Nothing was accepted from a log that reported it.
**Reproduce:** `safod-ambient-cc/scripts/audit_geophysics.py` (S1–S9) and
`safod-ambient-cc/scripts/audit_geophysics_supplement.py` (S10–S15). Job logs:
`notebooks/awd_clean/logs/audit_geo_40079562.out`,
`audit_geo2_40079138.out`, `audit_geo3_40080122.out`.

**Relation to `docs/AUDIT_CLAIMS.md`** (committed as `ad4664d` while this audit
was running). That document audits claim *provenance and wording*; this one
audits the *physics*. Where they overlap we agree — the p = 0.0005 resolution
floor, "consistent with" rather than "matches" for 1547 m/s, the 1702 flip being
bracketed rather than localised, and the src-211 figure importing src-400
statistics. **We disagree on one item.** `AUDIT_CLAIMS.md` §1 row 19 marks the
R±10 claim "UNSUPPORTED as evidence" and cites the rising downgoing/upgoing ratio
(1.64 → 1.97) as tension. Section E below resolves that tension: the ratio is a
*selected* statistic (§F) and carries no evidential weight, while two independent
quantitative metrics — usable range 506 m → 298 m and low-wavenumber energy share
20.2 % → 36.4 % — confirm the claim and its stated mechanism. I take
`AUDIT_CLAIMS.md`'s correction that **Nano virtual-source channels 10 and 40 are
in the air** (wellhead at channel 73, `nano_find_wellhead.txt`) and have applied
it below; it changes the verdict on claim F's Nano arm from "balanced" to "void".

Verdict key: **REFUTED** the claim as written is wrong · **WEAKENED** the claim
survives but is overstated or rests on an uncontrolled step · **CONFIRMED**
reproduced independently · **UNCHECKABLE** no data in this workspace can decide it.

---

## Summary table

| # | Claim | Verdict |
|---|---|---|
| A1 | An ambient arrival exists on the Deep fibre, p = 0.0005, 61/229 velocities clear | **CONFIRMED** |
| A2 | Its velocity is **1675 m/s** | **REFUTED** — direct picking gives 1443–1537 m/s on four gathers |
| A3 | causal/acausal **2.34** is evidence for that arrival | **REFUTED** — the grid-wide median is 2.007 and the ratio peaks at 825 m/s |
| A4 | Visible as moveout in a record section (`deep_section_depth_top211.npz`) | **WEAKENED** — visible in `deep_record_section.npz` (r = 0.979); the named product gives r = 0.137 |
| B1 | It matches the AWD **1547 m/s** mode | **WEAKENED** — 8.3 % apart, inside the 138 m/s scatter; band-matched value is 1592 m/s |
| B2 | 1547 m/s is an independent active-source measurement of a fluid-guided mode | **REFUTED** — it is a semblance pick confined to a 1300–1800 m/s tube-wave prior |
| B3 | ~1675 m/s is physically a tube wave | **WEAKENED** — exceeds plausible borehole-fluid velocity; 1443–1537 m/s does not |
| C | Not Lellouch's 3200 m/s; per-velocity p there is 0.20–0.36 | **CONFIRMED** — 0.1982 / 0.3677 / 0.2224 |
| D1 | Illumination asymmetry 0.348 (p = 0.005) in 2017 vs 0.040 (p = 0.73) in 2024–25 | **CONFIRMED** as arithmetic |
| D2 | 240-window scan: 11 significant vs 12.0 expected ⇒ no illumination | **WEAKENED** — the 240 p-values are non-uniform (KS p = 0.001), shifted high; the null is conservative and no power curve exists |
| D3 | Detectability exponent N^+0.042 shows nothing accumulates | **WEAKENED** — fitted at the 5950 m/s grid edge, not at the a-priori 3200 m/s |
| E | R±10 hurts on Deep; a spatial low-pass amplifies a low-k contaminant | **CONFIRMED** (mechanism measured: low-k share 20.2 % → 36.4 %); attribution refined |
| F | Directional f-k gives Deep 1.5–2× | **WEAKENED** — the ratio is selected to exceed 1; the direction label itself is correct |
| F′ | …Nano 0.985 (balanced) | **REFUTED / VOID** — that gather's virtual source is Nano channel 10, which is in the air |
| 1 | Units/differentiation dimensionally sound; gauge dismissed at 0.1–1.1 % | **CONFIRMED** (number) / **REFUTED** (stated formula, and its provenance) |
| 2 | 5–20 Hz and 15–30 Hz are adequately sampled at all three spacings | **CONFIRMED** — minimum λ/dx = 25 |
| 3 | Is the central negative weakened by the band? | **REFUTED for 2024–25**, **CONFIRMED for Deep** — see §3 |
| 4 | Geometry: columns, lead-in, near-vertical, reversal, TVD as moveout variable | **CONFIRMED**, with one new discrepancy (§4) |
| 5 | Tube-wave plausibility of 1675 vs 1547 | see B3; the 1675−1547 gap is **not significant** |
| 6 | ±12 ms gate, median across traces | **CONFIRMED** as harmless; a different bias found instead (§6) |
| 7 | Causal/acausal sign convention | **CONFIRMED** by synthetic round trip |

---

## A. The Deep-fibre arrival

Input: `awd_clean/ambient_transfer/deep_exact_stack/aggregate_deepA_src{400,800,1200}_ram0p1_cross_correlation_ordered_r0.npz`
(1200 windows, 300 × 60 s files = **5.00 h**, confirming the "5 hours" in the README).

### A1 — the arrival is real. CONFIRMED.

Recomputing `moveout_scores` from the stored `r_plus_minus_10_correlation`
reproduces the stored `causal_moveout_scores` to **8.9 × 10⁻¹⁶**. With an
independent seed and 4000 permutations (the product used 2000):

| src | peak score | at | per-velocity p | clearing own 95th pct | chance if independent |
|---|---|---|---|---|---|
| 400 | 3.3327 | 1675 m/s | 0.00025 | 59 of 229 | 11.5 |
| 800 | 2.2585 | 1525 m/s | 0.00025 | 38 of 229 | 11.5 |
| 1200 | 2.2520 | 1800 m/s | 0.00025 | 44 of 229 | 11.5 |

Two caveats on the p-value: 0.00025 = 1/4001 and the published 0.0005 = 1/2001
are both the **resolution floor**, so the honest statement is `p < 1/(N+1)`, not
`p = 0.0005`. And the "of 229" counts are not 229 independent tests, because the
score curve is heavily autocorrelated; the excess over 11.5 is nevertheless large.

It is **not** the zero-lag pedestal. Setting every gate to zero lag (v → ∞)
gives a score of **0.9565** for src 400, against a peak of 3.3327. Restricting to
offsets ≥ 250 m (9 traces, gates ≥ 149 ms from zero lag) still gives
3.7418 at 1675 m/s, p = 0.0012.

Independent corroboration on the return limb: `deep_section_depth_pair200out.npz`
(ch 308, outbound, TVD 201 m) and `deep_section_depth_pair200ret.npz` (ch 3092,
return limb, same TVD) both show a linear moveout with **r = +0.89 to +0.92** over
171 traces. Two physically different pieces of fibre at the same depth. This is
the strongest single piece of evidence in the project and it is currently
under-reported.

### A2 — the velocity is not 1675 m/s. REFUTED.

The project built `deep_record_section.py` explicitly so the moveout could be
*seen* rather than scored. Picking each trace's envelope maximum and regressing
picked lag on separation (`audit3.py` S14):

| gather | traces | sep window | velocity | intercept | r |
|---|---|---|---|---|---|
| `deep_record_section.npz` (dense, along fibre) | 172 | 100–450 m | **1443 m/s** | −0.007 s | +0.979 |
| `deep_record_section_pub.npz` (R±10) | 172 | 100–450 m | **1532 m/s** | +0.002 s | +0.960 |
| `deep_section_depth_pair200out.npz` (TVD, outbound) | 171 | 100–450 m | **1509 m/s** | +0.033 s | +0.916 |
| `deep_section_depth_pair200ret.npz` (TVD, return limb) | 171 | 100–450 m | **1485 m/s** | +0.032 s | +0.902 |

Theil–Sen agrees (1444–1651 m/s). Four gathers, two fibre limbs, two coordinate
axes: **1443–1537 m/s**, mean ≈ 1490 m/s. The published 1675 m/s is 9–16 % faster
than every direct measurement of the same wavefront.

**Why the scan is biased high.** `MAX_LAG_SECONDS = 0.35` and the published
geometry stops at 700 m offset. A gate at `offset / v` therefore falls outside the
recorded lag window whenever `v < offset / 0.362`:

| trial v | offsets inside ±0.35 s |
|---|---|
| 1100 m/s | 7 of 14 |
| 1500 m/s | 10 of 14 |
| **1675 m/s** | **12 of 14** |
| 1934 m/s and above | 14 of 14 |

The estimator is a **median over a changing subset**, and the subset that drops
out is always the far offsets — precisely the traces that carry the moveout. The
scan therefore prefers velocities fast enough to keep most offsets in the window.
Regressing the same 14 published offsets confirms the truncation: picks saturate
at exactly 0.350 s for the three longest offsets, and the fitted velocity inflates
to 1966–2113 m/s.

The peak velocity is also unstable under trivially defensible choices:
1675 → 2000 m/s on dropping offsets < 100 m (src 400); 1600 → 1700 m/s as the gate
half-width goes 50 → 4 ms; 1525 / 1675 / 1800 across three source channels
(sd 138 m/s).

**Recommendation:** report the arrival at **≈1.5 km/s (1443–1537 m/s by direct
picking on four gathers)**, and retire "1675 m/s", which appears hard-coded as
`V_MARK`/`V_ARR`/`--v-ref` in `deep_record_section.py:72`,
`deep_section_depth.py:53`, `deep_zerolag_vs_stack.py:56`, `fk_downgoing.py:67,145`
and every `jobs/fk_down_job.sbatch` line. Those figures currently draw a moveout
line that is measurably wrong.

### A3 — "causal/acausal 2.34" is not evidence. REFUTED.

Recomputed causal/acausal ratio for src 400 across the whole 300–6000 m/s grid:

* at the peak (1675 m/s): **2.338** — the quoted 2.34 is arithmetically right
* at 3200 m/s, where p = 0.198: **2.138**
* grid-wide: median **2.007**, 10th pct 1.223, 90th pct 2.286
* 76 % of all trial velocities have a ratio > 1.5
* the ratio is **largest at 825 m/s (3.321)**, not at the score peak

The gather's causal side is uniformly about twice its acausal side. Quoting the
value at 1675 m/s as though it discriminated that velocity is not supportable.
Same pattern at src 800 (peak 1.800, grid median 1.459, max at 1500 m/s) and
src 1200 (peak 1.466, grid median 1.440, max at 325 m/s). On the Nano gather the
same non-specificity appears: c/a = 2.580 at 2950 m/s and 2.609 at 3200 m/s.

**Provenance failure.** `docs/REPRODUCE.md:72` credits `scripts/deep_pervel.py`
with "c/a 2.34". That script never computes a causal/acausal ratio. No script in
either repository prints 2.34; it appears only as prose in `README.md:88`,
`docs/REPRODUCE.md:72` and a comment in `scripts/jobs/deep_ab.sbatch:22`.

### A4 — "visible as moveout in a record section". WEAKENED.

The claim names `deep_section_depth_top211.npz`. That product does **not** show
it: envelope-peak regression over separations > 100 m gives v = 9937 m/s,
r = **+0.137**. Its return-limb counterpart `deep_section_depth_top211ret.npz`
gives r = **−0.044** — nothing at all, which is the opposite of what the
return-limb control was set up to demonstrate. `deep_section_depth_src400.npz`
gives r = +0.340.

The moveout *is* clearly visible in `deep_record_section.npz` (r = 0.979) and in
the two `pair200` gathers (r ≈ 0.90). The claim should cite those.

Note also that the gate-scan statistic and direct picking disagree systematically
on the TVD sections: scan peaks at 1250–1475 m/s where picking gives 1485–1537 m/s.

---

## B. "It matches the 1547 m/s AWD mode"

### B2 — 1547 m/s is not an unconstrained measurement. REFUTED.

From `awd_clean/deep_dvv_frozen_trajectory.json`
(`"generated_by": "deep_dvv_injection_recovery.py --stage freeze"`):

```
slowness_search_range_s_per_m : [1/1800, 1/1300]   ->  1300 to 1800 m/s
outbound|15_30  1544.6 m/s   semblance 0.3303   role primary
return  |15_30  1549.7 m/s   semblance 0.2517   role primary
outbound|3_15   1591.8 m/s   semblance 0.2922   role secondary
return  |3_15   1565.2 m/s   semblance 0.1205   role secondary
```

"1547 m/s" is the arithmetic mean of the two 15–30 Hz legs. The search **could
not have returned a value outside 1300–1800 m/s**, and
`awd_clean/DEEP_DVV_STATUS.md:282-286` records that the window was chosen on
tube-wave grounds. A constrained estimate inside a tube-wave prior cannot then be
offered as evidence that the mode is a tube wave. `MANUSCRIPT_DISCUSSION.md:25`
already says "This is a proposed mechanism, not a measured one" — the ambient
write-up drops that hedge.

Two further points the ambient work should inherit:

* `MANUSCRIPT_RESULTS.md:82-88` records that per-aperture local picks span
  **1389–1562 m/s**, so 1547 is a single frozen global fit, not a tight number.
* `awd_clean/monitoring_sensitivity/LEGACY_RECONCILIATION.md:57,121` records that
  the independent reanalysis prefers **3–15 Hz at ~1600 m/s** for both branches
  and that "the preferred band does not" survive.

### B1 — is 1675 the same as 1547? WEAKENED.

| comparison | difference |
|---|---|
| ambient 1675 vs 1547.1 (15–30 Hz frozen mean) | +127.9 m/s, +8.3 % |
| ambient 1675 vs 1591.8 (3–15 Hz, the band that overlaps 5–20 Hz) | +83.2 m/s, +5.2 % |
| ambient 1675 vs ~1600 (clean reanalysis, 3–15 Hz) | +75 m/s, +4.7 % |

The scatter of the ambient estimate across three source channels is
sd = 138 m/s, which is **larger than the 128 m/s gap**. The two numbers are not
significantly different — and equally, they are not a demonstrated match.

If A2 is accepted, the picture improves rather than worsens: the corrected
ambient velocity **1443–1537 m/s straddles 1547 m/s**. The physical claim is
better supported than the project believes; the arithmetic backing it is wrong.

### B3 — is ~1675 m/s a physically possible tube wave? WEAKENED.

Low-frequency tube (Stoneley) wave in an open fluid-filled hole,
`v_t = v_f / sqrt(1 + ρ_f v_f² / μ)` with ρ_f = 1000, ρ_rock = 2500:

| v_f | Vs = 1000 | 1500 | 2000 | 2500 | 3200 m/s |
|---|---|---|---|---|---|
| 1450 | 1069 | 1237 | 1318 | 1361 | 1394 |
| 1500 | 1088 | 1268 | 1355 | 1402 | 1438 |
| 1600 | 1125 | 1326 | 1428 | 1483 | 1526 |

`v_t → v_f` only as μ → ∞, i.e. a perfectly rigid (well-cemented, steel-cased)
wall. **A tube wave cannot exceed the borehole-fluid velocity.** Fresh water at
in-hole temperature and pressure is ≈ 1480–1560 m/s; saturated brine reaches
≈ 1700–1800 m/s.

* 1547 m/s: comfortably attainable in a cased, cemented hole with water.
* **1675 m/s: requires a fluid faster than 1675 m/s (a brine) *and* a near-rigid
  wall.** Not impossible, but it is a stronger physical requirement than the
  project acknowledges, and 1443–1537 m/s requires none of it.

**UNCHECKABLE component.** The repository contains **no** borehole diameter,
fluid, casing or cement data. `deep_registration_forward_model.py:4-8,31-35` and
`NIGHTLY_STATUS.md:77` say so explicitly, and `deep_registration_forward_model.txt:8`
concludes "this result cannot distinguish fluid, casing, formation compliance, or
another guided mode". A quantitative tube-wave prediction for *this* hole cannot
be made here. The one contextual fact available is that Lellouch et al. describe
the main hole as carrying fibre "in a 0.9 mm steel tube cemented between casing
strings" (`MANUSCRIPT_METHODS.md:45-48`) — the Deep wireline hangs in whatever
fluid that cased hole contains, which nobody has recorded.

---

## C. Not Lellouch's 3200 m/s. CONFIRMED.

Recomputed per-velocity p at 3200 m/s from the aggregates, independent seed,
4000 permutations: **0.1982** (src 400), **0.3677** (src 800), **0.2224**
(src 1200). Matches the claimed 0.20–0.36. The velocity-scan peak is nowhere near
3200 for any Deep source channel.

---

## D. "No arrival on 2024–25; illumination is absent"

### D1 — the asymmetry numbers. CONFIRMED as arithmetic.

`products/interrogator_and_illumination_v2.txt` gives, at rank 0,
2024–25 |A| = 0.0399 (null95 0.1930, p = 0.7307) and 2017 |A| = 0.3478
(null95 0.2484, p = 0.0050). The README's "0.040 p = 0.73" and "0.348 p = 0.005"
are those rows. The 2017 windows are genuinely pre-event: arrivals picked at
4.100 s and 4.916 s of 5.00 s records, outside the 2.5 s analysis window.

Two asymmetries between the arms that are stated in the logs but not in the
README: the 2017 arm is **~5 s from two earthquake records** against 30 s per day
for 2024–25, and the 2017 control **loses significance at rank 4 and 8**
(p = 0.2718, 0.4489), i.e. the asymmetry it measures is itself low-rank.

### D2 — "11 significant vs 12.0 expected". WEAKENED.

Recomputed from `illumination_window_scan.npz`:

* 11 hits at α = 0.05 of 240 windows; Binomial(240, 0.05) mean 12.0, sd 3.38,
  95th pct 18. **P(X ≤ 11) = 0.459**, two-sided p = 0.917. The arithmetic is right.
* **But the 240 p-values are not uniform: KS D = 0.1243, p = 0.001**, and the
  departure is toward *large* p:

  | α | hits | expected |
  |---|---|---|
  | 0.01 | 3 | 2.4 |
  | 0.05 | 11 | 12.0 |
  | 0.10 | **18** | 24.0 |
  | 0.20 | **36** | 48.0 |

  A null that produces systematically inflated p-values is **conservative**, which
  means the scan is under-powered by an unquantified amount. Combined with the
  fact that the positive control at the same rank-2 setting clears only at
  p = 0.0299, the scan has **no demonstrated power** to detect illumination
  weaker than Lellouch's. "NO ILLUMINATED WINDOWS" should read "no window reached
  significance under a null that is measurably conservative".

### D3 — the convergence exponent. WEAKENED.

Refitting from `ambient_stack_convergence.npz` reproduces N^+0.042 (baseline) and
N^+0.019 (common-mode). But the stored curve carries **one velocity per stack
length, and the log records it as 5950–5975 m/s — the top edge of the 1500–6000
m/s grid.** The exponent therefore measures the convergence of the *grid maximum*,
which the project elsewhere declares a void statistic
(`deep_pervel.py:5-10`). **No stored array holds detectability at the a-priori
3200 m/s.** The test that matters — does detectability at 3200 m/s grow as
√N — has not been run on the 2024–25 archive.

### D4 — the same test has never been applied to the positive result. NEW.

`nano_long_stack.py` and `ambient_stack_convergence.py` both compute
**detectability = score / its own null 95th** and reject a result when the
exponent is flat. That statistic exists for the 2024–25 record (+0.042) and for
the 25 h Nano stack (+0.023). It **does not exist for the Deep 1675 m/s arrival**.

The only convergence data for Deep is the raw-score curve in
`products/deep_cc_steps.txt`. Refitting it:

```
N       1     2     4     8    16    32    64   128   239
score 0.959 1.931 2.277 2.075 2.297 2.364 1.862 2.762 2.941
raw-score exponent = N^+0.128
```

against **N^+0.158** for the 2024–25 baseline the project calls "FLAT: more data
does not help". The curve is also strongly non-monotonic. This does not refute
the arrival — the raw score is not detectability, and 1 h is a short lever — but
**the asymmetry is not defensible as it stands**: the project's own
disqualifying test is applied to its negatives and not to its positive.
Running `deep_zerolag_vs_stack.py`-style batch accumulation with a per-stack null
at 1500 m/s would close this in one job.

### D5 — the Nano 2950 m/s result is void, and separately null.

`docs/AUDIT_CLAIMS.md` §1 rows 1–2 establish that the virtual source used for the
Nano 2950 m/s claim is **channel 10, which is in the air**
(`nano_find_wellhead.txt`: entry at channel 73, "channels 0-72 … must not be
used"). That voids the headline. Two things are worth adding from the numbers:

**Even taken at face value, the 3 h scan was a deficit, not an excess.** Across
the five scanned source channels, **14 of 1145 velocity tests clear their own 95th
percentile — 1.22 %, against a nominal 5 %.** The a-priori p at 2950 m/s is
0.0474, 0.6384, 0.4239, 0.4090, 0.6309 across channels 10, 40, 80, 120, 160.
Šidák-correcting the smallest for five source positions gives **p = 0.2155**.

**Restricting to the three in-hole source channels makes it cleaner, not
weaker.** Channels 10 and 40 are in the air; 80, 120 and 160 are in the borehole,
and their a-priori p at 2950 m/s are **0.4239, 0.4090, 0.6309** — none remotely
significant, and their causal/acausal ratios at 2950 m/s are 1.010, 1.068, 1.053.
The 25 h stack (`nano_long_stack.txt`: p = 0.1621 at 2950 m/s, 0 of 229 clearing,
detectability exponent N^+0.023) was also run from channel 10 and so must be
repeated from channel ≥ 73 before it can be cited either way.

Direct picking on the dense Nano gather (`deep_record_section_nano.npz`,
592 traces) gives |r| ≤ 0.08 in every separation and lag window tested — no
moveout of any kind at 5–20 Hz — but that gather too has its source at channel 10
and is therefore uninformative about the cemented fibre.

---

## E. "The published R±10 sum hurts on Deep". CONFIRMED, attribution refined.

**The mechanism, measured directly.** A sum of 21 channels spaced dx is a spatial
boxcar with the Dirichlet response `sin(21 π f dx / v) / (21 sin(π f dx / v))`:

| geometry | 5 Hz | 10 Hz | 15 Hz | 20 Hz | first null |
|---|---|---|---|---|---|
| Deep, dx 2.0419 m, v = 1675 m/s | 0.973 | 0.896 | 0.775 | **0.622** | **39.1 Hz** |
| 2024–25, dx 1.0209 m, v = 3200 m/s | 0.998 | — | — | **0.971** | 149.3 Hz |

k = 0 — any static spatial pattern — is passed at exactly 1.000 in both cases.
So on Lellouch's geometry the operator costs 3 % of the arrival, and on Deep it
costs **38 % at the top of the band, with its first null at 39 Hz, inside the
useful bandwidth.** The claimed mechanism is exactly right.

**Measured on the gathers.** Share of gather energy at |k| < 1/200 m⁻¹:

| gather | low-k share | offset where smoothed gate-SNR falls below 1.5 |
|---|---|---|
| dense, single channel | **20.2 %** | never within 749 m |
| dense + common mode removed | 21.4 % | 506 m |
| published R±10 | **36.4 %** | 492 m |
| published R±10 + common mode removed | 36.6 % | **298 m** |

The R±10 sum **nearly doubles** the low-wavenumber energy share, 20.2 → 36.4 %.
The README's "~450 m → ~250 m" pair is dense+cm vs pub+cm; my independent metric
gives **506 m → 298 m**, a factor 1.70 against the claimed 1.8. Confirmed.

**Refinement.** The README's table now carries a blank row for "published R±10,
no common mode", but `deep_record_section_pub.npz` exists and scores 492 m —
indistinguishable from dense+cm at 506 m. The honest statement is that **both
operations cost range and they compound**: 749 → 492 (R±10) → 298 (R±10 +
common mode). That row can be filled in from the stored product.

**Resolving the tension the README now flags.** `README.md` §5 and
`docs/AUDIT_CLAIMS.md` §1 row 19 both note that R±10 lowers the moveout energy
(3.366 → 2.057) but *raises* the downgoing/upgoing ratio (1.64 → 1.97), and treat
that as evidence against the claim. It is not evidence either way: §F shows the
ratio is chosen to be the larger of two reciprocals, so it measures the imbalance
of a gather that R±10 has already made noisier, not the arrival's strength.
The two metrics that are not selected — usable range and low-k energy share —
both move in the direction the claim predicts, and the array response predicts
the size of the effect from first principles. **The claim stands; it is the
downgoing/upgoing counter-evidence that should be withdrawn.**

One structural caveat: at 1675 m/s an offset of 700 m maps to 0.418 s, beyond
`MAX_LAG_S = 0.35`. **No configuration can test range beyond ~605 m** in these
products, so "~450 m" is measured against a ceiling of 605 m, not against 750 m.

---

## F. "Directional f-k gives Deep 1.5–2×, Nano 0.985". WEAKENED.

**The numbers reproduce**, and the direction label is physically right, but the
ratio is close to guaranteed by construction.

`fk_downgoing.separate(keep="auto")` chooses which wavenumber quadrant to *call*
downgoing by measuring envelope energy along the assumed moveout trajectory — the
same trajectory the reported `moveout_energy` metric then uses. Forcing both
choices:

| input | keep = positive-k | keep = negative-k | reported |
|---|---|---|---|
| `deep_record_section.npz` | 0.608 | **1.644** | 1.644 |
| `deep_record_section_cm.npz` | 0.639 | **1.564** | 1.564 |
| `deep_record_section_pubcm.npz` | 0.507 | **1.972** | 1.972 |
| `deep_record_section_nano.npz` | 1.016 | **0.985** | 0.985 |

The two rows are exact reciprocals. The script keeps the half the data say is
larger, so "downgoing/upgoing > 1" carries no information by itself; only the
*magnitude* of the imbalance and the input-level null do. (The Nano row is the
one case where the selection statistic and the reported statistic disagreed,
which is why 0.985 < 1 was reported at all.)

**What rescues it.** I derived the correct quadrant independently. A synthetic
gather with `t = +z/v` (downgoing from a source at z = 0) puts
9.50 × 10⁷ of its energy at k·f < 0 against 3.30 × 10⁶ at k·f > 0 — a factor 29.
In `numpy`'s `fft2` convention **downgoing is k·f < 0**, which is
`fk_downgoing`'s `mask_neg`, which is `keep = "negative-k"` — the choice reported
for **all four** products. So the label is correct and the Deep/Nano contrast
(1.56–1.97 against 0.985) is a real, correctly-signed contrast.

**Recommendation:** fix the convention a priori from the FFT convention, delete
`keep="auto"`, and report the ratio without selection. The result will not change.

### F′ — the Nano arm is void, not "balanced". REFUTED.

`deep_record_section_nano.npz` has `source_channel = 10`. `nano_find_wellhead.txt`
places the borehole entry at **channel 73** on three agreeing diagnostics
("channels 0-72 look like SURFACE/AIR fibre and must not be used"). The 0.985 is
therefore a cross-correlation from a virtual source **in the air**, and it says
nothing about the cemented fibre's downgoing/upgoing balance. It must not be
reported as the Nano counterpart to the Deep 1.5–2×. The correct experiment —
`fk_downgoing` on a Nano gather with the source at channel ≥ 73 — has not been
run.

Two further caveats on F: the Nano arm is computed at **5–20 Hz**
(`deep_cc_steps.BAND`, inherited through `steps.correlate`) with `v_ref = 2950`,
which is a 30–60 Hz velocity — see §3; and the Deep arm uses `v_ref = 1675`, which
§A2 shows is 9–16 % too fast. Re-running the Deep arm at 1490 m/s is a one-line
change and should be done before the number is quoted again.

---

## 1. Units and scaling. CONFIRMED (differentiation) / mixed (gauge).

**Differentiation.** `ambient_lellouch2019_exact_stack.py:602-605`:

```python
rate[:, 0] = 0.0
np.subtract(raw[:, 1:], raw[:, :-1], out=rate[:, 1:])
rate *= np.float32(fs)
```

A first difference multiplied by the sample rate is `d(phase)/dt`, i.e. one clean
factor of fs, no `dx` (correct — the derivative is temporal, not spatial), and no
doubling. Same code in `deep_cc_steps.differentiate`,
`deep_record_section.py:153`, `deep_section_depth.py:119`. Verified equal to the
engine to a relative 7.2 × 10⁻⁸ by `deep_cc_steps.py --verify`.

The absolute scale is in any case irrelevant downstream: the running-absolute-mean
divides each trace by its own moving |amplitude|, so the `× fs` is cosmetic. Note
this means the products are **not** calibrated strain rate and no amplitude in
them has physical units.

**Gauge length. The number is right; the formula and its provenance are not.**

The DAS gauge response is a boxcar average of length L along the fibre,
`sinc(L/λ) = sin(π f L / v) / (π f L / v)`. The repository states it as
`sinc(pi k L)` (`fig7c_negative/DOCUMENT.md:170`, `:762`, and
`sections/02_data_methods.md:20-23`, `sections/03d_illumination.md:78-79`). With
k = 1/λ that evaluates `sin(π² L/λ)/(π² L/λ)` — **wrong by a factor of π in the
argument**, and it has a first null at λ = π L rather than λ = L. **REFUTED as
written.**

The quoted 0.1–1.1 % is nevertheless correct, as the *differential* between the
two gauge lengths at 3200 m/s:

| f | 10 m (2017) | 16.335 m (2024–25) | difference |
|---|---|---|---|
| 5 Hz | 0.999598 (0.040 % loss) | 0.998929 (0.107 % loss) | **0.067 %** |
| 10 Hz | 0.998394 | 0.995719 | 0.268 % |
| 20 Hz | 0.993587 (0.641 %) | 0.982943 (1.706 %) | **1.064 %** |

So "0.1–1.1 % at 3200 m/s" is right to the stated precision, and the cross-epoch
dismissal stands. Three problems remain:

1. It is **computed by no code**. `grep -rn "np.sinc"` over `awd_clean/*.py`
   returns nothing; the string is a hard-coded literal in
   `lellouch2017_window_audit.py:201`, `ambient_apparent_velocity_census.py:307-308`,
   `ambient_fixed_k_test.py:17-18`, `build_ambient_fk_qc_notebook.py:46,1049,1153`
   and six markdown files. `fig7c_negative/REPRODUCE.md:236` lists it with **no
   producing script**, uniquely in that file.
2. `DOCUMENT.md:169` says "We quantify the gauge-length difference in section 3.4".
   **Section 3.4 is "Adaptive f-k".** The quantification does not exist.
3. It is a *differential* at one velocity in one band, and it is quoted elsewhere
   as if it were a general dismissal. The **absolute** losses are much larger
   where it matters: Nano 2026 at 2950 m/s loses 4.6 → **17.4 %** across 30–60 Hz,
   and the 2024–25 fibre at 2950 m/s would lose 4.5 → **17.2 %** if that band were
   ever analysed. Deep at 1547 m/s over 15–30 Hz loses 1.6 → 6.3 %.

The manuscript's separate first-null argument is correct: 2950/16.459 = 179 Hz and
1547/10 = 155 Hz (`MANUSCRIPT_METHODS.md:38-43`), both far above the analysis
bands.

---

## 2. Wavenumber and aliasing. CONFIRMED — nothing is aliased.

Nyquist wavenumber is 1/(2 dx); a wave is aliased only if λ < 2 dx.

| fibre | dx | 2 dx | worst case tested | λ | λ/dx | k/k_Nyq |
|---|---|---|---|---|---|---|
| Deep | 2.0419 m | 4.08 m | 1547 m/s @ 30 Hz | 51.6 m | **25.3** | 0.079 |
| Deep | 2.0419 m | 4.08 m | 1675 m/s @ 20 Hz | 83.8 m | 41.0 | 0.049 |
| Nano 2026 | 1.26606 m | 2.53 m | 2950 m/s @ 60 Hz | 49.2 m | **38.8** | 0.052 |
| 2024–25 | 1.0209 m | 2.04 m | 3200 m/s @ 20 Hz | 160 m | 156.7 | 0.013 |
| 2024–25 | 1.0209 m | 2.04 m | 2950 m/s @ 60 Hz | 49.2 m | 48.2 | 0.042 |

The tightest case in the entire project is 25 samples per wavelength, 12.6× the
Nyquist limit. **Both 5–20 Hz and 15–30 Hz are adequately sampled at all three
spacings, with an enormous margin; so would 30–60 Hz be.** Spatial aliasing is not
an available explanation for anything here, and `deep_timeseries.py:26-27` states
the same conclusion correctly.

The operators applied to these data are far more band-limiting than the sampling:
the R±10 sum has its first null at 39 Hz on Deep, and the 16.459 m gauge is 13
channels wide at the 2026 Nano spacing, so adjacent Nano channels are not
independent measurements.

---

## 3. The frequency band. REFUTED for 2024–25; CONFIRMED for Deep.

This was flagged as likely the most important item. It splits.

### The central negative is *not* weakened by the band. REFUTED.

The worry is that 5–20 Hz misses the cemented fibre's mode, which the manuscript
places "near 2950 m s⁻¹, strongest at 30–60 Hz". But the project's own
measurement is dispersive, and it is dispersive in the direction that rescues the
test. `awd_clean/nano_mode_identification.txt` reports, by band:

```
15-25 Hz: 3300 [3250,3350] m/s
35-45 Hz: 2950 [2950,2950] m/s
45-60 Hz: 2950 [2950,2950] m/s
60-80 Hz: 2950 [2950,2950] m/s
slowness trend 0.483 [0.414, 0.549] us m^-1 Hz^-1
```

Extrapolating that measured trend from 2950 m/s at 60 Hz down into the ambient
band:

| frequency | with 0.414 | 0.483 | 0.549 µs m⁻¹ Hz⁻¹ |
|---|---|---|---|
| 20 Hz | 3102 | 3128 | 3154 m/s |
| 12.5 Hz | 3132 | 3164 | 3196 m/s |
| 5 Hz | 3162 | **3201** | 3238 m/s |

**In the 5–20 Hz band the cemented fibre's own active-source dispersion predicts
3.10–3.24 km/s — i.e. Lellouch's 3200 m/s**, and the 15–25 Hz direct measurement
of 3300 m/s independently agrees. The ambient test was therefore run at a
band-appropriate target velocity. 2950 m/s is a 30–60 Hz number that should never
have been used as the a-priori ambient target in the first place (see D5).

Two residual caveats, both real but secondary: (i) the *amplitude* of the mode
may still be lower at 5–20 Hz even if its velocity is right, and nothing in the
tree measures that; (ii) `ambient_coherence_spectrum.txt` shows inter-channel
coherence length is 3.6–4.5 m in every band from 0.5 to 50 Hz, reaching 50 m only
at 0.15–0.60 Hz which carries 0.0 % of the power — so the archive looks
unpromising at *every* band, not just 5–20 Hz.

### For the Deep fibre the band is genuinely mismatched. CONFIRMED.

The ambient Deep result is 5–20 Hz; the active-source Deep mode is quoted at
15–30 Hz. The overlap is 15–20 Hz only. The band-matched active-source value is
the 3–15 Hz leg, **1591.8 m/s**, and the independent reanalysis prefers 3–15 Hz at
~1600 m/s for both branches (`LEGACY_RECONCILIATION.md:57,121`). Note also that
`deep_tube_validation.txt` shows the outbound leg's beamforming power is **higher
at 3–15 Hz (0.23474) than at 15–30 Hz (0.21768)**, so "strongest at 15–30 Hz" is
not uniformly supported by the product that is cited for it.

### The band is structurally frozen. NEW.

`ambient_lellouch2019_exact_stack.py` — the reference correlator — exposes 30
command-line arguments and **not one of them changes the band**;
`OUTPUT_BAND_HZ = (5.0, 20.0)` at line 55 is a module constant. Every downstream
script hard-codes the same pair. `deep_timeseries.py` (written 2026-08-20) is the
first script with a `--band` switch and states the problem correctly in its
docstring, but **it has produced no output files** in either repository.

---

## 4. Geometry. CONFIRMED, with one new discrepancy.

### Column identification, verified without reference to `safod_geometry.py`

`SAFOD_Phase2_GeoReferenced_Channels.xlsx`: 3201 data rows, columns
A, C–M (there is no column B, so no column can be mis-assigned).

* **C/D and E/F are 3.280840 to 1.3 × 10⁻¹¹** — the ft/m ratio, as claimed.
* A is 0…3200 in unit steps: the channel index.
* **The decisive test:** integrating `cos(G) · dD` along the outbound leg
  reproduces column F to a maximum error of **0.0060 m over a 2547.5 m interval
  (0.0002 %)**. F is TVD and G is inclination in degrees. Confirmed.
* Integrating `sin(G) · dD` reproduces M to 18.6 m against a 1112.5 m departure
  (1.7 %, the residual being azimuth curvature): M is horizontal departure.
* I/J are UTM E/N (720807–721472, 3983660–3984556); K/L are lat/lon
  (35.974–35.982 N, −120.552 to −120.545) — Parkfield. Correct.
* The cached `deep_channel_depth.npz` is identical to a fresh parse in both repos.

### The docstring's structural claims. All CONFIRMED.

* surface lead-in: 222 channels (0–210 at the head, 3190–3200 at the tail)
* channel 211: MD 2.05 m, TVD 2.05 m, inclination 0.05° — first in-hole
* near-vertical (outbound, inc < 5°): **channels 211–949, contiguous, TVD 2.0–1513.0 m**
* reversal at **channel 1700, MD 3052.9 m, TVD 2549.5 m**

### Is ΔTVD the right moveout variable? Yes, and 5° is tight enough.

For a **vertically propagating body wave**, ΔTVD is exactly right. For a
**borehole-guided mode** the correct variable is ΔMD, and the two differ by
cos(inclination). Over channels 211–949 the inclination is median 1.80°, mean
1.99°, max 4.98°:

* cos(4.98°) = 0.996222 → worst-case single-channel bias **0.379 %**
* the 5° cutoff admits at most **0.382 %** velocity bias
* integrated, ch 400 → 949: ΔMD 1124.85 m vs ΔTVD 1123.87 m, ratio 0.99913 — a
  guided wave read on a TVD axis appears **0.09 % slow**

Both are two orders of magnitude below the 9–16 % discrepancy found in §A2, so
the choice of axis is not the source of that problem. **The 5° cutoff is more than
tight enough.** It would not be outside the near-vertical section: at the 49.6°
and 54.3° of channels 1200 and 1600 the same confusion would cost 35–42 %, which
is why `safod_geometry.py`'s exclusion of those source positions is correct.

### NEW: the spreadsheet and the interrogator disagree on channel spacing by 0.343 %.

```
|diff(MD)| along the outbound leg : 2.048911 m   (a single unique value)
HDF5 Acquisition/SpatialSamplingInterval : 2.04190469 m
ratio 1.003433
```

Every along-fibre offset in the project (`geometry()` in the engine,
`deep_record_section.py`) is built from the HDF5 value; every depth is built from
the spreadsheet. The two coordinate systems are on scales that differ by 0.343 %.
The effect on any velocity here is ≈ 5 m/s and does not change a conclusion, but
it should be resolved and stated rather than left implicit, because it means
`offsets_m` and `tvd` are not the same metre.

### NEW: `aggregate_deepA_src1600` is not interpretable.

Source channel 1600 with 50–700 m offsets places receivers at channels 1624–1943,
which **crosses the fibre reversal at 1700**. That gather mixes outbound and
return-limb receivers and its acausal score exceeds its causal score. It should be
excluded, not read as part of a "sign flip at channel 1702" trend. (For the same
reason, the "flip at 1702" claim in `build_ambient_cc_steps_notebook.py:440`
cannot be supported by the stored aggregates, whose source positions are
98/400/800/1200/1600/2000/2400.)

---

## 5. Tube-wave plausibility. See B3.

Summary: 1547 m/s is physically comfortable; 1675 m/s is at or above the fluid
velocity and requires a brine; **1443–1537 m/s, the directly measured value, is
the most comfortable of the three.** The 1675 − 1547 = 128 m/s gap is smaller than
the 138 m/s scatter across source channels and is **not significant**. No borehole
fluid/casing/diameter data exist in this workspace, so tube-wave *identity* is
UNCHECKABLE here.

---

## 6. The moveout statistic. Gate and median are fine; the truncation is not.

**±12 ms gate, 5–20 Hz.** The band's dominant periods are 50–200 ms and the
envelope width is ≈ 1/BW = 67 ms, so a 24 ms gate samples the arrival's peak
without averaging over it. Sweeping the half-width on the src-400 aggregate:

| half-width | peak | at | p | clearing |
|---|---|---|---|---|
| 4 ms | 3.2846 | 1700 m/s | 0.0017 | 61 |
| 8 ms | 3.2795 | 1700 m/s | 0.0017 | 59 |
| **12 ms** | 3.3327 | 1675 m/s | 0.0017 | 58 |
| 20 ms | 3.3591 | 1625 m/s | 0.0017 | 58 |
| 30 ms | 3.2848 | 1650 m/s | 0.0017 | 60 |
| 50 ms | 3.1269 | 1600 m/s | 0.0017 | 62 |

The detection is insensitive to the gate; the velocity drifts 1700 → 1600 m/s,
a further ±3 % on the estimate. **12 ms is defensible.**

The gate does not bias toward a velocity, but its *velocity tolerance* scales as
v²Δt/x, so the effective bin is much coarser at high velocity and short offset:
±673 m/s at (1675 m/s, 50 m) and ±2458 m/s at (3200 m/s, 50 m). That is why the
"clearing range 1100–2750 m/s" is so wide and should not be read as a resolution.

**The median across traces is not pathological.** Mean gives 3.4594 at 1675 m/s
(same velocity), 20th percentile gives 2.4801 at 1575 m/s. The arrival is on most
traces, so the median is a reasonable robust choice.

**What is pathological is the NaN bookkeeping** — see §A2. `np.nanmedian` over a
subset that shrinks from 14 to 7 traces as the trial velocity falls is the
mechanism that biases the velocity estimate high. Two fixes, either sufficient:
raise `MAX_LAG_SECONDS` above 700/1400 ≈ 0.50 s, or restrict the scan to
velocities at which all 14 offsets are in range (v > 1934 m/s) and report the
result from the dense section instead.

**NEW code discrepancy.** The engine's `moveout_scores`
(`ambient_lellouch2019_exact_stack.py:1009-1010`) returns NaN when the gate falls
outside the lag window. `deep_cc_steps.moveout_curve:169` instead does
`k = argmin(|lags - x/v|)`, which **clamps to the nearest edge sample**. The two
are not the same estimator. This is why the `crop_then_band` ablation in
`products/deep_cc_steps.txt` "peaks" at 3.9106 at **300 m/s** — at that velocity
every gate is clamped to lag = ±0.35 s, where a filter applied after cropping has
its edge transient. That row says nothing about band/crop order and should be
withdrawn or recomputed with NaN semantics.

---

## 7. Causal/acausal convention. CONFIRMED.

Synthetic round trip through the exact code path
(`conj(rfft(src)) * rfft(rcv)` → `irfft` → `fftshift`), with the receiver
waveform being the source delayed by 0.120 s:

```
peak of the correlation is at lag +0.1200 s
```

**Positive lag = energy propagating from the virtual source to the receiver.**
With a source above its receivers that is downgoing, so `moveout_scores(sign=+1)`
is correctly the causal/downgoing branch and `sign=-1` the acausal branch. The
`fftshift` centring (`correlation_from_spectrum`, midpoint = `n_fft // 2`) is
consistent with the `lags` array it returns. No sign error found in the ambient
chain.

The f-k half of the convention is also correct — see §F: downgoing is k·f < 0 in
`numpy`'s `fft2` convention, verified on a synthetic gather (energy ratio 29:1),
and that is the quadrant every product keeps.

Given the project's documented history of sign errors (`BRANCH_LAG_SIGN`,
`dvv_core.py`), the synthetic in `ambient_lellouch2019_exact_stack.synthetic_validation`
should be extended to assert the f-k quadrant too, so the convention is pinned by
a test rather than by a data-driven choice at runtime.

---

## What to do next, in priority order

1. **Retire "1675 m/s."** Replace it with the directly measured 1443–1537 m/s
   throughout `README.md`, `docs/REPRODUCE.md`, and the `V_MARK` / `V_ARR` /
   `--v-ref` constants in five scripts. Every figure currently overlays a moveout
   line 9–16 % too fast.
2. **Withdraw "causal/acausal 2.34" as evidence**, or replace it with the
   grid-wide profile, which shows the causal side is uniformly ~2× the acausal.
3. **Run the detectability-convergence test on the Deep arrival** at a fixed
   1500 m/s, with the null rebuilt at each stack length. That is the one test the
   project uses to kill results and has not applied to this one.
4. **Fix the moveout estimator**: NaN-vs-clamp in `deep_cc_steps.moveout_curve`,
   and raise `MAX_LAG_SECONDS` so the 700 m offset is reachable below 1934 m/s.
5. **Correct the gauge-length text**: the response is `sinc(L/λ)`, the 0.1–1.1 %
   figure is a *differential*, and `DOCUMENT.md`'s forward reference to
   "section 3.4" points at the wrong section. Add three lines of code that
   compute it.
6. **Re-run every Nano product from a source channel ≥ 73.** The 3 h scan, the
   25 h stack, the dense record section and the `fk_downgoing` Nano arm were all
   driven from channel 10, which is in the air. Until then the Nano fibre has no
   ambient result of any sign, and "Nano 0.985 (balanced)" must not be quoted.
   Restore the R±10-without-common-mode row in README §5 from
   `deep_record_section_pub.npz`, and withdraw the downgoing/upgoing ratio as
   counter-evidence there (§E).
7. **State the 1300–1800 m/s constraint** wherever 1547 m/s is cited, and quote
   the band-matched 3–15 Hz value (1592 m/s) when comparing with a 5–20 Hz
   ambient result.
8. **Exclude `aggregate_deepA_src1600`** from any interpretation; its receiver
   span crosses the fibre reversal.
9. **Add a power statement to the illumination scan.** The p-values are
   significantly non-uniform toward large values (KS p = 0.001), so the null is
   conservative and the negative is weaker than "11 vs 12.0 expected" implies.
10. **Resolve the 0.343 % channel-spacing disagreement** between the spreadsheet
    (2.048911 m) and the HDF5 header (2.0419047 m).
