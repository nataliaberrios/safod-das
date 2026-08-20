# Addendum to `AMBIENT_FIG7C_STATUS.md` — for Natalia's sign-off

**Written 2026-08-19 by the assistant. `AMBIENT_FIG7C_STATUS.md` itself was left
completely untouched, deliberately.** That file is Natalia's own document. It is
not edited here, not annotated in place, and no banner was inserted into it. This
separate file exists so that the parts of it affected by later work are flagged
without anyone editing her text on her behalf. Nothing below is a correction that
has been *applied*; everything below is a flag awaiting her decision.

`git log --follow AMBIENT_FIG7C_STATUS.md` will show no commit from this pass.

---

## 0. First, a discrepancy in the instruction that she should see

This addendum was requested with the instruction that **"section 4 contains
assistant-produced census numbers that are now withdrawn."** That does not match
the file. Counting the `##` headings:

| # | Heading | Contains census numbers? |
|---|---|---|
| 1 | Question and answer | no |
| 2 | Exact input and baseline | no |
| 3 | Full-day branch results | no |
| 4 | **Evaluation of the three disputed claims** | **no** |
| 5 | Validation gates | no |
| 6 | Interpretation and next test | no |
| 7 | Provenance | no |

**There are no energy-census numbers anywhere in `AMBIENT_FIG7C_STATUS.md`.** By
`##` count, section 4 is "Evaluation of the three disputed claims" (common-mode
removal, Equation 6, the R±10 sum), and every claim in it still stands — the
sqrt(21) ≈ 4.6 arithmetic included. Rather than guess which block was meant and
flag the wrong thing, this addendum flags **everything in the file that later work
touches**, in file order, in §2 below. Sections 1 and 4 come out clean.

The withdrawn census numbers she may be thinking of are real, but they live in
**other** files. See §3.

---

## 1. What stands, unchanged and reinforced

State this plainly first, because it is most of the document.

- **The headline answer — "Answer for 20 December 2024: no" — stands**, and is now
  much better supported than when written. Six complete days give a minimum
  p = 0.1345; Fisher's combination across days gives p = 0.524; a coherent
  four-day 96.0-hour, 23,036-window stack gives p = 0.9184. Independently, **0 of
  181 velocities clear the per-velocity null** — at 3,200 m/s, 2.752 against a
  threshold of 2.811 — so the negative does not rest on the familywise correction
  at all.
- **Every number in the section 3 branch table stands.** An independent chunk
  re-sum reproduces the stored aggregates to max abs difference 0.000e+00. Two
  reading caveats have been added since, both in §2 below; no value changed.
- **All five validation gates in section 5 stand.**
- **Section 4's three verdicts stand**, including the judgement that common-mode
  removal must not be inserted into the paper baseline and that Equation 6's
  estimator is underspecified.
- **Section 6's warning was right, and a later assistant claim contradicting it
  was withdrawn.** She wrote: *"The result at 5.85 km/s should not be interpreted
  as a high-velocity arrival: it lies near the flat-moveout boundary and survives
  receiver scrambling."* That is correct and is now confirmed mechanically:
  `corr(trial velocity, score) = +0.976`, moving the grid cap moves the "peak" with
  it (3,500 → 3,475; 4,000 → 3,925; 8,000 → 7,700; 20,000 → 18,775 m/s), and a
  moveout-free control reproduces the observed curve to within 3.5 % at every
  velocity. A later assistant-written sentence in `FIG7C_MULTIDAY_RESULT.md` and in
  the QC notebook — *"this archive does contain a statistically detectable coherent
  component"* — went the other way and has been **withdrawn 2026-08-19**. Her
  hedge was the correct call.

---

## 2. What in her file is affected, in file order

Six items. None invalidates a measurement; four are supersessions and two are
newly-understood caveats on how a number should be read.

### 2.1 Header, line 3 — "supersedes every earlier Figure 7c verdict in this repository"

**Status: was true when written, now out of date in scope. Her call.**

Three later documents now carry parts of this question, and by this repo's own
convention the later file wins on the topic it covers:

| file | written | owns |
|---|---|---|
| `AMBIENT_FIG7C_STATUS.md` (hers) | 2026-08-14 | the **single-day 2024-12-20** operator, branches and gates — still authoritative |
| `FIG7C_MULTIDAY_RESULT.md` | 2026-08-14, reconciled 08-19 | the **six-day and coherent-stack** result |
| `AMBIENT_LOWK_MECHANISM.md` | 2026-08-19 | the **static fixed-k contaminant** (mechanism of the processing failures) |
| `AMBIENT_CC_LITERATURE_REVIEW.md` | 2026-08-19 | the **illumination** finding (the binding constraint) and the adaptive-f-k result |

Suggested wording, if she wants it: "supersedes every earlier *single-day* Figure
7c verdict in this repository." No number changes.

### 2.2 Section 3 table, "Best velocity (m/s)" column — read as scan ceilings

**Status: values correct; interpretation now pinned down.** The 5,850 / 5,925 /
4,550 / 3,650 entries are where the velocity grid stops, not velocity estimates,
for the reason in §1 above. Her section 6 already says this for the 5.85 km/s
figure specifically; the finding generalises to the whole column. The p values are
unaffected — the receiver-order permutation preserves the gate geometry and carries
the identical bias (97 % of null curves also peak at ≥ 5,500 m/s).

### 2.3 Section 3 table, "Common mode, median of all 900 channels" branch — NEW CAVEAT

**Status: p = 0.9220 unchanged; the estimator is now known to be contaminated.**
This is the one item worth her attention on the science rather than the wording.

Two later measurements bear on it:

- **Channels 0–22 carry ~97 % of the array's 5–20 Hz energy.** They are the
  uncemented surface lead-in.
- **They are also the only place with a stable cross-day spatial pattern** —
  |corr| 0.8426 against a control of 0.3889, i.e. an instrumental fingerprint.
  Inside the analysed aperture, channels 23–708, the dominant pattern is **not**
  stable (0.0188 against a control 95th percentile of 0.1282), so it is not an
  instrumental fingerprint there. Source: `interrogator_and_illumination_v2.txt`.

So a median taken over **all 900** channels estimates a common mode that is
dominated by 23 channels the analysis aperture deliberately excludes. The branch's
negative conclusion is unaffected in direction, and the Figure 7c *baseline* is
untouched (it has no common-mode step at all, correctly, per her section 4.1). But
if this branch is ever re-quoted as "common-mode removal does not reveal the
target", the honest version restricts the estimator to channels 23–708 or states
the caveat. The 900-channel branch was not re-run.

### 2.4 Section 6 — the "clean next test" has been run

**Status: superseded, all three parts.** She proposed: freeze the operator, apply
it to independently chosen complete days, then a predeclared multi-day convergence
sequence. All three were done.

| her proposed test | result | product |
|---|---|---|
| independently chosen complete days | six days, min **p = 0.1345**; Fisher **p = 0.524** | `FIG7C_MULTIDAY_RESULT.md` |
| predeclared multi-day stack | four-day coherent 96.0 h, 23,036 windows: **p = 0.9184** | `ambient_lellouch2019_multiday_stack.txt` |
| convergence sequence | detectability **N^+0.042** (baseline), **N^+0.019** (common-mode removed), against **N^+0.50** required | `ambient_stack_convergence.txt` |

The convergence exponent is the sharpest version of the negative: stacking 24×
more data should have improved detectability 4.9-fold and improved it by about 4 %.

One caution she should know about before reading that table: **one selected
day-pair reaches p = 0.039, and it is a selection artefact, not a result.** The two
days were chosen *because* they scored lowest; the pair peaks at the scan ceiling
with causal/acausal 0.97. It must not be quoted as a partial reproduction.

### 2.5 Section 6 — "F–K-assisted correlations remain a separate extension"

**Status: stands, and the verdict has since come in — negative.** Per
`Ambient_FK_QC_workflow.ipynb`: the 2.5–4.5 km/s fan produces an apparent ridge
but **fails the pre-filter channel-scramble gate**, so it is not accepted as an
independently recovered physical arrival. Her requirement that such a result carry
matched input-level controls is exactly what rejected it.

That verdict was re-tested from a second direction and held. Adaptive f-k (Isken
et al. 2022) was run in five configurations; none passed. The configuration with
**no** prior static removal reached p = 0.0060 — by amplifying the fixed-k pedestal
6.6-fold, while carrying the worst pedestal diagnostic of the five and peaking at
5,900 m/s. The mechanism generalises and is worth recording next to her section 6:
**a receiver-order null permutes the finished gather, so any operator applied
before the gather is formed sits outside its own null.** The cleanest
configuration (common mode + rank-2 + AFK) drove the pedestal to −0.036, moved the
peak into the physical range at 3,300 m/s, and gave **p = 0.9392**.

### 2.6 Section 7 Provenance — "Advisor notebook ... v8"

**Status: stale version marker only.** The notebook is now **v11**
(v8 → v10 → v11). Contents and location unchanged. Nothing else in the provenance
table is affected; the four commit hashes and three SLURM job IDs were not
re-checked against the cluster but nothing in this pass touched them.

---

## 3. The withdrawn census numbers — where they actually are

For completeness, since the instruction pointed at them. They are **not** in her
file. The energy census (`ambient_fk_energy_census.py`) is withdrawn because it had
no geometric baseline: 98.4 % of in-band (f,k) cells lie below 1,500 m/s by
construction, so "81–99 % of energy below 1,500 m/s" is at or below the white-noise
expectation. Measured as density per cell the data are 2.6–2.8× *enriched* at
body-wave velocities. Three further defects: the k = 0 column silently dropped 22 %
of band energy; there was no spatial taper, so the fan measurement is largely
Dirichlet leakage from the k ≈ 0 peak; and any run without the `_ch23-896` suffix
is determined by the 23 lead-in channels. Separately, "the downgoing share never
exceeds 49.9 %" is simply **false** — the tree's own outputs include 51.6 %.

Where those numbers appear, and what was done in this pass:

| location | action |
|---|---|
| `FIG7C_MULTIDAY_RESULT.md` | withdrawal banner already present since 2026-08-14; left in place, and a second banner added for the two claims withdrawn on 08-19 |
| `AMBIENT_FIG7C_MULTIDAY.md` | **marked SUPERSEDED**; its "census fan %" column and the section built on it flagged as void |
| `build_ambient_fk_qc_notebook.py` §11 | withdrawal text already present; retained and extended |
| `ambient_fk_energy_census_*.txt` | left alone — derived products, regenerable, and the script's own caveats now cover them |

---

## 4. Other claims withdrawn on 2026-08-19, for her awareness

None of these come from her document; they are listed so she is not caught out by
a number she may have seen in a commit message or an older notebook cell.

| withdrawn claim | why |
|---|---|
| "81–99 % of energy below 1,500 m/s implies absence" | no geometric baseline |
| "downgoing never exceeds 49.9 %" | false; 51.6 % exists in this tree's own output |
| "2024-25 decorrelates in 4 m vs 26 m in 2017" | the two arms were processed differently |
| "a coherent component at 5,850 m/s" | scan ceiling; a moveout-free control reproduces it |
| "tapering buys no discriminating power" | unmatched real and synthetic paths; corrected ratio 1.51–2.01, not 0.63 |
| "the contaminant and target are unresolved at this aperture" | withdrawn **in part** — the aperture arithmetic stands, but the target **is** resolved above ~12 Hz |
| "a manned field deployment supplied the 2017 illumination" | speculation; no site log, no operational record, no statement in the paper |
| "our picker recovers Lellouch's published Figure 9 at r = 0.948" | overstatement. r = 0.948 is `corr(depth, velocity)` for our picks on **his released traces** — our picks' monotonicity, not agreement with his published curve |

Also ruled out as explanations, rather than left open: **gauge length** (10 m vs
16.335 m is a 0.1–1.1 % effect at 3,200 m/s) and **the interrogator inside the
analysed aperture** (§2.3).

---

## 5. What she needs to decide

1. Whether to narrow the line-3 supersession claim to "single-day" (§2.1). Wording
   only.
2. Whether the 900-channel common-mode branch should be re-run over channels
   23–708 (§2.3). This is the only item that could change a number in her table.
3. Whether section 6's "next test" paragraph should be replaced by a pointer to
   `FIG7C_MULTIDAY_RESULT.md`, or left as the historical record of what was
   proposed before it was run (§2.4).
4. Whether the v8 marker in Provenance should be bumped to v11 (§2.6).

Until she says otherwise, `AMBIENT_FIG7C_STATUS.md` remains as she wrote it, and
this file is the only record of the above.
