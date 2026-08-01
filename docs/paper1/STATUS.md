# Paper 1 — status and roadmap

_Single source of truth. Supersedes any loose status file; where numbers
disagree, this file wins._

- Consolidated 2026-07-21 (author).
- Audited against the live notebook 2026-07-31 — see
  [`notebook_audit_2026-07-31.md`](notebook_audit_2026-07-31.md). Audit findings
  are folded in below and marked **[audit]**.

Live analysis notebook: `notebooks/JULY24_fml_coupling.ipynb`.

---

## The paper in one sentence

Repeated active-source (AWD) borehole DAS at SAFOD produces coherent wavefields
over a finite usable aperture, including a candidate borehole-guided arrival,
but cannot resolve De Fazio-scale (5×10⁻⁴) tidal velocity change — we quantify
that observability limit.

**Venue:** SRL technical note. Not BSSA — there is no new Earth inference; this
is measurement-system characterization plus a candidate mode plus a formal null.

**Scope:** Version A (wireline-only observability). Version B (dual-fiber
coupling) is stronger but blocked on fiber-identity confirmation; parked.

---

## Results — done

| Result | Value | Confidence |
|---|---|---|
| Shallow aperture repeatability cliff | CC 0.987→0.636, SNR 22.5→4.8 dB, lag scatter 0.43→1.2 ms across ~374 m | Clean, defensible |
| Pre-registered tidal null | amp 3.4×10⁻⁴ ± 2.3×10⁻⁴, p_null = 0.51, **95% UL = 7.98×10⁻⁴ (1.6× De Fazio)** | Best-executed piece; robust to dropping first 6 h |
| dv/v floor (depth-median) | 2.87×10⁻³ (5.7×) | LOO reference, amplitude-balanced windows |
| dv/v floor (common-mode removed) | **1.54×10⁻³ (3.1×)** | Headline sensitivity |
| Sensitivity lineage | 65× → 5.7× → 3.1× → 1.6× | All methodology, same data |
| Candidate guided arrival | ~1.4–1.5 km/s, 3 bands, semblance 3.8–7.3× null95 | Keep as "candidate", not definitive tube wave |
| Source is the limiter | σ_α ≈ 0.30 ms common-mode; 39% amplitude CV (does not couple to timing); compaction drift matches observed plate settlement | Cross-validated vs Otway/SOV literature |
| Deep-fiber within-burst repeatability | CC 0.054, NRMS 137%, lag scatter 10.9 ms | Clean negative; deep fiber unusable per-shot |

**Permutation screen on the guided mode — now closed [audit].** `TW_N_NULL` and
`N_CONFIRM_NULLS` are both already set to **199** (cells 158, 162). Observed
semblance 0.326 vs shuffled-depth 99th percentile 0.084, `p_perm = 0.005`. That
is the 1/200 floor for 199 nulls, not the old 1/21. This was previously listed as
outstanding work; it is complete.

### Open mechanism question

α–ε residual correlation r = +0.59; split-aperture halves +0.73 / +0.77 but
cross-half −0.16. **Unresolved** — likely in-window contamination. If so,
1.54×10⁻³ may not be a hard floor. This is the one scientific loose end that
could change a headline number.

---

## Known defects in the notebook [audit]

Ordered by effort-to-value. Only item 2 can move a reported number.

1. **Dead line breaks the headline-null robustness figure.** Cell 175 — the
   trend-sensitivity table and null-vs-data figure — raises a `ValueError` on:

   ```python
   dfz = c_tr[0] + (Glin[:,1:]@c_tr[1:] - c_tr[0]) * 0   # baseline
   ```

   `dfz` is assigned once, multiplied by zero, and never read. Delete the line.
   It is the only error in the notebook, and it sits on the last executed cell,
   at exactly the result the paper leads with.

2. **`t_hours` is one element short.** Cell 155 prints `epoch_stacks has 49
   epochs; contaminated has 49; t_hours has 48`. Anything pairing `t_hours`
   against a 49-length series truncates an epoch or misaligns the time axis by
   one half-hour slot. **Fix before rebuilding the locked table.**

3. **Cell 132 still present** — the `SEP_ZMAX_M = 130.0` short-lever OLS fit
   (σ_α = 4.57 ms). Already fenced with `# DEMO / DO NOT QUOTE`, so it is not
   circulating as a competing result, but deleting it closes the two-σ_α
   ambiguity for good. The other cell named in the original blocker (the
   quarantined ch-374 zero-scatter fit) is **already gone** — absent from both
   `JULY24` and `JULY23`.

4. **Stale prose.** The summary (cell 173) says "48 half-hour epochs" and "46 of
   48 clean epochs". The analysis ran on **49**, with **47 of 49** clean. The
   `48`s elsewhere are a superseded file-count estimate (cell 15, self-labelled
   provisional) and cosmetic 0-indexed display (`epoch 00/48 … 48/48`) — both
   harmless.

5. **Internal contradiction on depth.** Cell 115 says the AWD signal "does not
   penetrate below ~100m", citing CC 0.027 at 633 m. That evidence supports "no
   usable signal at 633 m", not a 100 m limit, and it contradicts the 373–374 m
   usable interval established in cells 39, 58, and 173. Restate.

---

## Work remaining

### Analysis cleanup — gates everything

1. Apply defects 1–5 above. Items 1, 3, 4, 5 are edits, not re-runs.
2. **[BLOCKER] Clean-kernel re-run on 49 epochs**, then rebuild
   `paper1_locked_results.md` from that single run. The notebook currently holds
   three apertures and two σ_α values, so no number is provably from one dataset.
3. **Pick one aperture and delete the others.** Candidates: 127–374 m (what the
   dv/v chain runs on), 76–373 m (four-criteria with the `SIGMA_T_MAX` gate),
   74.7–492.5 m. Recommend the four-criteria aperture, with the dv/v chain
   brought into agreement. **This is the largest genuine analysis task left.**

### Writing — only Methods exists

4. Write Results (4 subsections → figure set), Intro, Discussion, Conclusions,
   Data & Resources.
5. **Read Niu et al. (2008, Nature)** — same site, pre-DAS active-source velocity
   monitoring, currently missing from the prior-art memo. Highest-value citation;
   get it in before writing the Intro.
6. Assemble ~5 figures against the reconciled numbers.

### External blocker

7. **Confirm fiber identity with Ettore.** Cell 149 labels a block "WIRELINE
   (DEEP)" and then says "compare wireline 373 m" — the wireline/cemented labels
   are tangled. Blocks Version B and the use of both words throughout.

---

## Critical path

```
defects 1-5 (edits)  →  clean 49-epoch re-run  →  reconcile aperture
  →  rebuild locked table  →  read Niu (2008)  →  write Results/Intro/Discussion
  →  figures  →  submit
```

**Fastest next action:** delete the dead line in cell 175 and fix `t_hours`.
That restores the headline-null figure and removes the only defect capable of
changing a number, in well under an hour.

---

## Companion documents

These live outside the repo. Bring them in under `docs/paper1/` so the writing
and the analysis version together:

- `paper1_methods_draft.md` — Methods §2, drafted
- `paper1_skeleton_SRL.md` — SRL skeleton, A/B decision, figure set
- `paper1_prior_art_memo.md` — prior-art and novelty positioning
- `paper1_locked_results.md` — to be rebuilt after the clean re-run
