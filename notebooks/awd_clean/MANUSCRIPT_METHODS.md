# Manuscript methods — whole-paper draft

Growing draft following the writing sequence in the project plan (§21: Methods
first). The Deep-specific subsections live in
[`DEEP_DVV_METHODS_DRAFT.md`](DEEP_DVV_METHODS_DRAFT.md) and are cross-referenced
rather than duplicated here, so there is one source per number.

Every value below is transcribed from `awd_manifest.csv`,
`canonical_epoch_stacks_paired_deep_all.npz`, or a named script. Anything I could
not source is marked **[NEEDS SOURCE]** rather than guessed.

Status: §2 Experiment and instrumentation drafted. Remaining Methods subsections
listed at the end.

---

## 2. Experiment and instrumentation

### 2.1 Site and fibres

Two distinct borehole DAS installations at SAFOD recorded the same repeated
surface source. They are not interchangeable observations of one wavefield: the
installations differ in coupling and geometry, and each recovers a different
coherent mode.

The **Nano** fibre is cemented and shallow. It provided 732 channels at 1.266 m
channel spacing, an aperture of 927 m. Data were written as protobuf
(`timeseriesScaler_*.pb`).

The **Deep** fibre is wireline-deployed and reverses at channel 1702, so the
3200 recorded channels of 2.0419 m spacing trace an outbound limb of 3475 m and a
return limb of 3059 m along 6532 m of fibre. Acquisition parameters are recorded
in the filenames (`SAFOD-Deep-10mGL-1000HzFs-2mChDualPulse_*.h5`): **10 m gauge
length, 1000 Hz sample rate, 2 m channel spacing, dual-pulse interrogation.**

Both systems sampled at 1000 Hz.

**[NEEDS SOURCE]** interrogator make and model for each system; Nano gauge length
and pulse configuration; the physical depth of the Deep turnaround; casing and
completion details. Channel coordinate is distance along fibre throughout — the
coordinate-to-depth mapping remains provisional and no result here depends on it.

### 2.2 Source and survey

The source was an accelerated weight drop (AWD) held at a fixed position for the
survey, approximately 15 m horizontally from the Nano wellhead. Drops were fired
in bursts; drops separated by more than 60 s are assigned to different bursts
(`build_manifest.py`, `BURST_GAP_S = 60`).

The survey comprised **49 bursts at a median cadence of 30.0 min (range
27.6–31.2 min), spanning 23.96 h** from 2026-06-16 23:47 to 2026-06-17 23:44 UTC.
Note that 24 h at half-hour cadence gives 48 *intervals* and therefore 49
*bursts*; there are no gaps exceeding 1.5× the median cadence.

Each drop is windowed 0.5 s before to 3.0 s after its GPS time.

### 2.3 Drop inventory and coverage

Four different drop counts are legitimate at different stages, and they must not
be used interchangeably:

| Stage | Nano | Deep | Both |
|---|---:|---:|---:|
| GPS drops listed in `p26.cc9.txt` | — | — | 989 |
| Rows in `awd_manifest.csv` | 988 | 926 | — |
| Full 0.5/3.0 s window fits inside the data file | 970 | 875 | — |
| **Used by the burst stacks** | 970 | 875 | **859 common** |

The **analysed dataset is 859 drops common to both fibres, across 46 bursts.**

Two attritions account for the difference, neither of them quality rejection.

**The Deep fibre stopped recording after burst 45 (22:16:23 UTC).** Bursts 46, 47
and 48 — at 22:46, 23:16 and 23:44 UTC on 2026-06-17 — carry full Nano coverage
and zero Deep coverage. This alone is why the analysis uses 46 bursts and not 49.
Consequently, **restricting the analysis to drops recorded on both fibres costs
the Deep analysis nothing**: a Deep-only analysis could not have used those three
bursts either. The 46-burst count is set by Deep's own coverage, not by the
pairing requirement.

The 988→970 and 926→875 reductions are **window truncation at file boundaries**:
a drop is retained only when its full pre- and post-window falls inside the
recording file (`paired_stack_job_deep_all.py:239`). Deep loses proportionally
more because its 3.0 s post-window is long relative to the `.h5` segmentation.
The inventory closes: 859 common + 111 Nano-only = 970, and 859 + 16 Deep-only =
875.

Because the missing bursts fall at the end of the survey, the analysed series is
truncated at one end to 21.98 h (held-out bursts) or 22.49 h (all bursts) rather
than the nominal 24 h. This matters only for time-series fits over the survey
duration, where it slightly worsens the degeneracy between a diurnal signal and
instrumental drift.

### 2.4 Burst stacking

Drops within a burst are stacked to form one waveform per burst per channel.
Bursts are never silently truncated to an expected drop count; the retained count
per burst is carried through as a stacking weight. Individual drops vary
substantially in timing, amplitude, waveform similarity and signal-to-noise
ratio, whereas burst stacks are highly repeatable — the quantitative repeatability
hierarchy is given in §3.5.

---

## Remaining Methods subsections

Following the plan's §21 order. Sourcing status noted so nothing inherits
unearned confidence.

| Subsection | Source | Status |
|---|---|---|
| 3.1 Preprocessing and coordinate conventions | `PREPROCESSING.md` | not drafted |
| 3.2 Nano mode identification | `nano_mode_identification.py` | not drafted, **not verified by me** |
| 3.3 Deep guided-mode identification and split-sample validation | `deep_tube_validation.py`, `DEEP_DVV_METHODS_DRAFT.md` | drafted in the Deep file |
| 3.4 Moveout correction and beam construction | `DEEP_DVV_METHODS_DRAFT.md` | drafted in the Deep file |
| 3.5 Repeatability metrics | Nano repeatability analysis | not drafted, **numbers not verified by me** |
| 3.6 Nano injection–recovery | `nano_dvv_injection_recovery.py` | not drafted; parameters verified |
| 3.7 Deep local-delay-gradient estimator | `DEEP_DVV_METHODS_DRAFT.md` | drafted in the Deep file |
| 3.8 Blinding, preregistration, software audit | `DEEP_DVV_PREREGISTRATION.md`, `DEEP_DVV_STATUS.md` §7 | drafted in the Deep file |
| 3.9 Reliability definitions | `DEEP_DVV_STATUS.md` §2 | drafted in the Deep file |
| 3.10 Influence diagnostics | `deep_dvv_influence.py` | drafted in the Deep file |
| 3.11 Paired-leg estimator | `deep_dvv_paired_legs.py` | drafted in the Deep file |
| 3.12 Predefined Nano–Deep comparison | `DEEP_DVV_STATUS.md` §2, §3 | drafted in the Deep file |
| Small-signal tidal benchmark | `deep_dvv_tidal_fit.py`, `safod_tides.ipynb` | **placement undecided** — see `DEEP_DVV_STATUS.md` §5 |

The Deep half of Methods is therefore essentially complete. The gaps are the
**Nano** subsections (3.2, 3.5, 3.6) and preprocessing (3.1), plus a decision on
where the tidal benchmark lives.
