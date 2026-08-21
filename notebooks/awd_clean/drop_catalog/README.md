# AWD weight-drop catalogue — node timing → DAS detection

Clean, self-contained products for the June 2026 SAFOD active-source (AWD)
survey. Built 2026-08-20. Everything here is regenerable with one command:

```bash
conda activate das
python notebooks/awd_clean/drop_catalog/build_drop_catalog.py
```

Survey window: **2026-06-16 23:47 → 2026-06-17 23:47 UTC** (24.0 h),
local **2026-06-16 16:47 → 2026-06-17 16:47 PDT**.

---

## Provenance — read this first

**The drop times were delivered, not derived here.** The cross-correlation picks
came as `p26.cc9.txt` (node 453009664) and `p26.cc4.txt` (node 453001432): the
same 989 weight drops picked independently on two standalone seismic nodes
~10–20 m apart. No script in this repository produces those files — every
reference reads them. `p26.cc9.txt` first appears on disk 2026-06-29; the rest of
the `Check shots/` delivery was copied in 2026-08-02.

**What is this project's work** is everything downstream of that:

| Step | Script | Product |
|---|---|---|
| 1. Validate the delivered timing against a second node | `faultzone/timing_files_compare.py` | timing uncertainty |
| 2. Intersect drop times with DAS file coverage | `awd_clean/build_manifest.py` | `awd_manifest.csv` |
| 3. Test whether each drop is *detectable* on Nano | `awd_clean/nano_hierarchical_repeatability.py` | `nano_drop_repeatability.csv` |
| 4. Test whether a drop is *visible at all* on Deep | `awd_clean/deep_drop_visibility.py` | `deep_drop_visibility.{txt,png}` |
| 5. Join into one catalogue | `drop_catalog/build_drop_catalog.py` | this directory |

The distinction step 3–4 enforces, in `deep_drop_visibility.py`'s own words:

> **Coverage is not detection.**

A drop being inside a DAS file only means the recorder was running. Whether the
drop is *in* the data is a separate measurement.

---

## Products

| File | Rows | Contents |
|---|---|---|
| `awd_drop_catalog.csv` | 989 | one row per delivered pick: UTC + local time, burst/drop id, node CC, two-node offset, Nano/Deep coverage, Nano detection metrics and flag |
| `awd_burst_summary.csv` | 49 | per burst: start UTC + local, duration, drop count, coverage counts, detected count, median beam SNR |
| `timing_uncertainty.txt` | — | the two-node comparison in full |

## The notebook

`AWD_drop_catalog.ipynb` is the readable version of all of the above: the
numbers, then three figures. It is **committed with its outputs** so it can be
read without running anything — the same exception `manuscript/README.md` makes,
and for the same reason. It is 0.4 MB, not a dashboard.

It has **one switch**, in the first code cell:

| `REBUILD_PRODUCTS` | What happens |
|---|---|
| `False` *(default)* | load the committed CSVs and plot. Seconds. |
| `True` | re-run `build_drop_catalog.py` from the node picks first, then plot. Needs `Check shots/` on disk. |

The figures are identical either way; that is what the switch is for.

The notebook is **generated** — edit `build_drop_notebook.py`, never the
`.ipynb`. To rebuild and re-execute with figures embedded:

```bash
sbatch awd_clean/drop_catalog/exec_drop_nb_job.sh
```

That job does both halves, because they need different interpreters: the build
step needs the **system `python3`** (it has `nbformat`), and execution needs the
**`das` kernel** (which does not). Execution goes through
`manuscript/execute_notebook.py`, which exits non-zero if a cell that draws a
figure captured none — `run_nb_cells.py` cannot catch that.

Figures: (1) every drop across the 24 h with Nano detection status; (2) six
drops from six different bursts overlaid, with the delivered drop time marked by
a red dotted line; (3) one full burst as an image, same marker. Figures 2 and 3
read from `../nano_hierarchical_repeatability.npz`, which is gitignored — they
need that file present on Sherlock.

---

## Headline numbers

**Timing (delivered picks, validated here).** The two nodes agree to a median
**−1.00 ms**, MAD **0.50 ms**, and 99.5 % of drops agree within 10 ms. Picks are
quantised to 0.5 ms.

Both nodes sit metres from the source, so they share the same unknown
source-to-node travel time. This is *relative* pick quality, **not** the absolute
offset — that is the ~90 ms measured against the check shot
(`faultzone/digitize_checkshot.py`).

**Coverage (988 of 989 drops enter the manifest; one falls outside DAS coverage).**

| | count | share |
|---|---:|---:|
| bursts | 49 | 19–23 drops each |
| Nano coverage | 988 | 100.0 % |
| Deep coverage | 926 | 93.7 % |
| both fibres | 926 | 93.7 % |

**Detection — the part that matters, and the two fibres differ sharply.**

*Nano:* **409 of 988** drops meet both criteria (leave-one-out signal NCC > 0.90
**and** beam SNR > 10 dB). Median LOO signal NCC **0.889** against a noise-arm NCC
of **0.288**; median beam SNR **7.7 dB**, 90th percentile **25.7 dB**. The
leave-one-out delay is median 0.000 s with a 10–90 range of ±1 ms — consistent
with the node timing above.

*Deep:* largely **not** detected. Stacked drop/sham peak-SNR ratio has median
**1.017** — indistinguishable from no signal across the fibre as a whole. 895 of
3200 channels exceed ratio 1.5, in a few localised sections.

> **Corrected 2026-08-20.** An earlier version of this file quoted the strongest
> responding section as fibre **200–664 m, ratio 34.2** with no further comment,
> and stated that no depth registration exists for Deep. Both need fixing.
>
> Registration **does** exist: `safod_geometry.py`, from
> `SAFOD_Phase2_GeoReferenced_Channels.xlsx`. It puts the surface lead-in at
> channels **0–210**, the first in-hole channel at **211**, and the reversal at
> channel **1700** (MD 3053 m, TVD 2549 m), after which the fibre retraces the
> same depths on a return limb.
>
> At dx = 2.0419 m, fibre 200–664 m is channels **98–325** — so roughly the first
> half of that section is **surface lead-in, fibre that never enters the hole**.
> Its 34.2 ratio is therefore not straightforwardly a formation response, and
> should not be quoted as one. The second section, fibre 706–1060 m (channels
> 346–519), is genuinely in-hole and near-vertical.
>
> `deep_drop_visibility.py` has been re-run with the registration applied; it now
> reports TVD and limb per section and excludes the lead-in from its in-hole
> summary. See `deep_drop_visibility.txt` for the current numbers, which
> supersede any Deep figure quoted here.

---

## Caveats

- The **detection thresholds** (NCC > 0.90, SNR > 10 dB) are set in
  `build_drop_catalog.py` for reproducibility. They were **not** pre-registered,
  so treat the 409 count as descriptive, not as a tested result.
- **78.8 %** of node picks have Max_CC < 0.9. Median CC is 0.860 on node 9 and
  0.741 on node 4. The picks are usable but not uniformly high quality.
- 5 drops disagree between nodes by ≥ 10 ms, 4 by ≥ 100 ms. Those rows are in the
  catalogue and are worth excluding from any timing-critical stack.
- The Deep numbers come from a 30-drop stack against a 29-window sham arm, not
  all 926 covered drops.

## Raw inputs are not committed

`Check shots/` is gitignored (`awd_clean/.gitignore:3`), so `p26.cc9.txt` and
`p26.cc4.txt` are **not** in the repository. They are the provenance for every
number above. On Sherlock:

```
notebooks/awd_clean/Check shots/p26.cc9.txt     (also copied to notebooks/p26.cc9.txt)
notebooks/awd_clean/Check shots/p26.cc4.txt
```

Whether they can be committed is a question for whoever supplied them.
