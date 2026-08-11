# What code produced the paper

Every number in the manuscript, traced to the script that made it and the order
things must run in. If you can only read one thing to answer "what did we
actually do", read this.

Environment for all of it:

```bash
ml gcc/12.4.0
source ~/miniconda3/etc/profile.d/conda.sh
conda activate das
cd /home/groups/ettore88/nberrios/safod_das_git/notebooks/awd_clean
```

Nothing below needs a GPU. Anything heavy has an `.sbatch` beside it; submit with
`-p serc`.

---

## The pipeline, in dependency order

Stages 1–2 build the data. Stage 3 is the Nano analysis, stage 4 the Deep. They
are independent of each other once stage 2 exists.

### 1. Manifest — which drops exist, on which fibre

```bash
python build_manifest.py
```

| | |
|---|---|
| Reads | `p26.cc9.txt` (989 GPS drop times), Nano `.pb` and Deep `.h5` coverage on `$OAK` |
| Writes | `awd_manifest.csv` — 988 rows, one per drop |
| Key setting | `BURST_GAP_S = 60` — drops more than 60 s apart start a new burst |
| Gives the paper | 49 bursts; Nano available on 988 drops, Deep on 926; Deep stops after burst 45 |

### 2. Burst stacks — the file everything else reads

```bash
sbatch rebuild_stacks.sbatch      # runs ../paired_stack_job_deep_all.py
```

| | |
|---|---|
| Reads | `awd_manifest.csv`, raw `.pb` and `.h5` on `$OAK` |
| Writes | **`canonical_epoch_stacks_paired_deep_all.npz`** (2.7 GB) |
| Contents | `nano_stacks` (49, 732, 3500), `deep_stacks` (49, 3200, 3500), `n_common`, `begtimes_str`, `fs`, `dx_nano`, `dx_deep` |
| Key line | `paired_stack_job_deep_all.py:239` — a drop is kept only if its full 0.5 s / 3.0 s window fits inside the file. This is the 988→970 Nano and 926→875 Deep attrition |
| Gives the paper | 859 drops common to both fibres, across 46 bursts |

**Everything downstream reads this one file.** If it is rebuilt, every result
below must be rerun.

### 3. Nano analysis

```bash
python nano_mode_identification.py                 # minutes
sbatch nano_hierarchical_repeatability.sbatch      # reads raw .pb again, not the stacks
sbatch nano_dvv_injection_recovery.sbatch          # or: python nano_dvv_injection_recovery.py --stage all
```

| Script | Writes | Gives the paper |
|---|---|---|
| `nano_mode_identification.py` | `nano_mode_identification.{csv,txt,npz,png}` | apparent speeds by band (2900–3300 m/s), positive ordering P = 1.000, dispersion trend 0.483 µs m⁻¹ Hz⁻¹ |
| `nano_hierarchical_repeatability.py` | `nano_drop_repeatability.csv`, `nano_stack_convergence.csv`, `nano_hierarchical_repeatability.{npz,png,txt}` | signal/noise NCC 0.889 / 0.288, across-burst 0.976, substack convergence, beam SNR spread |
| `nano_dvv_injection_recovery.py` | `nano_dvv_{blind_trials.npz, blind_truth.csv, recovery.csv, summary.csv, injection_recovery.npz/.png/.txt}` | Nano null threshold 3.25×10⁻³, reliable level 1×10⁻², lever arm 0.121 s |

Note `nano_hierarchical_repeatability.py` re-reads the raw `.pb` files rather
than the stacks, and uses **988 drops across all 49 bursts** — not the 46-burst
common set. That is deliberate and is why its counts differ from everything else.

### 4. Deep analysis

Run in this order; each depends on the one before.

```bash
python deep_tube_validation.py                              # ~10 min, establishes the mode
sbatch deep_dvv_injection_recovery.sbatch                   # ~16 min, 96 GB — the main result
python deep_dvv_synthetic_validation.py                     # seconds
sbatch <wrap> deep_dvv_influence.py                         # ~2 min
python deep_dvv_paired_legs.py                              # seconds
python deep_dvv_paired_legs.py --population allbursts       # seconds
python deep_dvv_tidal_fit.py                                # ~1 min
python deep_dvv_tidal_fit.py --population allbursts         # ~1 min
```

| Script | Writes | Gives the paper |
|---|---|---|
| `deep_tube_validation.py` | `deep_tube_{candidates.csv, burst_repeatability.csv, validation.npz/.txt}`, 3 PNGs | the mode exists: permutation p = 0.002 in all four leg/band tests, 100% of bursts positive-slowness |
| `deep_dvv_injection_recovery.py` | `deep_dvv_frozen_trajectory.json`, `deep_dvv_blind_truth.csv`, `deep_dvv_{recovery,summary,controls,nano_comparison}.csv`, `.npz`, `.png`, `.txt`, plus ~4.5 GB of blinded gathers in `$SCRATCH/deep_dvv_blind/` | frozen trajectories 1544.6 / 1549.7 m/s; thresholds 1.84×10⁻³ and 3.05×10⁻³; reliable levels 5×10⁻³ and 1×10⁻²; timing controls |
| `deep_dvv_synthetic_validation.py` | `deep_dvv_synthetic_validation.{csv,txt}` | injection is faithful: physical and pipeline paths agree to 5 sig figs, both legs pass |
| `deep_dvv_influence.py` | `deep_dvv_influence.{csv,txt}` | no single aperture or burst carries the result |
| `deep_dvv_paired_legs.py` | `deep_dvv_paired_{legs.csv, comparison.csv, legs.txt}` (+ `_allbursts`) | ρ = 0.121; paired scatter 7.9–8.8×10⁻⁴; detection level unmoved at 5×10⁻³ |
| `deep_dvv_tidal_fit.py` | `deep_dvv_tidal_fit.{csv,txt,png}` (+ `_allbursts`) | tidal 95% upper limit 1.03×10⁻³ at 46 bursts, clean nulls |

The four stages of `deep_dvv_injection_recovery.py` can be run separately —
`--stage freeze|inject|recover|summarize` — which is how the blinding is
enforced: `recover` never reads the truth table.

### 5. Tide theory (context, not a result)

`safod_tides.ipynb` in `notebooks/` — degree-2 forcing, the Niu expected-response
scale, and the Model A / Model B comparison. Smoke-test with
`python run_nb_cells.py safod_tides.ipynb` from `notebooks/`.

---

## If you had to explain this to someone in two minutes

> We dropped a weight repeatedly for 24 hours and recorded it on two fibres in
> the same borehole — one cemented, one on wireline. We first worked out what
> coherent arrival each fibre sees, and showed each arrival is real by splitting
> the bursts in half, picking the arrival on one half, and confirming it on the
> other against randomised controls.
>
> Then we asked how small a velocity change each fibre could detect. Rather than
> estimate that from theory, we took the real data, injected velocity changes of
> known size into it, and tried to recover them blind — the recovery code never
> saw the true answer. The smallest change recovered reliably is that fibre's
> sensitivity.
>
> The cemented fibre resolves 1%. The deep fibre's outbound branch resolves 0.5%,
> its return branch 1%. The interesting part is why the deep fibre isn't far
> better: its arrival travels 11.8× longer, which should help a lot, but its
> timing wobbles 7.6× more from burst to burst, and that cancels most of the
> advantage.
>
> Finally we fitted an Earth-tide model to the measurements. No tidal signal, and
> the upper limit is about 18× above what a tide should produce here. Not because
> the instrument is bad — because 24 hours isn't long enough to tell a daily
> signal apart from instrument drift.

---

## Things that will trip you up

- **`canonical_epoch_stacks_paired_deep_all.npz` is 2.7 GB** and is gitignored.
  It is not in the repo; it must be rebuilt or copied.
- **The blinded Deep gathers live on `$SCRATCH`**, which purges after 90 days of
  inactivity. They are regenerable by rerunning `--stage inject`.
- **Two burst-parity conventions exist.** `deep_tube_validation.py` uses
  `epoch % 2 == 1` for discovery; `deep_target_burst_repeatability.py` uses the
  opposite and calls it by the same name. The Deep dv/v work follows
  `deep_tube_validation.py`. Check before reusing either.
- **`cc_tools.py` defines several functions twice.** The second definition wins.
  Not used by this paper, but it is in the tree.
- **Numbers that look wrong but aren't:** 49 bursts vs 46 analysed (Deep stopped
  after burst 45); 988 vs 970 vs 859 drops (window truncation, then common-fibre);
  988 drops / 49 bursts in the Nano repeatability only.
