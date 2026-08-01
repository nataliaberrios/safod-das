# SAFOD DAS sanity reproduction of Lellouch et al. 2019 Fig 7c

This directory is a clean, minimal pipeline that asks one question:

> Can the SAFOD DAS 2024–2025 dataset recover a body-wave Green's function from
> ambient noise cross-correlation when run with the Lellouch et al. 2019 recipe
> exactly as the paper describes it?

If the answer is yes, the existing pipeline can be repaired toward this canonical
configuration and the rest of the project (day/night, temporal monitoring, etc.)
can rest on a validated CC. If one day fails, the next controlled checks are:
stack more days, move the virtual source/top channel, and compare against the
2017 data used by Lellouch et al. before deciding that the 2024 data cannot
support this body-wave ambient-noise workflow.

## Target day

**2024-10-23** (Wednesday). Chosen because:
- The week 2024-10-22 to 2024-10-28 has zero M ≥ 1 NCEDC events within 30 km of
  Parkfield (lowest in your processed range).
- Wednesday → cultural noise is active.

Override with `SANITY_DATE=YYYY-MM-DD` if you want a different day. Use
`SANITY_DATES` for multi-day stacks:

```bash
export SANITY_DATES=2024-10-22:2024-10-28          # inclusive range
export SANITY_DATES=2024-10-22,2024-10-23,2024-10-24
```

By default, the CC stage uses **all continuous files** for the selected UTC
date(s), matching the paper's description more closely than the earlier
daytime-only sanity run. Use `SANITY_HOURS=day` or `SANITY_HOURS=night` only
when explicitly testing source timing.

## Faithful Lellouch parameters

| Parameter            | Lellouch 2019           | This pipeline           |
|----------------------|-------------------------|-------------------------|
| Bandpass             | 5–20 Hz                 | 5–20 Hz                 |
| Window               | 30 s, 50% overlap       | 30 s, 50% overlap       |
| Time normalization   | running-AM (Bensen 2007) | running-AM, **0.1 s window** |
| Spectral whitening   | not described in paper  | **OFF by default** (optional knob `SANITY_WHITEN=true` for phase-only) |
| Geometry (a) source  | top channel             | `SANITY_SOURCE_CH` (= `ch_start` by default) |
| Receivers            | along array             | all channels in `[ch_start, ch_end)` |
| Per-pair stacking    | R ± 10, pre-shifted at 3200 m/s | post-processed in `sanity_plot.py` |
| Data length          | 1 day                   | 1 day by default; multi-day via `SANITY_DATES`; all continuous files unless `SANITY_HOURS` is set |

The previous pipeline used `TN_WINDOW = 10.0` s (~100× too long for body-wave CC),
included event-triggered ~20 s files in the "continuous" stack, and never applied
the adjacent-channel pre-shift stack. Those are all corrected here.

## What's intentionally **not** the same as the paper

| Item               | Paper           | Here              | Why it's OK |
|--------------------|-----------------|-------------------|-------------|
| Sample rate        | 2500 Hz         | 500 Hz (after `Desample=20`) | Nyquist 250 Hz → 5–20 Hz CC band fully usable |
| Gauge length       | 10 m            | 16 m              | In 5–20 Hz, the gauge-length sinc term is close to flat for 3200 m/s body waves and 1500 m/s guided/tube waves; not expected to be a first-order obstacle |
| Year               | 2017            | 2024–2025         | Same fiber and well, but cultural noise patterns may have shifted |

## Files

- `preflight.py` — Diagnostic on 1 hour of daytime raw data. Produces RMS-per-channel,
  spectrogram, time-distance plot, and **F-K spectrum**. The F-K plot is the existence
  check: if there is no body-wave-velocity ridge in the wavefield, no CC trick will
  recover one.
- `sanity_cc.py` — Canonical Lellouch CC for one or more days. Saves npz files to
  `/oak/stanford/groups/ettore88/nberrios/sanity_v1/` with date, channel, and
  source-channel tags.
- `sanity_plot.py` — Loads the daily npz, applies R±10 pre-shifted stacking at
  3200 m/s, plots the Fig 7c equivalent and the F-K of the stacked CC.
- `run_sanity.sbatch` — SLURM submission. `STAGE` env var picks the script.

## Pushing this code to Sherlock

This `sanity/` directory was written locally inside your repo working tree but
has not been committed or pushed yet. From your laptop:

```bash
cd ~/Documents/Claude/Projects/DAS/safod-das
git status                                    # confirm sanity/ shows up as new
git add sanity/
git commit -m "Add canonical Lellouch sanity reproduction pipeline"
git push
```

Then on Sherlock:

```bash
ssh nberrios@login.sherlock.stanford.edu
cd /home/groups/ettore88/nberrios/safod-das
git pull
ls sanity/                                    # confirm files arrived
```

If `git pull` complains about local changes in the parent repo, stash or commit
them first — this `sanity/` directory is self-contained and won't conflict with
the existing `stack_daily.py` workflow.

## Running on Sherlock

The SLURM script `run_sanity.sbatch` matches the conventions of your existing
`run_daily_array.sbatch` (gpu partition, `cuda/11.7.1`, `py-pytorch/2.0.0_py39`,
the `das` conda env at `/home/users/nberrios/miniconda3/envs/das/bin/python`, and
the `DAS-utilities/python` directory on `PYTHONPATH`).

### Step 1 — Pre-flight diagnostic (fast, ~few minutes)

```bash
cd /home/groups/ettore88/nberrios/safod-das/sanity
sbatch --export=ALL,STAGE=preflight run_sanity.sbatch
```

Watch the job:

```bash
squeue -u nberrios
ls -lt logs | head
tail -f logs/safod_sanity_<JOBID>.out
```

When it finishes, the four PNGs land here on Sherlock:

```
/home/groups/ettore88/nberrios/safod-das/sanity/preflight_out/
  rms_per_channel.png
  spectrogram.png
  timedist_one_file.png
  fk_one_file.png             <-- THE existence check
```

Copy them back to your laptop to look at:

```bash
# from your laptop
mkdir -p ~/Documents/Claude/Projects/DAS/sanity_results
rsync -avz nberrios@login.sherlock.stanford.edu:/home/groups/ettore88/nberrios/safod-das/sanity/preflight_out/ \
  ~/Documents/Claude/Projects/DAS/sanity_results/preflight_out/
```

**Decision point.** Look at `fk_one_file.png`. The dashed reference lines show the
slope expected for body waves (3200 / 1600 m/s), tube waves (1500 m/s), and surface
modes (500 m/s) in (k, f) space. Body-wave energy shows up as a *steep* ridge near
the f-axis; tube/surface modes show up as *shallow* ridges near the k-axis.

- **Steep ridge inside the 5–20 Hz cyan band → proceed to Step 2.**
- **Energy concentrated along the 1500 m/s line or shallower → stop and talk to
  Ettore** before burning cluster time on CC. The wavefield doesn't have what we
  need, and no preprocessing trick will recover a body-wave Green's function.

### Step 2 — Canonical CC for one day (~45–90 min depending on file count)

```bash
sbatch --export=ALL,STAGE=cc run_sanity.sbatch
```

The output lands at a tagged path such as:

```
/oak/stanford/groups/ettore88/nberrios/sanity_v1/sanity_cc_2024-10-23_hoursall_ch150-800_src150.npz
```

This stays on `oak` (not in the repo) because `.npz` files are heavy and shouldn't
go through git.

### Step 3 — Post-process and plot

Can run on Sherlock or locally. On Sherlock:

```bash
sbatch --export=ALL,STAGE=plot run_sanity.sbatch
```

Plots land at:

```
/home/groups/ettore88/nberrios/safod-das/sanity/plot_out/
  cc_raw_<TAG>.png        # before pair-stacking
  cc_stacked_<TAG>.png    # Lellouch Fig 7c equivalent
  cc_fk_<TAG>.png         # F-K of stacked CC
```

Then rsync them back the same way as the preflight outputs.

To run plotting locally instead (faster iteration on visualization choices):

```bash
# from laptop, after rsyncing the npz back
rsync -avz nberrios@login.sherlock.stanford.edu:/oak/stanford/groups/ettore88/nberrios/sanity_v1/ \
  ~/Documents/Claude/Projects/DAS/sanity_results/sanity_v1/

cd ~/Documents/Claude/Projects/DAS/safod-das/sanity
SANITY_NPZ=~/Documents/Claude/Projects/DAS/sanity_results/sanity_v1/sanity_cc_2024-10-23_hoursall_ch150-800_src150.npz \
SANITY_PLOT_OUT=~/Documents/Claude/Projects/DAS/sanity_results/plot_out \
python sanity_plot.py
```

### Inspecting and debugging jobs

```bash
# what's queued / running
squeue -u nberrios

# wall time + memory used for completed jobs
sacct -j <JOBID> --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem

# logs (most recent)
ls -lt logs | head
cat logs/safod_sanity_<JOBID>.out
cat logs/safod_sanity_<JOBID>.err
```

If the job dies on import, paste the first ~50 lines of `.err` — most often it's
the PYTHONPATH or env activation mismatching.

### Re-running with overrides

```bash
# different day
export SANITY_DATE=2024-10-24
sbatch --export=ALL,STAGE=cc run_sanity.sbatch

# multi-day stack
export SANITY_DATES=2024-10-22:2024-10-28
export SANITY_HOURS=all
sbatch --export=ALL,STAGE=cc run_sanity.sbatch

# source-timing tests, only after the all-hours run
export SANITY_HOURS=day
sbatch --export=ALL,STAGE=cc run_sanity.sbatch
export SANITY_HOURS=night
sbatch --export=ALL,STAGE=cc run_sanity.sbatch

# tighter channel range after looking at the RMS plot
export SANITY_CH_START=180
export SANITY_CH_END=750
sbatch --export=ALL,STAGE=cc run_sanity.sbatch

# source-channel sweep inside the same receiver aperture
export SANITY_CH_START=150
export SANITY_CH_END=800
export SANITY_SOURCE_CH=180
sbatch --export=ALL,STAGE=cc run_sanity.sbatch

# test phase-only whitening on top of strict-Lellouch (paper does NOT specify this)
export SANITY_WHITEN=true
sbatch --export=ALL,STAGE=cc run_sanity.sbatch

# different pre-shift velocity for the adjacent-channel stack in plot stage
export SANITY_VAVG=2800
sbatch --export=ALL,STAGE=plot run_sanity.sbatch
```

## Reading the output

`sanity/plot_out/cc_stacked_<DATE>.png` is the Lellouch Fig 7c equivalent. A
successful body-wave Green's function recovery will show:

- A clear arrival following the dashed `+3200 m/s` reference line on the causal
  (positive-lag) side, sloping outward from the source channel.
- Substantially weaker zero-lag stripe and weaker arrivals in any other direction.

`cc_fk_<DATE>.png` cross-checks the result by showing what apparent velocities
are dominant in the recovered correlations.

## Override knobs

```bash
export SANITY_DATE=2024-10-24
export SANITY_CH_START=150       # uphole/in-well boundary (verify with rms_per_channel.png)
export SANITY_CH_END=800
export SANITY_SOURCE_CH=150
export SANITY_VAVG=3200.0        # pre-shift velocity for adjacent-channel stack
export SANITY_PAIR_HALF=10       # R±10 -> 21-channel stack
sbatch --export=ALL,STAGE=all run_sanity.sbatch
```

## Controlled follow-up tests before abandoning 2024 body-wave CC

1. **Stack more days using all continuous files.** Keep preprocessing fixed and
   increase the number of files without imposing a day/night source assumption:

   ```bash
   export SANITY_DATES=2024-10-22:2024-10-28
   export SANITY_HOURS=all
   export SANITY_CH_START=150
   export SANITY_CH_END=800
   export SANITY_SOURCE_CH=150
   sbatch --export=ALL,STAGE=cc run_sanity.sbatch
   sbatch --export=ALL,STAGE=plot run_sanity.sbatch
   ```

2. **Compare source timing.** After the all-hours result, explicitly test whether
   daytime cultural noise helps or hurts the stack:

   ```bash
   for hours in day night; do
     export SANITY_DATES=2024-10-22:2024-10-28
     export SANITY_HOURS=$hours
     export SANITY_SOURCE_CH=150
     sbatch --export=ALL,STAGE=cc run_sanity.sbatch
   done
   ```

3. **Move the virtual source/top channel.** Keep the same channel aperture, but
   test whether the chosen source channel is suppressing the result:

   ```bash
   for src in 150 180 200 250; do
     export SANITY_SOURCE_CH=$src
     export SANITY_HOURS=all
     sbatch --export=ALL,STAGE=cc run_sanity.sbatch
   done
   ```

   Plot each output by passing `SANITY_NPZ=/oak/.../<tag>.npz` or by letting
   `sanity_plot.py` pick the most recent `.npz`.

4. **Run the same code on 2017 data.** Use the Lellouch-era data/manifest path
   through `SAFOD_CSV` and, if needed, a matching `SANITY_CH_START`,
   `SANITY_CH_END`, and `SANITY_SOURCE_CH`. This separates implementation
   reproducibility from 2017-vs-2024 wavefield differences.
