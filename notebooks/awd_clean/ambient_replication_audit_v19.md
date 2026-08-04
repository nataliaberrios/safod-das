# Ambient-noise interferometry replication audit (v19)

## Purpose

This is a diagnosis of the earlier attempt to reproduce Lellouch et al. (2019), not a proposal to repeat the same unconstrained cross-correlation experiment. The published result is specific: a one-day ambient-field interferometry stack produced a downgoing P-wave packet near 3.2 km/s and a depth-dependent travel-time pick. The authors’ public release contains 14 reference correlograms (`050m.npy` through `700m.npy`), each 1,901 samples long, described as stacks of seven one-day correlations; it does not provide the original continuous ambient records, so the exact 2017 computation cannot be rerun from that repository alone. The paper used two geometries: a fixed virtual source at the top of the array with receivers every 50 m, and constant 50 m source–receiver separations. The reported cross-correlations were filtered at 5–20 Hz, normalized to suppress earthquake transients, and evaluated from seven separate one-day stacks. See the [published paper](https://doi.org/10.1029/2019JB017533) and the authors' [released SAFODDAS products](https://github.com/ariellellouch/SAFODDAS).

## What the previous notebook actually did

The saved notebook `06_updated_cc_safod_nat.ipynb` provides a useful diagnostic record:

| Item | Previous implementation | Why it matters |
|---|---|---|
| Virtual-source geometry | Original channel 200 correlated with every channel in 150:800 | This is neither the paper's top-channel geometry nor its 50 m constant-offset geometry. A source at channel 200 also makes the dominant zero-lag autocorrelation sit inside the displayed array. |
| Spatial sampling | `dx = 1.0` m hard-coded | The archive index reports `dCh = 1.020952` m. The error is small for a 3.2 km/s diagnostic, but it should not be hidden in a velocity estimate. |
| Array/acquisition | 2024--2025 archive: 900 channels, native index `fs=10000` Hz, `dCh=1.020952` m, median gauge length 16.335 m; the processing call reports `fs=500` Hz | The released Lellouch data are a different acquisition: 800 channels, 250 Hz, approximately 1 m spacing. Reproduction requires the same 2017 records or an explicit instrument-transfer comparison. |
| Stacking design | 70 selected short files for day and 70 for night | The paper's result is a one-day continuous stack, repeated over seven different days for uncertainty. A random subset of short files is not equivalent unless its time coverage and weighting are demonstrated. |
| Preprocessing | Median-over-channel removal, 5--20 Hz zero-phase bandpass, running-absolute-mean normalization; spectral whitening was available but the saved stack used `None` | The running-absolute-mean step is consistent with the paper, but the unwhitened branch and median removal need to be treated as explicit alternatives, not as the paper's unique recipe. |
| Correlation scaling | FFT correlation divided by `npts`, with no pairwise RMS normalization | The saved arrays reach 1.55 in magnitude. Those values are not conventional bounded correlation coefficients, so amplitude/SNR comparisons with published correlograms are invalid. Timing can still be tested, but only after a normalized implementation is added. |
| QC and failure handling | The stacking loop catches exceptions and continues; the run reports 70/70 processed | This is not the primary failure here, but a replication result needs a per-file audit and explicit rejection reasons rather than a single count. |
| Geometry of the claimed moveout | F--K displays and overlays used an arbitrary 3.2 km/s reference | A reference line is not a recovered arrival. The paper picked the three-sample correlation maximum in 50 m blocks and compared seven day stacks. |

## Most likely reasons for the non-reproduction

1. **Different data set and installation.** The strongest explanation is that the 2024--2025 archive is not the 2017 800-channel, 250-Hz data set used by Lellouch et al. A failure to reproduce a waveform from a different interrogator, gauge length, channel count, and deployment is not evidence that the method fails.

2. **Wrong virtual-source geometry.** The previous source at channel 200 leaves the source autocorrelation in the middle of the profile. Lellouch's diagnostic geometry starts at the top channel or holds the separation at 50 m. This changes the expected travel-time intercept, the amount of near-zero-lag contamination, and the visual appearance of the moveout.

3. **Insufficiently matched stacking.** The prior day/night products are stacks of selected short files, while the published figure is explicitly a one-day ambient stack and the uncertainty is based on seven days. Ambient source direction, transient contamination, and file weighting can change the result even with identical filtering.

4. **Correlation normalization is not yet replication-safe.** The current FFT routine computes an unnormalized cross-correlation. Its output exceeding one is a direct diagnostic that it should not be called a normalized correlation coefficient. This can obscure a weak coherent packet beneath source-channel energy and prevents a like-for-like SNR comparison.

5. **Path and metadata ambiguity.** The CSV stores file paths without the `/data/SAFOD/` component, whereas the files currently exist under `/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/`. The previous notebook nevertheless processed 70/70 files, so this is not proven to be the cause of the negative result; it is a reproducibility hazard that must be resolved in the audit script.

6. **The test was aimed at the wrong observable.** A broad F--K panel from one virtual source is not the published ambient-interferometry observable. The decisive test is a causal/anti-causal correlogram for the paper's two geometries, with a 3.2 km/s moveout pick and day-to-day uncertainty.

## What is already ruled out

The prior work did not simply fail to make a pretty figure. It established that the saved source-200 stacks are dominated by near-zero-lag source autocorrelation and do not, by themselves, show the published 3.2 km/s travel-time ridge at useful amplitude. That is a valid result about that processing branch. It is not yet a controlled reproduction test because the acquisition, geometry, stacking, and correlation normalization differ from the published experiment.

## Controlled replication test

The minimum fair test is now fixed before looking at results:

1. Use the exact 2017 records/released cross-correlated products if they are available; otherwise label the experiment a transfer test on the 2024--2025 archive.
2. Correct every file path through an explicit existence check and save a per-file manifest.
3. Implement both published geometries: (a) top-channel source to receivers sampled every 50 m, and (b) 50 m source–receiver separation sliding down the array.
4. Compute causal and anti-causal correlations with pairwise RMS normalization, then repeat the same calculation with one-bit/phase-only normalization and with the current unnormalized scaling as a diagnostic.
5. Stack a complete UTC day, reject or flag gaps and transients, and repeat for seven independent days where coverage permits.
6. Pick arrivals using the paper's local three-sample quadratic maximum; fit travel time versus separation; report the slope, intercept, bootstrap uncertainty, and a permutation null that shuffles receiver positions within each day.
7. Compare the recovered slope and packet polarity with the released Lellouch correlograms before drawing any conclusion about the 2024--2025 installation.

## Current conclusion

The earlier five-month effort should be described as an **inconclusive reproduction attempt with a documented processing branch**, not as a failure of ambient interferometry and not as an unexplored analysis opportunity. The leading explanation is a mismatch between the published 2017 experiment and the 2024--2025 archive, compounded by a different virtual-source geometry and unnormalized correlation scaling. A controlled transfer test can distinguish “method does not transfer to this installation” from “previous implementation did not reproduce the published observable.”
