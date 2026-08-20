### 4.7 More data does not help (Figure 8)

{{FIGURE:fig8_convergence}}

Section 4.2 reported that our coherent four-day stack (p = 0.9184) sits further
from significance than our best single day (p = 0.1345). That comparison is
suggestive but it varies the days as well as the stack length, so it does not
isolate the effect of adding data. Here we isolate it. On one contiguous day,
2024-12-20, we stack 1, 2, 4, 8, 16 and 24 hourly chunks of the same day's
correlations and fit `score ~ N^b`, in both the baseline and the
common-mode-removed branch.

The quantity to read is **not** the raw score. A raw score can rise merely because
a repeatable contaminant accumulates coherently, which is exactly what section 4.4
says this field contains. We therefore report **detectability**: the score divided
by the 95th percentile of its *own* receiver-order null, rebuilt from the same
stacked data at every stack length. A contaminant raises the observation and its
null together and leaves detectability flat, whereas averaging suppresses
incoherent noise as `1/sqrt(N)`, so detectability grows as `N^+0.50` if and only if
there is a coherent arrival to accumulate.

| chunks | windows | baseline detectability | baseline p | c-m removed detectability | c-m removed p | c-m pedestal |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 240 | 0.858 | 0.9005 | 0.834 | 0.8691 | +0.915 |
| 2 | 480 | 0.950 | 0.4113 | 0.754 | 0.9580 | +0.813 |
| 4 | 960 | 0.858 | 0.9800 | 0.765 | 0.9125 | +0.869 |
| 8 | 1,920 | 0.966 | 0.2234 | 0.744 | 0.8946 | +0.892 |
| 16 | 3,840 | **1.014** | **0.0170** | 0.769 | 0.6807 | +0.946 |
| 24 | 5,759 | 0.971 | 0.1614 | 0.921 | 0.8836 | **-0.219** |

| fitted exponent | baseline | c-m removed | stacking works |
|---|---:|---:|---:|
| raw score | N^+0.158 | N^-0.004 | N^+0.50 |
| **detectability** | **N^+0.042** | **N^+0.019** | N^+0.50 |

Both branches are flat. Twenty-four times the data buys a factor of 1.14 in
detectability where a coherent arrival would buy a factor of 4.9.

**The one threshold crossing, reported as what it is.** At 16 chunks the baseline
reaches detectability 1.014 at p = 0.0170, which in isolation reads as a detection.
At 24 chunks, with 50 % more data, it falls back to 0.971 at p = 0.1614. We report
the trajectory rather than either endpoint because the shape is the diagnostic: a
genuine arrival does not un-detect itself when data is added. It is the same
failure mode as the selected day-pairs of section 4.2 and the unfiltered adaptive
f-k configuration of section 4.6 -- a threshold crossed by a statistic whose
fluctuation exceeds the effect being sought.

**Why the flat curve indicates absence rather than contamination.** A flat
detectability curve admits two explanations: there is no arrival to accumulate, or
there is one and the limiting contaminant accumulates exactly as it would. The
common-mode-removed branch discriminates. At the full 24-chunk stack its pedestal
diagnostic is **-0.219**, so by the criterion of section 3.2 the repeatable
component is suppressed and the statistic is finally measuring moveout rather than
proximity to the zero-lag lobe; consistently, its peak leaves the scan ceiling for
3,575 m/s, inside the physical fan. With the contaminant demonstrably suppressed,
detectability is still 0.921 at p = 0.8836 and the exponent is still N^+0.019.
Contamination does not account for the flat curve. Absence does.

This is the quantitative form of the comparison with Behm (2016), who obtained
robust borehole interferograms from 30 s of noise under adequate illumination. If
30 s suffices when the wavefield contains the arrival, then a day that yields no
growth in detectability is not an under-sampled recording.

**Limits.** The scan covers one day at hourly granularity, so it measures
within-day convergence; the four-day stack of section 4.2 covers the across-day
case and also degrades. The 24-chunk row uses 5,759 windows, the full contiguous
day of section 3.1, so the endpoint is the same quantity reported there. The null
is a receiver-order permutation of the finished gather and therefore contains no
pre-correlation operator -- the limitation identified in section 4.6 -- so this
test bounds stack length, not operator choice.
