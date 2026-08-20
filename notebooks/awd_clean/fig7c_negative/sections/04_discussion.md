## 5. Discussion

### 5.1 Why borehole body-wave illumination is fragile

Interferometry reconstructs a Green's function only from sources in the
stationary-phase zone of the path in question. For two channels on a **vertical**
line the direct P path between them is vertical, so the stationary-phase sources
lie on the extension of that line: a narrow cone directly above the wellhead,
radiating steeply downward. Noise arriving at oblique incidence contributes nothing
to that path and merely adds to the background against which the arrival must be
detected.

This is a much narrower requirement than the one governing surface-wave
interferometry, where sources anywhere along the receiver azimuth contribute, and
it explains an asymmetry in the DAS literature that is otherwise easy to misread.
Ambient-noise DAS studies overwhelmingly report *surface-wave dispersion*
inversions, and they succeed routinely. Body-wave retrieval on a vertical array is
a different problem with a far more restrictive source requirement, and the
frequency with which the surface-wave case succeeds should not be read as evidence
that the body-wave case is robust.

What fills that cone is local and typically cultural: machinery and vehicles on the
wellhead pad, drilling, pumps, work on an access road, wind coupling into the
casing. Ocean microseism and regional earthquakes do not, because they do not
arrive at near-vertical incidence. Illumination is therefore a property of *site
activity*, and site activity is not a property of the rock.

The consequence for reproduction is direct and, we think, under-appreciated. Same
fibre, same borehole, same rock guarantees the same Green's function; it does not
supply the sources needed to reconstruct it. A published ambient-noise body-wave
result is a statement about a wavefield during a recording window, not a
transferable property of an installation. Behm (2016) is the informative contrast:
his success came in a producing oil field under water-flooding, about the noisiest
wellhead environment obtainable, and even there the illumination was diagnosed from
the interferograms rather than assumed.

### 5.2 The stacking trap

The natural response to a weak or absent arrival is to stack more. Under adequate
illumination that is unnecessary -- Behm obtained robust interferograms from 30 s --
and under inadequate illumination it is counterproductive in a specific way that our
results display.

Stacking suppresses incoherent noise as `1/sqrt(N)`, but it does **not** suppress a
coherent artefact. The static pattern of section 4.4 is perfectly coherent, so it
accumulates at full strength while any genuine arrival, if present at all, grows no
faster. We measured this directly: the pedestal diagnostic is `-0.381` for a single
day but `+0.951` once four days are coherently stacked. The artefact wins as the
stack lengthens, which is why our 96-hour stack (p = 0.9184) is worse than our best
single day (p = 0.1345).

Section 4.7 puts a number on it at fixed geometry. Stacking 1 to 24 hourly chunks
of one day grows detectability as `N^+0.042` without common-mode removal and
`N^+0.019` with it, against `N^+0.50` for a coherent arrival accumulating against
incoherent noise. The second of those exponents is the load-bearing one: it is
measured with the pedestal diagnostic at `-0.219`, so the coherent artefact has been
suppressed and the curve is still flat. That distinguishes the two readings of a
non-improving stack -- an artefact masking a signal, or no signal -- in favour of
the second.

The practical rule this suggests: **a stack that degrades with added data is
evidence about the wavefield, not a reason to add more.** Convergence behaviour
should be reported as a matter of course in ambient-noise studies, because it
distinguishes a weak signal from an absent one at no additional cost.

### 5.3 Selection regions and the manufacture of apparent arrivals

Our first F-K result looked positive: a fixed 2.5-4.5 km/s fan produced a
convincing ridge. It failed when we required the real-data result to survive an
**input-level** null, in which the surrogate is built before the operator runs so
that the operator is inside the null. Under a pre-filter channel scramble the ridge
persisted, which means the fan was largely manufacturing it.

This is a general hazard for narrow selection regions applied to low-wavenumber-
dominated data, and it is worth stating as a methodological point rather than as a
local mishap. At a 700 m aperture the 2.5-4.5 km/s fan is barely more than one
wavenumber resolution cell wide at the low end of a 5-20 Hz band and sits within
about 1.4 cells of k = 0 (section 4.5). A filter that narrow, applied that close to
a dominant zero-wavenumber peak, is capable of shaping the pedestal into something
that resembles moveout. Post-filter nulls will not catch it, because they inherit
the operator.

We therefore recommend that velocity-fan results in this regime be reported with an
input-level null, and that the fan's width be compared explicitly against the
array's wavenumber resolution.

### 5.4 What would be required

For this site, the ambient body-wave route is closed at 5-20 Hz on this archive,
and no processing development would change that: the wavefield does not contain a
net downgoing component to recover. Three routes remain open.

**Earthquakes rather than ambient noise.** Lellouch's Figure 9 velocity model was
mainly earthquake-derived, and the scientific target -- a velocity profile -- is
therefore still reachable. Our archive has 206 catalogued events cached, and an
independent channel-depth registration already matches the 2005 check shot to 0.2 %.
This is the route we recommend for the velocity objective.

**A controlled source, or a documented noisy interval.** If a downgoing wavefield
is required, it can be supplied. The asymmetry measurement of section 3.3 provides
an inexpensive way to identify intervals in which site activity happens to provide
it, should any future archive contain them.

**A feasibility gate before acquisition.** The measurement needs seconds of data,
no velocity assumption, no correlation stack and no filter. Run on our archive it
would have returned a null immediately, at a cost of minutes rather than the
substantial compute and analysis this study consumed. We would run it first next
time, and we suggest others do the same.

### 5.5 Value of the negative

Three things here generalise beyond one figure.

1. **A reproduction failure with a working positive control is informative.**
   Because the same code returns the published velocity range from the original
   author's data (depth-velocity r = 0.948) and physically impossible velocities
   from ours, the null is attributable to the recording rather than to the implementation. Reproduction studies without such a control cannot make
   that separation, and we would encourage the pattern.
2. **Illumination should be measured, not assumed.** The asymmetry statistic is
   cheap, needs no velocity assumption, and is diagnostic. Reporting it alongside
   ambient body-wave results would let readers judge whether a positive result
   reflects favourable illumination or a robust method.
3. **A single site can yield both a positive and a null.** The 2017 positive and
   the 2024-2025 null come from one fibre in one borehole. This bounds how far any
   ambient body-wave demonstration can be read as evidence of general feasibility,
   which matters directly for proposals that assume such measurements can be made
   routine for monitoring.
