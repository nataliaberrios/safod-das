# What each plain-look figure shows, in plain language

> **Note:** the figures were reorganised into subfolders. The `figNN_*` files
> described here now live in `diagnostic/`, `deep_fig*` in `deep/`, `vs_fig*` in
> `virtual_source/`. For the current first-look set see `basic/` and the
> `README.md` beside it.

Companion to the PNGs in this folder. Nothing here changes the figures; it just
explains them. If you only read one section, read "How to read a record section"
and then figures 2, 6 and 9.

---

## How to read a record section

Most of these figures are the same kind of picture, so it is worth getting it
once.

- **Vertical axis = distance along the fiber, in metres.** Zero is the top. It
  increases downward, so the picture is oriented like the borehole.
- **Horizontal axis = time, in seconds, with 0 = the weight drop.** Negative time
  is before the drop.
- **Colour = ground motion at that place and time.** Red one way, blue the other,
  white is nothing. The colours have no meaning beyond sign.

Each horizontal row of pixels is one channel — one measurement point on the
fiber. A seismic wave shows up as a **diagonal streak**, because it reaches
shallow channels before deep ones. The slope of that streak is the speed: steep
means fast, shallow-angled means slow. That slope is called *moveout*, and it is
the main thing to look for.

Two more things that recur:

- **"trace-normalised"** means every row was rescaled to its own maximum, so weak
  and strong channels are equally visible. Good for seeing shape, useless for
  comparing strength between channels.
- **Nano** is the cemented shallow cable (~930 m). **Deep** is the wireline cable
  that goes much further and folds back on itself.

---

## Figure 1 — one drop, nothing done to it

`fig01_raw_record_section.png`

The rawest possible view: one weight drop, no filtering, no cleanup. Top row
Nano, bottom row Deep. Left is true amplitude, right is trace-normalised.

**What to notice.** In the Nano panels you can just make out a diagonal streak in
the top ~100 m in the first half-second. Everything below that looks like static.
That is honest — raw DAS really is mostly noise, and the signal only emerges once
you filter (figure 2).

**The Deep panels look broken, and that is the finding.** They show flat
horizontal stripes with no time structure at all. That is because each Deep
channel sits on a large constant offset that swamps everything else, so
normalising by the peak just shows the sign of that offset. It is the reason the
reader removes a linear trend before doing anything. Do not try to read physics
off the bottom row — read it as "Deep must be detrended first."

---

## Figure 2 — the same drop in six frequency bands ⭐

`fig02_frequency_bands.png`

The single most useful picture here. Same Nano drop six times, each filtered to a
different frequency range, from unfiltered through to 100–250 Hz.

**What to notice, panel by panel:**

- *unfiltered* — static. No visible signal.
- *1–5 Hz* — a broad smear covering the whole array at once, starting **before**
  the drop. Arriving everywhere simultaneously means it is not a wave travelling
  down the hole; it is background noise shared by every channel.
- *5–20 Hz* — a clear diagonal arrival appears, reaching ~300–400 m.
- *20–50 Hz* — the cleanest version. Sharp, steep, traceable to ~400 m.
- *50–100 Hz* — visible only in the top ~200 m.
- *100–250 Hz* — nothing.

**Why it matters.** This is the picture that justifies working at 20–50 Hz. You
are not taking that band on faith or inheritance — you can see the signal live
there and nowhere else.

The vertical stripes at the far left and right edges of some panels are filter
artifacts, not data. Ignore the outer ~0.1 s.

---

## Figure 3 — where the noise and the signal sit in frequency

`fig03_noise_psd.png`

Left column: spectra — how much energy there is at each frequency — for a few
depths, plus the pre-drop noise as a dashed line. Right column: the same noise
information as an image, one row per channel, so you can see how the noise
changes with depth. Top row Nano, bottom row Deep.

**What to notice.** The dotted vertical lines mark 60 Hz and its harmonics, where
electrical pickup would appear. It was measured and **there is no 60 Hz line in
this data** (0.2 dB, i.e. nothing) — worth knowing, because it is a standard thing
to be asked about.

On the right-hand image, a horizontal streak that is much brighter than its
neighbours is a noisy channel. Figure 11 chases those down.

---

## Figure 4 — five minutes of data, all frequencies

`fig04_spectrogram.png`

A spectrogram at three depths. Horizontal axis is time through the whole raw
file, vertical axis is frequency, colour is energy. Cyan lines mark where the
manifest says drops happened.

**What to notice.** Each drop should appear as a brief vertical stripe at a cyan
line. This is the sanity check that the timestamps line up with the data, and
it is also where you would spot anything odd going on in the background — a
machine turning on, a persistent tone, a gap.

---

## Figure 5 — which channels are misbehaving

`fig05_channel_qc.png`

Four per-channel health checks along the fiber, from one whole 5-minute file.

- *RMS* — how loud each channel is. Spikes are channels much noisier than their
  neighbours.
- *DC offset* — the constant level each channel sits on.
- *lag-1 autocorrelation* — a whiteness measure. Near zero means the channel is
  seeing pure noise; away from zero means it is seeing something with structure.
- *flat-topping* — a clipping check.

**What to notice.** No channels are dead, and there is no clipping. Beyond about
100 m the whiteness measure flattens out at about −0.32 and stays there, which
means that over a quiet 5 minutes the deeper fiber is sitting at its own noise
floor. A few channels spike well above their neighbours — but see figure 11,
because those turned out not to be permanent.

---

## Figure 6 — the same drop, twenty times ⭐

`fig06_drop_repeatability.png`

Twenty drops from one burst, drawn on top of each other, at three depths
(150 m, 350 m, 550 m). Left column unfiltered, right column 20–50 Hz. Each
panel's title gives the drop-to-drop scatter for that depth.

**What to notice.** In the 20–50 Hz column at 150 m the twenty traces lie almost
exactly on top of each other — that is excellent repeatability, visible without
computing anything. At 350 m they still follow the same shape but spread apart.
At 550 m they no longer agree at all.

**Why it matters.** That is your usable depth range, shown rather than argued.
The unfiltered column shows why filtering is not optional: the same twenty traces
are indistinguishable static.

The scatter at the very left and right edges of each panel is filter ringing, not
real variability.

---

## Figure 7 — what each cleanup step actually does ⭐

`fig07_processing_ladder.png` (and `fig07b_ladder_check.png`)

Before you see any data, the file reader silently does four things: removes a
linear trend, applies a taper, bandpass-filters, and subtracts the across-channel
median. This figure applies them one at a time so you can see each one's effect.
The title reports the overall signal level after each step.

**What to notice — this is the punchline.** The numbers go
`0.418 → 0.418 → 0.418 → 0.104 → 0.104`.

Only step 4, the bandpass, changes anything. Detrending does nothing measurable
here, the taper weight at this drop's position in the file is exactly 1.000 so it
does nothing at all, and the median removal does nothing you can see. Three of
the four steps are no-ops for this section.

The bottom-right panel is the common-mode trace that step 5 removes — it is
small.

`fig07b` is a cross-check that the hand-built version matches what the reader
actually returns. It is a debugging figure; you can skip it.

---

## Figure 8 — the two fibers are not measuring the same quantity

`fig08_units.png`

Nano records *strain rate*; Deep records *strain*. They differ by a time
derivative. This shows a Deep record both ways, as sections, as one trace, and as
spectra.

**Why it matters.** If you ever compare amplitudes or waveform shapes between the
two fibers, they must first be put in the same units, or the comparison is
measuring the instrument rather than the earth. Three scripts in the repo already
do this conversion; the paired-stack and record-section paths do not.

---

## Figure 9 — with and without common-mode removal ⭐

`fig09_common_mode_before_after.png`

"Common-mode removal" means subtracting, at every instant, the median across all
channels — the idea being that anything every channel sees at once is instrument
noise rather than ground motion. It is on by default. This shows the data both
ways: before, after, and the difference, plus traces, spectra, and how much each
channel loses.

**What to notice.** At 20–50 Hz it barely does anything — the removed part is
**0.9%** of the signal, and the before/after traces sit on top of each other.

**Why it matters.** You do not have to worry about this step distorting the AWD
results, because it is nearly a no-op in the band you work in. It matters much
more at low frequency, which is what `nano_fiber_end_coherence.py` already found.
The caveat worth remembering: common-mode removal would also delete a real
arrival that hit every channel simultaneously, so it is an assumption, not a
cleanup.

---

## Figure 10 — the taper problem

`fig10_taper_audit.png`

The reader applies a smooth fade-in/fade-out (a "taper") across each whole file.
Drops are timed by the survey, not by file boundaries, so whether a drop lands in
the faded part is pure luck.

**What to notice.** Panel A shows the fade curve with the drops plotted on it —
lots land on the slopes. Panel C shows the resulting amplitude bias drifting
across the survey with a bump in the middle, which is not physical. Panels E/F
show that for about a third of drops the gain changes *within* the 3.5 s window,
which distorts waveform shape rather than just scaling it.

**Why it matters, and why not to panic.** It affects amplitude only. Timing,
correlation and velocity results cannot be touched by it, and those are your
headline results. See `awd_clean/GUIDE_EDITS_TODO.md` for exactly what is and is
not affected.

---

## Figure 11 — are the noisy channels broken, or real?

`fig11_bad_channels.png`

A noisy channel can be an optical problem or a real mechanical one. The
distinguishing test is **gauge length**: the instrument averages over 16.5 m
(13 channels), so anything real must be at least that wide. Anything narrower
cannot be ground motion.

**What was found.** Nine channels flagged, all between 10 and 27 m — the very top
of the fiber. Six are narrower than the gauge length, so instrumental. And they
are intermittent, appearing in only two to four of the six files sampled.

**Important correction:** the dramatic spikes you can see at 137 m and 510 m in
figure 5 do **not** persist across the survey. They were in one file. There is no
permanent bad channel near 510 m.

Ignore the "neighbour correlation" panel — that test does not work as designed
and returns nothing meaningful.

---

## Deep-fiber figures

`deep_fig01` … `deep_fig06`. Same ideas as above, for the wireline cable.

One thing to know first: **the Deep fiber folds back on itself at channel 1702.**
It goes down the hole and comes back up, so channel number is not depth. The
figures split it into an outbound leg and a return leg and show them separately.
`deep_fig01` also folds them together — if the depth mapping is right, the two
legs should agree.

`deep_fig06` tests whether the taper problem changes what the Deep record looks
like. Caveat: the synthetic gain in that figure was applied to the short excerpt
rather than to the whole file, which creates an edge artifact the real reader
never produces. The tilt in the lower-left curve is real; the absolute levels
are not. That figure is due a fix.

---

## Virtual-source figures

`vs_fig01`, `vs_fig02`, `vs_fig03`.

Different idea from everything above. Instead of looking at the wavefield from
the weight drop, you correlate one channel against all the others. Because the
source is the same in both, it cancels — and what is left behaves as if a source
sat at that channel, down in the borehole. This is what the advisor was asking
for.

`vs_fig03` is the one to look at: correlation and deconvolution side by side,
plus wiggle traces. The clean diagonal streak leaving the green line is the
redatumed arrival, travelling at roughly 3 km/s. Deconvolution is visibly sharper
than correlation.

**Why it matters.** It removes the weight drop's own variability, which
`docs/paper1/STATUS.md` already identifies as the thing limiting your
sensitivity.

---

## Glossary: acquisition words, translated

Written for someone comfortable with the maths but not with seismic-acquisition
vocabulary. Nothing here is deep; the words are just borrowed from a field with
its own dialect.

**Channel.** One measurement point along the fiber. The DAS interrogator reports
strain over a short segment every ~1.27 m (Nano) or ~2.04 m (Deep). Think of it
as a spatial sampling grid: the fiber is a 1-D array of ~732 sensors.

**Trace.** The time series recorded by one channel. A row of the data matrix.

**Gather.** A set of traces plotted together, usually as an image with time on
one axis and position on the other. "Record section" and "shot gather" are the
same thing for our purposes. It is just the 2-D field u(z, t) rendered as a
picture.

**Moveout.** How arrival time varies with distance. For a wave of speed v
travelling along the array,

        t(z) = t0 + z / v

a straight line in the (time, distance) plane — which is why arrivals look like
diagonal streaks. "Apparent"
speed because you only measure the component along the fiber; a wave crossing
the fiber obliquely looks faster than it is.

**Slant stack.** Integrate the wavefield along those straight lines for a whole
range of trial speeds v. This is a Radon transform. The output tells you which
speeds carry coherent energy, without picking any arrival by hand.

**Semblance.** Normalised coherence along a trial path. Take the N trace
amplitudes a1 ... aN sampled along the path, then

        S  =  (a1 + a2 + ... + aN)^2  /  ( N * (a1^2 + a2^2 + ... + aN^2) )

Cauchy-Schwarz puts S between 0 and 1. S = 1 means every trace is identical
(perfectly in phase). Independent noise gives S about 1/N. So with 300 traces,
S = 0.45 is enormous and S = 0.05 is marginal.

**F-K.** The 2-D Fourier transform of the gather: frequency f against wavenumber
k. A plane wave of apparent speed v maps onto the straight line f = v * k, so
selecting a wedge in the (f, k) plane selects a range of speeds. The sign of k
gives the direction of travel; "signed F-K" just means keeping the two
directions separate.

**Strain vs strain rate.** DAS measures how much a fiber segment stretches.
Some interrogators report strain, others report its time derivative, the strain
rate. Nano reports rate; Deep reports strain. In the frequency domain the two
differ by a factor of i*omega -- an amplitude slope rising across the band, plus
a 90-degree phase rotation. Harmless within one fiber, wrong if you compare the
two without converting.

**Taper.** A smooth window multiplying the record, fading in at the start and out
at the end, so the FFT does not see a step discontinuity at the edges. A Tukey
window with alpha = 0.4 is flat over the middle 60% and cosine-tapered over the
outer 20% at each end. Standard and harmless — unless, as here, it is applied
over a different span than you assumed.

**Common mode.** The part of the signal identical on every channel at a given
instant, estimated as the across-channel median. Usually laser and interrogator
noise. Subtracting it is standard, but note it also removes any real arrival with
no moveout.

**Detrend.** Subtract a least-squares straight line from each trace. Removes DC
offset and slow drift, which on the Deep fiber are large enough to swamp
everything else.

**Virtual source / redatuming.** Cross-correlate the trace at A with the trace at
B. The source term is common to both and cancels; what survives behaves as if a
source sat at A and was recorded at B. You have moved the source from the surface
to a point inside the borehole without moving any equipment.

**Deconvolution interferometry.** Same idea, but divide B by A in the frequency
domain instead of multiplying by A's conjugate. That cancels the source spectrum
as well as the source timing, at the cost of blowing up noise wherever A is
small -- hence the "water level", a small constant added to the denominator to
stop it dividing by nearly zero.

**SNR in dB.** Ten times the base-10 logarithm of a power ratio. 0 dB means
signal equals noise, 10 dB means ten times the power, 20 dB a hundred times.

**Burst.** A group of weight drops fired close together in time. 989 drops in
this survey, grouped into 49 bursts of roughly 20.

**Stack.** Average repeated records. Coherent signal adds in proportion to N
while incoherent noise adds only as the square root of N, so signal-to-noise
improves as sqrt(N). This is why 859 stacked drops show structure that a single
drop does not.
