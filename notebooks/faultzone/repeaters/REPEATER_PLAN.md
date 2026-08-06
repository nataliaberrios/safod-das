# Identifying repeating earthquakes on SAFOD DAS — the full procedure

Written 2026-08-05. Supersedes the ad-hoc sequence of scripts that produced the
current 8-pair set. The aim is a defensible catalog, not a fast one.

---

## The question, in Bill Ellsworth's words

> *"This is only a start, obviously, as you will want to find events that are close
> in magnitude and location, **which will then need to be verified as either
> repeaters or neighbors**."*

That distinction is the whole problem. Two events can be

- 1–8 m apart in the double-difference catalog,
- within 0.1 magnitude units,
- and correlate at 0.99 on a borehole seismometer,

and still be **neighbours** — adjacent patches on the same fault — rather than the
same patch rupturing twice. Only the latter measures creep, because only the latter
gives a recurrence interval that means anything.

Nothing done so far in this project distinguishes the two. The 8 pairs currently
called "confirmed" are confirmed as *similar*, which is a weaker statement.

His other guidance, from the same thread:

- DDRT catalog, NCEDC, 2024/05/01–2026/04/01, `delta=35.982,-120.544,0,15` → **329
  events** (`DDRT_DAS.txt`)
- *"Spots with multiple colors are the ones I would focus on"* — multiple dates at
  one location, on the date-coloured cross-section
- *"Larger events would be better too"*
- *"or even to add the brief recording period in 2017"*

---

## What already exists, and what is wrong with it

| asset | state |
|---|---|
| Ellsworth's DDRT catalog + R script | authoritative starting point, `awd_clean/repeating/` |
| Prior screen: 6325 candidate pairs, 1200 with coverage | good; coverage audit used only one of two manifests |
| My 8 HRSN-confirmed pairs | **all 8 appear in the 6325**, separations 1–8 m — mutual validation |
| `correlate_perchannel.py` | validated detector: per-channel CC → SNR-weighted correlogram stack, fast path asserted against explicit loop |
| G0 registration | wellhead ≈ channel 23; channels 23–896 uniformly coupled |
| Moveout correction | +27 dB; `moveout_test.slant_scan` fits apparent slowness per event |
| `hrsn_extend.py` + 96 cached HRSN events | independent-instrument confirmation |
| `cache_all/` | 206 DAS event windows already extracted |

**The central deficiency: 128 well-screened covered candidate pairs exist and I
have tested 8.** The prior pilot tested 7 more, and *none of them* were pairs that
later confirmed — their ranking (compactness × magnitude) does not order waveform
similarity, exactly as Gao et al. 2021 warns.

---

## Phase A — the candidate set, built correctly

**A1. Coverage audit against both manifests, at file-interval resolution.**
`SAFOD_2024_2025.csv` spans 2024-05 → 2025-07-25; `SAFOD_vertical_2026_01_28.csv`
spans 2025-05-18 → 2025-09-24; they overlap. My own event list used *day-level*
coverage, which is too coarse — a day can be listed while the specific minute is a
gap. Three of my four "uncovered" pairs fall inside manifest 1's span, so those are
real recording gaps, not bookkeeping.
*Output:* every DDRT event with an exact containing interval, and every pair with
both events covered.

**A2. Two candidate streams, kept separate.**
- *Location-first* (Bill's): close in location and magnitude, prefer larger events.
- *Similarity-first*: all covered pairs correlated regardless of catalog position.

They must be kept apart because each is the other's control. A pair found both ways
is far stronger than one found either way. The prior work did only the first; my
`correlate_all.py` did only the second.

## Phase B — detection with the validated chain

Per-channel CC → SNR-weighted correlogram stack (Lellouch et al. 2021), channels
23–896, moveout-corrected using each event's own fitted slowness, bands 2–8, 4–16
and 10–40 Hz, ranked by the **weakest** band so a pair must survive all three.

Threshold from the empirical null: 6×MAD, with the **time-reversed acausal floor**
measured over the full pair set rather than n=2. No imported absolute threshold.

## Phase C — repeaters or neighbours: the verification Bill names

Five tests. A pair is a repeater only if it passes all that apply; each is reported
separately so a reader can disagree with any one of them.

**C1. Source-area overlap.** Separation < source radius at a stated stress drop.
r ≈ 18 m at M1.0, Δσ = 3 MPa; use Abercrombie 2014's *measured* Parkfield drops
rather than assuming 3 MPa. Current pairs are 1–8 m apart, so they pass — but the
DD relative precision is comparable to r, so this test is weak on its own and must
be stated as such.

**C2. DAS moveout consistency — the test only this instrument can do.**
The fibre measures apparent slowness along the borehole directly. Two events from
the *same patch* must arrive with the *same* incidence geometry, hence the same
apparent slowness, to within measurement error. Neighbouring patches a source
radius apart differ slightly but systematically. `slant_scan` already fits this per
event; it needs a finer slowness grid (the current grid quantises to ~250 m/s) and
a proper error estimate. **This is a geometric co-location test independent of
waveform similarity, and no seismometer network can make it as directly.**

**C3. Differential S–P across HRSN.** Same patch ⇒ identical S–P at every station.
Scale-free, needs no absolute timing, and is the standard check
(Uchida & Bürgmann 2019 list it among the "multiple data constraints" required).

**C4. Coda coherency decay.** Neighbours decorrelate with lapse time faster than
repeaters, because their scattered paths diverge. Measure CC as a function of lapse
window — already implemented in `coda_window_survey.py`.

**C5. Magnitude and moment consistency.** ΔM ≤ 0.3 (Waldhauser & Ellsworth 2002).
Three of the current 8 fail this; C1–C4 decide whether that is disqualifying or
whether the magnitudes are simply uncertain.

## Phase D — independent confirmation on HRSN

Everything above is DAS. HRSN is the independent instrument and the array the
Parkfield repeater catalogs were built on. Convention: mean CC ≥ 0.9 at ≥5 common
stations (Waldhauser & Schaff 2008). Report DAS and HRSN verdicts separately; do
not merge them into one number.

## Phase E — the 2017 extension (Bill's suggestion)

No 2017 data on `$OAK` — the only pilot-fibre holding is March–April 2026. Lellouch's
public archive (`github.com/ariellellouch/SAFODDAS`, 800 ch × 1 m, 250 Hz) is event
snippets, so it supports testing *specific* events, not a search. If any Phase C/D
sequence also appears in 2017, the baseline goes from 14 months to **nine years**,
which is the difference between a marginal recurrence estimate and a real one.

Also check the sequences against Waldhauser–Schaff relocations (1984–2019); a
sequence active then and now inherits decades of history.

## Phase F — recurrence → creep

Only after A–E. Two slip estimates side by side, because they disagree by design:
Nadeau & Johnson 1998's empirical scaling (now in `planpapers/01_repeaters/`) and a
circular-crack model. Note that N&J was *calibrated against geodetic creep at
Parkfield*, so using it to measure creep here is partly circular — say so.

---

## Standing rules for this work

1. **No threshold is imported.** Every cut is calibrated against a null measured on
   the same instrument, band and processing. This project has twice reported a
   false negative from an imported threshold.
2. **Never correlate at zero lag.** Catalog origin times carry 0.1–0.5 s; at
   5–20 Hz that drives identical waveforms to CC ≈ 0. Use `dvv_core.bulk_align`.
   This error has occurred twice.
3. **Counts are quoted with their criteria attached.** "7 pairs" and "83 pairs"
   were the same analysis at different cuts.
4. **Controls run before the result is believed, not after.** Five times in this
   project a control caught something that already looked like a result.
5. **Read-only outside `nberrios/` and the shared SAFOD data.**

## Open items needing a person

- The earlier student's repeater catalog — not found by search; needs a name.
- Whether the SAFOD Phase-2 georeferenced-channel file resolves the ±25 m depth
  datum (ask, don't read).
- Ellsworth: are these sequences in his Parkfield catalogs, and what recurrence
  history do they have?
