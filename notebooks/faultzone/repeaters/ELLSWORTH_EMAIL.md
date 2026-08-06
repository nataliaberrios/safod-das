# Bill Ellsworth's guidance on the repeater search — verbatim record

Recovered 2026-08-06 from the session transcript and filed here because it had
already been compacted out of context once and was being paraphrased from memory.
Ellsworth is Professor Emeritus of Geophysics at Stanford and a co-author of the
foundational Parkfield repeater and double-difference literature; this is the
expert guidance the repeaters workstream is supposed to be following.

---

## 2026-05-20 17:08 — William Ellsworth to Natalia Berrios-Rivera, "Search for repeaters"

> Hi Natalia,
>
> I took a few minutes this afternoon to see what potential repeaters might have
> been recorded by the DAS between May 2024 and April 2026. I didn't have the exact
> operation dates (or when data is actually available as I know that there are
> gaps), but it would be easy enough to do with the correct dates or even to add
> the brief recording period in 2017.
>
> I first went to https://ncedc.org/ncedc/catalog-search.html and selected the
> Double Difference Catalog in readable format using the start and stop times of
> 2024/05/01,00:00:00 and 2026.04.01,00:00:00 and the additional search parameter
> delta=35.982,-120.544,0,15. This goes in the box and all of the other defaults
> were set to blank. It found 329 events. They are in the attached text file
> DDRT_DAS.txt.
>
> I then read the file into r and ran the script in 'DAS repeater search.src.txt'.
> It uses a helper function date.to.decimal(). Both are attached.
>
> It made two plots (attached). One a cross section along the fault with
> earthquakes sized by magnitude and colored by depth (shallow is red, violet is
> deep) and the other a distance time plot with the same symbols. There are a
> number of promising potential repeaters.
>
> This is only a start, obviously, as you will want to find events that are close
> in magnitude and location, which will then need to be verified as either
> repeaters or neighbors.
>
> I hope that these files will prove to be useful. Let me know if you have any
> questions or want to discuss any of this.
>
> Regards,
>
> Bill

## 2026-05-20 17:37 — follow-up

> Forgot one more function… distaz.
>
> If don't use r but want to try it, you will first need to run these commands to
> make the functions available:
>
>     source('distaz.src')
>     source('date.to.decimal.src.txt')
>
> Note that you will either need to be in the same directory as the two functions
> or give their full pathname.

## 2026-05-21 — after Natalia's reply

> Hi Natalia,
>
> Here's a better cross section, this time colored by date. Spots with multiple
> colors are the ones I would focus on. **Larger events would be better too.**
>
> Regards,
>
> Bill

---

# What this actually settles

## 1. There is no ready-made Parkfield repeater catalog for 2024-2026

This is the important one. Ellsworth did **not** point to an established repeater
catalog covering the DAS period. He built candidates himself, from the NCEDC
**Double Difference (DDRT) catalog**, and said explicitly that they *"need to be
verified as either repeaters or neighbors."*

So the question "should we build our own catalog or use an established one" is
already answered: the expert built his own. Our approach is his approach plus the
verification step he named. Earlier worry in this project about circularity from
self-defined families was over-weighted — Nadeau & McEvilly and Waldhauser also
built theirs from waveform cross-correlation, so an "established" catalog is not an
independent definition of a repeater either.

What an external catalog would still add: decades of prior recurrence (for any
creep/recurrence section) and independent relative relocations to validate DAS-derived
offsets against. The second can be substituted with HRSN double-difference relative
location from data already cached in `hrsn_cache/`.

## 2. A selection criterion we have not used

> *"Here's a better cross section, this time colored by date. Spots with multiple
> colors are the ones I would focus on."*

Colour is **date**, so a multi-colour spot is one location occupied at several
different times. That is a **space-time clustering criterion applied to the DDRT
cross-section**, independent of waveform similarity. Everything in this project so
far has selected on similarity (HRSN beta, DAS CC). Bill's criterion is a genuinely
separate prior and gives a second, non-circular way to nominate candidates.

Worth implementing: it can nominate candidates that similarity ranking misses, and
agreement between the two selections is itself evidence.

## 3. "Larger events would be better too" — said twice, effectively

Independently confirmed by `hf_snr_test.py`: M3.17 at 1.8 km has SNR > 3 on 855
channels at 80-100 Hz, while M0.65 at 2.4 km is dead by 70 Hz. The confirmed pair
set used so far is M0.65-1.86, i.e. concentrated at the **weak** end. Both the
expert's advice and the instrument's own response say to weight toward larger,
nearer events.

## 4. Search geometry

`delta=35.982,-120.544,0,15` — centre 35.982, -120.544, radius 0-15 km, over
2024/05/01 to 2026/04/01, 329 events. Note this centre is the one recorded elsewhere
in this project as **1243 m from the surveyed wellhead** (the notebook/MH029 value
35.974204, -120.552141). Immaterial for a 15 km radius search, but the DDRT event
coordinates should be referenced to the surveyed collar, not to this centre.

## 5. Two standing offers

- *"it would be easy enough to do with the correct dates or even to add the brief
  recording period in 2017"* — he will redo the search given the **actual DAS
  operational dates and gaps**. That is the email to send, and it is a concrete
  request rather than a general question.
- *"Let me know if you have any questions or want to discuss any of this."*

Note his window ends **2026/04/01**; DAS coverage now extends past that, so a rerun
would add events. The 2017 offer is moot — `METHODS_STATUS` records that 0 of 26
sequence locations had an event in that 10-day window.

## Attachments — FOUND, in `notebooks/awd_clean/repeating/`

Not in `repeaters/`, which is why an earlier search missed them:

    DDRT_DAS.txt                  329 events, the DDRT extract
    DAS repeater search.src.txt   Bill's R script
    date.to.decimal.src.txt       helper
    distaz.src                    helper
    cross_section.png             the plot he refers to
    distance-time.png

### What the R script actually does

    junk1 <- delaz(35.982,-120.544,junk$lat,junk$lon)
    junk$n <- junk1$dist*cos(junk1$az*pi/180)
    junk$e <- junk1$dist*sin(junk1$az*pi/180)
    junk$r <- junk$n*cos(45*pi/180) + junk$e*sin(45*pi/180)
    junk$t <- -junk$n*sin(45*pi/180) + junk$e*cos(45*pi/180)

Three things worth pinning down, because this project has reimplemented the
geometry more than once:

1. **The fault strike is taken as 45 deg.** `t` is the SE (along-strike) component
   and `r` the NE (fault-normal) component. Bill labels `t` "Distance along strike".
2. **The reference point is called MH030** in his own comment — the monument this
   project's notes elsewhere flag as *not* the surveyed collar (the notebook /
   MH029 value 35.974204, -120.552141 sits 1243 m away). His centre is fine for a
   15 km radius search; it is not the right origin for source-relative geometry.
3. **The script only plots.** There is no candidate-selection algorithm in it.
   "There are a number of promising potential repeaters" and "spots with multiple
   colors" are conclusions Bill drew by eye from the two figures. So the
   space-time clustering criterion has to be written from scratch — it is his
   judgement, not his code.

---

# A prior screening effort exists, and its negative result is now explained

`notebooks/awd_clean/repeating/` also holds a full screen from 2026-08-02/03 —
`repeater_catalog_screen.py`, `repeater_candidate_pairs.csv` (+ ranked, + with
coverage), `repeater_event_waveform_coverage.csv`, `repeater_waveform_extract.py`,
`repeater_catalog_report.md`, `repeater_pilot_waveform_correlation.md`. It is not
referenced from `repeaters/` and was not used by any of the `g*` work.

It reports:

- 6325 candidate pairs (horizontal sep <= 1.0 km, depth sep <= 1.0 km, dM <= 0.50)
- 146/329 events with exact ambient-archive coverage; 1200 covered pairs
- a 10-pair waveform pilot: **median channel-wise max correlation 0.139-0.177, with
  0.0 % of channels exceeding 0.5**, concluding "this pilot does not confirm
  repeating earthquakes"

**That negative result is a processing artefact, and this project has since
diagnosed the cause.** The pilot computed *channel-wise* correlation on
individual channels with no moveout correction. Single DAS channels carry roughly
a 15 dB penalty against a clamped geophone (Lellouch et al. 2020), so per-channel
CC around 0.15 is exactly what a real repeater should give under that processing.
`METHODS_STATUS` sections 2.4 and 17 quantify it: a flat stack across the 0.31 s
moveout destroys ~27 dB of array gain, and correcting it takes the same class of
pair from CC 0.680 to **0.956**, at parity with HRSN's 0.954.

So the awd_clean pilot should be read as *"per-channel correlation without array
processing does not detect repeaters"*, which is true and consistent, rather than
as evidence against repeaters in this dataset.

### Two numbers that disagree between the trees

| quantity | `awd_clean/repeating/` | `repeaters/` |
|---|---|---|
| catalog events with DAS coverage | 146/329 | 208/329 (`phaseA_coverage.py`) |

Both are computed against the same 329-event DDRT list, so the definitions of
"covered" differ — most likely exact-segment versus file-interval overlap. Flagged
rather than silently resolved, per the repo convention. Resolve before either
number is quoted in a paper.

### What to reuse rather than rebuild

`repeater_candidate_pairs_ranked_with_coverage.csv` is a 6325-pair space-and-
magnitude screen with coverage already intersected. That is a far larger candidate
pool than the 10 HRSN-confirmed pairs this project has been recycling, and it is the
natural input to Bill's "spots with multiple colors" criterion. Read it; do not
modify anything under `awd_clean/`.
