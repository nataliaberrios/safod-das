# AWD plain-look figures

Five folders. **Start with `basic/`.**

| folder | what it is | made by |
|---|---|---|
| **`basic/`** | **the standard first-look set for the Nano (cemented) fiber — 12 figures. Use these.** | `awd_clean/basic.py` |
| **`basic_deep/`** | **the same set for the Deep (wireline) fiber — 11 figures** | `awd_clean/basic_deep.py` |
| `diagnostic/` | the dense multi-panel QC set, plus the taper audit and bad-channel check, plus the CSV/NPZ products | `awd_clean/plain_look.py`, `taper_audit.py`, `bad_channel_check.py` |
| `deep/` | the same looks for the Deep (wireline) fiber, hairpin split at channel 1702 | `awd_clean/deep_plain_look.py` |
| `virtual_source/` | correlation and deconvolution gathers — source removed, arrival recovered at 2948 m/s | `awd_clean/awd_virtual_source.py` |
| `superseded/` | two earlier attempts at a presentation set. Kept for history; do not use. | `awd_clean/look.py`, `plain_look_simple.py` |

## `basic/` — what each one is

| file | shows |
|---|---|
| `basic01_gather_one_drop` | shot gather, a single weight drop, 20–50 Hz |
| `basic02_gather_stack` | the same but 859 drops stacked |
| `basic03_traces` | stacked traces every 50 m |
| `basic04_spectrum` | amplitude spectrum, signal window against pre-drop noise |
| `basic05_fk` | f–k spectrum, with 2975 and 1500 m/s lines drawn |
| `basic06_amplitude_vs_distance` | RMS signal and noise against distance |
| `basic07_snr_vs_distance` | SNR in dB against distance |
| `basic08_moveout` | relative traveltime, fitted apparent velocity **2948 m/s** |
| `basic09_bands` | one drop through six bands, unfiltered to 100–250 Hz |
| `basic10_spectrogram` | a whole 5-minute raw file at 149 m, drop times marked |
| `basic11_spectrum_vs_distance` | amplitude spectrum against distance |
| `basic12_survey_fold` | drops per burst across the 24-hour survey |

## Numbers worth knowing before quoting any of these

- **2948 m/s** apparent velocity from `basic08` (Nano), over 56–598 m. Independent
  of, and within 1% of, the 2975 m/s established elsewhere in the project.
- **1576 m/s** from `basic_deep/dbasic08`, over 1013–3275 m along the outbound leg.
  Sits just above the 1.4–1.56 km/s quoted for the Deep slow-mode candidate, from
  a plain cross-correlation with no wedge or scan. Only 43 points cleared the
  6 dB SNR cut and all lie beyond 1000 m, so the shallow leg is not sampled.
- **CC 0.95 → 0.18** between 150 m and 550 m for drops within one burst: the
  source repeats very well shallow and not at all deep.
- **SNR 28.3 dB at 460 m** in `basic07`. This is the **859-drop** stack. The
  figure guide's 7.6 dB at the same distance is for **burst** stacks of ~19
  drops; √(859/19) ≈ 16 dB accounts for the difference. Do not put the two
  numbers side by side without saying which stack each refers to.
- **Distance along fiber, not depth**, on every axis. The Nano fiber is the same
  cemented main-hole fiber Lellouch et al. used, so the depth convention is
  inherited — see the figure guide — but the axes here are the raw fiber
  coordinate.

## Caveats

- `basic/` is the **Nano** (cemented) fiber; `basic_deep/` is the **Deep** (wireline)
  fiber. They are not interchangeable: Deep is differentiated to strain rate to
  match Nano, uses the 3-15 Hz working band rather than 20-50 Hz, and shows the
  outbound leg only because the fiber hairpins at channel 1702.
- All single-drop figures use burst 18, chosen mechanically as the fullest burst
  contained in one raw file, not by how it looks.
- Record sections have ~25 instrumental bad channels interpolated over. They are
  narrower than the 16.5 m gauge length so they cannot be ground motion; left in
  they draw stripes across the whole panel. `diagnostic/fig11_bad_channels.png`
  is the evidence.
- `deep/deep_fig06_taper_ramp_test.png` has a known artifact: the synthetic gain
  is applied to the 3.5 s excerpt rather than the whole file, which creates a
  filter edge transient the real reader never produces. The tilt is right, the
  absolute levels are not.

`FIGURES_EXPLAINED.md` in this folder describes the `diagnostic/` set in detail
and carries a glossary of the acquisition vocabulary.
