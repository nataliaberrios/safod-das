"""Build the audited advisor-facing Figure 7c reproduction and F-K QC notebook."""
from pathlib import Path
import nbformat as nbf

HERE = Path(__file__).resolve().parent
OUT = HERE / "Ambient_FK_QC_workflow.ipynb"
EXACT_SUMMARY = HERE / "ambient_transfer" / "lellouch2019_exact_stack" / "aggregate_2024-12-20_src23_ram0p1_cross_correlation_ordered_r0.json"
NOTEBOOK_VERSION = "v10" if EXACT_SUMMARY.is_file() else "v7"
nb = nbf.v4.new_notebook()
nb.metadata = {"kernelspec": {"display_name": "Python [conda env:das]", "language": "python", "name": "das"}}

def markdown(text):
    nb.cells.append(nbf.v4.new_markdown_cell(text))

def python(text):
    nb.cells.append(nbf.v4.new_code_cell(text))

markdown(f"""# Ambient-noise Figure 7c reproduction and F–K QC — {NOTEBOOK_VERSION}

This is the single advisor-facing notebook for the ambient-noise analysis. It follows one linear chain:

1. reproduce the unfiltered Lellouch et al. (2019) Figure 7c calculation from raw 2024 data;
2. test the measured gather against receiver-order, broadband-white-noise, and channel-scramble controls;
3. show the raw → normalized → correlation → F–K progression as a separate project extension;
4. compare velocity bands and signed coordinate-direction branches; and
5. place the available 5-hour, one-day, and multi-day products on a stack-duration axis.

The notebook is a decision document, not a gallery. A processing choice passes only if the measured result survives controls using the same segmentation, receiver stack, display band, and velocity-selection rule.

**v8 result.** The paper-faithful one-day rerun is complete. It uses all 1,440 files and 5,759 continuous 30-s windows, the literal unshifted R±10 sum, and no F–K filter. The measured baseline peaks at 5.85 km s⁻¹ rather than 3.2 km s⁻¹ and does not exceed the receiver-order scan-max null (p = 0.147). Its gather remains almost unchanged after receiver channels are scrambled before preprocessing (flattened correlation 0.9976), so the dominant structure is not ordered moveout. This rejects a Figure 7c reproduction for the matched 20 December 2024 day; it does not invalidate F–K filtering as a standard directional tool or rule out other dates and longer stacks.
""")

markdown(f"""{NOTEBOOK_VERSION}""")

markdown("""## Current exact-reproduction status

**Authoritative status: completed for the matched one-day test; Figure 7c is not reproduced.** The raw-to-chunk operator is `ambient_lellouch2019_exact_stack.py`; the relevant audited commits are `9c70083` (paper operator), `cc07080` (all-channel common-mode and post-average Equation 6 sensitivities), and `bb48942` (filter the complete correlation before lag cropping). Chunk array `38988141` and final v4 aggregate array `38993456` completed with zero exit codes.

The baseline was fixed before its result was examined: provisional wellhead channel 23; receiver centers at 50–700 m offsets; the same source paired with R−10 through R+10; time differentiation; centered running-absolute-mean normalization; continuous 30-s windows every 15 s; ordinary cross-spectrum averaging; an unshifted R±10 sum; full-correlation 5–20 Hz filtering followed by lag cropping; and exactly 5,759 unique windows. There is no F–K filter, linear detrend, common-mode subtraction, input bandpass, spectral whitening, per-window correlation normalization, or imported 3.2-km-s⁻¹ alignment in this baseline.

| Branch | Why it exists | Status relative to the paper |
|---|---|---|
| source 23, RAM 0.1 s, ordinary correlation | primary Figure 7c reproduction | baseline; RAM duration remains an explicit ambiguity |
| source 0 | tests the current channel-origin assumption | coordinate sensitivity |
| RAM 5 s | tests the former project value | preprocessing sensitivity |
| 900-channel median common-mode subtraction | evaluates the quoted opinion and zero-lag concern | unreported sensitivity |
| post-average stabilized source-power division | evaluates one operational interpretation of Equation 6 | underspecified sensitivity |
| iid broadband white noise | verifies that the complete operator does not create significant ordered moveout | full-pipeline null |
| receiver channels permuted before preprocessing | destroys spatial ordering while retaining measured temporal content | full-pipeline null |

Every branch contains 1,440 files and 5,759 windows. None clears the familywise receiver-order test at α = 0.05. The baseline and channel-scramble gathers are 0.9976 correlated, while all-channel common-mode subtraction removes the dominant structure but does not reveal a significant 3.2-km-s⁻¹ ridge. Everything below an explicit **legacy/withdrawn** banner is retained only as audit history and does not determine this conclusion.
""")

python("""from pathlib import Path
import json
from IPython.display import display, Image
exact_dir = ROOT/'ambient_transfer'/'lellouch2019_exact_stack' if 'ROOT' in globals() else Path.cwd()/'ambient_transfer'/'lellouch2019_exact_stack'
if not exact_dir.exists() and Path.cwd().name != 'awd_clean':
    exact_dir = Path.cwd()/'awd_clean'/'ambient_transfer'/'lellouch2019_exact_stack'
validation_json = exact_dir/'synthetic_validation.json'
validation_png = exact_dir/'synthetic_validation.png'
assert validation_json.is_file() and validation_png.is_file(), (validation_json, validation_png)
display(Image(filename=str(validation_png), width=1100))
print(json.dumps(json.loads(validation_json.read_text()), indent=2))
""")

markdown("""**Figure QC-1. Independent validation of the exact Figure 7c operator before measured-data interpretation.** (a) Power spectral density of two deterministic-seed, independent Gaussian raw-input channels spans the complete 0–50 Hz Nyquist interval of this reduced synthetic test; the nearly zero fitted log-power slope demonstrates that the control is broadband rather than pre-shaped into the 5–20 Hz analysis band. (b) Their normalized cross-correlation remains near zero over ±1 s, and the zero-lag Pearson coefficient is −0.00179, demonstrating that the chosen noise channels do not already share a coherent waveform. (c) A plane wave with known 500 m s⁻¹ apparent velocity is propagated toward increasing channel coordinate and passed through differentiation, centered RAM, 30-s/15-s segmentation, cross-spectrum accumulation, R±10 summation, full-correlation bandpass filtering, and lag cropping. The recovered median lag error is exactly 0 at the 0.01-s sampling interval. Machine-checked gates additionally recover 5,759 unique full-day windows, verify source-versus-summed-receiver equivalence to 3.3×10⁻⁸ relative error, and reproduce an independently constructed full-filter-then-crop reference with zero relative error. These tests validate implementation sign, segmentation, input-noise construction, R±10 linearity, and filter order; they do not assert that the measured archive contains the published P-wave arrival.""")

markdown("""## Read this first: exactly which data are used

**This notebook does not load the June 2026 AWD files from `ActiveJune2026/Nano/`.** It analyzes the separate 2024–2025 continuous recording from the cemented fiber in the SAFOD main hole. AWD enters only once, in Figure 5, as an independent propagation-direction check; it is not included in any ambient-noise correlation stack.

The measured-data loading chain is

`SAFOD_2024_2025.csv` → corrected mounted HDF5 path → `Acquisition/Raw[0]/RawData` → channel × time array → preprocessing/correlation → frozen NPZ product → notebook figure.

The archive index is

`/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv`.

Rows in that transferred index contain the older prefix `/oak/stanford/groups/ettore88/data/SAFODAS1-harddrive-transfer/`. The loader inserts the missing `/SAFOD/` directory component and checks that the resulting HDF5 file exists. This is a path correction only; it does not alter timestamps or data values.

The primary development day is **20 December 2024 UTC**. It contains 1,440 one-minute files from 00:00:06 to 23:59:06 UTC. The five-hour matched and surrogate tests use the first 300 files, from 00:00:06 through 04:59:06 UTC. Seven other complete days are held out for seasonal validation. Exact file names, header values, dates, figure inputs, and derived-product paths are printed in Section 0 below.

**Raw versus derived inputs.** HDF5 files are the measured DAS data. NPZ and JSON files under `ambient_transfer/` are frozen intermediate results computed from those HDF5 files by the named scripts. They are loaded here so an advisor can inspect the full QC argument without recomputing eight complete days and 20 surrogate realizations during notebook execution.
""")

markdown("""## Frozen analysis contract

The corrected Figure 7c baseline uses provisional wellhead channel 23 as the virtual source and repeats the complete calculation with channel 0 as a coordinate sensitivity. Receiver centers are placed at approximately 50-m increments from 50 to 700 m and each output is the literal sum of correlations to R−10 through R+10 using the same source. Phase/strain is differentiated continuously, normalized with a centered running absolute mean, divided into 30-s windows every 15 s across file boundaries, correlated without per-window amplitude normalization, simply stacked, and filtered at 5–20 Hz only after correlation. Figure 7c contains no travel-time shifts. Earthquake-containing records remain in the archive; temporal normalization suppresses their amplitude leverage but does not remove their waveforms.

Lellouch et al. (2019) used F–K filtering for earthquake S-wave isolation, not in the ambient-interferometry subsection. Our ambient F–K work is therefore an explicitly labelled extension. Channel coordinate increases downhole. A synthetic check maps F×K<0, evaluated at positive lag, to increasing-coordinate propagation and F×K>0, evaluated at negative lag, to decreasing-coordinate propagation. The known surface-AWD arrival independently validates that physical mapping at 25–60 Hz; its 5–20 Hz energy is nearly directionally balanced, so the ambient-band physical labels retain that caveat.
""")

markdown("""## Published basis and project-specific extensions

The workflow deliberately separates what is reported in the literature from what is being tested here.

| Choice | Basis | Status in this notebook |
|---|---|---|
| Differentiate DAS strain in time | [Lellouch et al. (2019)](https://doi.org/10.1029/2019JB017533) convert the 2017 SAFOD strain records to strain rate. | Reproduced. |
| Correlate short, overlapping windows and stack them | Lellouch et al. use 30-s windows with 15-s overlap at SAFOD; [Seats et al. (2012)](https://doi.org/10.1111/j.1365-246X.2011.05263.x) analyze why overlapping-window correlation can improve ambient-noise convergence. | Reproduced, then convergence-tested. |
| Reduce transient amplitude leverage before correlation | Lellouch et al. use running-absolute-mean normalization; [Bensen et al. (2007)](https://doi.org/10.1111/j.1365-246X.2007.03374.x) describe temporal normalization as a way to suppress earthquake and instrumental transients. | Reproduced in form. The baseline uses **0.1 s**, following Bensen et al.'s recommendation of approximately half the maximum 0.2-s period in the 5–20 Hz band; the former 5-s choice is rerun as a sensitivity. Neither duration is reported by Lellouch et al. |
| Interpret a stacked cross-correlation as an inter-receiver response | [Shapiro and Campillo (2004)](https://doi.org/10.1029/2004GL019491) demonstrate coherent wave emergence from ambient-noise correlations; [Wapenaar and Fokkema (2006)](https://doi.org/10.1190/1.2213955) derive the source-boundary and medium assumptions behind Green's-function retrieval. | A hypothesis to test, not assumed from appearance alone. |
| Separate signed apparent-slowness branches in F–K space | Frequency–wavenumber analysis estimates slowness and propagation direction ([Rost and Thomas, 2002](https://doi.org/10.1029/2000RG000100)); VSP literature uses directional filtering to separate upgoing and downgoing wavefields ([Suprajitno and Greenhalgh, 1985](https://doi.org/10.1190/1.1441973)). | Valid directional operator after synthetic sign calibration. |
| Use F–K or coherence filtering to denoise DAS | [Isken et al. (2022)](https://doi.org/10.1093/gji/ggac229) apply adaptive F–K denoising to borehole DAS and explicitly note distortion of relative waveform amplitudes. [Ehsaninezhad et al. (2024)](https://doi.org/10.1093/gji/ggae134) compare conventional DAS ambient interferometry with a coherence-enhanced workflow and evaluate stack convergence. | Supports before/after controls and cautious amplitude interpretation; neither study validates treating a hard, pre-correlation velocity fan as independent evidence for that velocity. |
| Apply a 2.5–4.5 km/s F–K fan **before ambient correlation** | This is **not** the ambient workflow reported by Lellouch et al.; their F–K filter is used in the earthquake/VSP analysis. | Project-specific extension requiring white-noise, channel-scramble, and time-shift falsification tests. |
| Estimate the minimum stack duration empirically | Bensen et al. emphasize temporal repeatability as QC; Seats et al. compare short stacks with a long-term reference. | Project-specific combinatorial convergence test, reported separately for unfiltered and fan-filtered sections. |

The literature therefore justifies the operations and the need for QC; it does **not** validate this dataset's 3.2 km/s ridge. That decision comes only from the real-data controls and surrogate tests below.
""")

markdown("""### Lellouch reproduction fidelity audit

| Item | Lellouch et al. (2019) | This notebook | Consequence |
|---|---|---|---|
| Physical array | Cemented fiber in the SAFOD main hole | Same physical fiber; source channel 23 is the provisional G0 wellhead and channel 0 is a predeclared sensitivity | Geometry is explicit; absolute top-channel registration remains uncertain at the tens-of-metres level. |
| Recording epoch | June–July 2017 | 2024–2025 archive; matched example is 20 December 2024 | The ambient source field and instrument response can differ. |
| Interrogator sampling | 2,500 Hz, approximately 1-m channel spacing, 10-m gauge length | 500 Hz and 1.020952-m header spacing; current gauge/interrogator metadata must be taken from the 2024 acquisition records | The method is reproduced, not the acquisition hardware. |
| DAS observable | Strain differentiated in time to strain rate | Same | Direct reproduction. |
| Temporal normalization | Running absolute mean; duration unreported | Centered 0.1-s baseline plus 5-s sensitivity | The operation is reproduced; the duration dependence is measured rather than hidden. |
| Correlation segmentation | 30-s windows, 15-s overlap | Same | Direct reproduction. |
| Virtual-source geometries | Fixed top source for Figure 7c; separate 50-m constant-offset Figure 7d | The active goal and exact rerun target Figure 7c only | Figure 7d is a later velocity-profile stage and is not required to decide whether Figure 7c reproduces. |
| Nearby-receiver stack | Figure 7c uses the simple R−10:R+10 sum; shifts are introduced only after Figure 7c supplies an average velocity | The exact baseline computes only the unshifted Figure 7c sum; the old imported-velocity panel is legacy | This prevents circularly forcing 3.2 km/s into the reproduction target. |
| Display band | 5–20 Hz after correlation stack | Same | Direct reproduction. |
| F–K use | Reported for earthquake/VSP processing, not the ambient-interferometry subsection | Applied before ambient correlation as the tested extension | All F–K conclusions require matched controls and input-level nulls. |

The source for the 2017 acquisition and processing entries is [Lellouch et al. (2019)](https://doi.org/10.1029/2019JB017533). Values for the present archive are read from the current HDF5 headers and frozen products.
""")

markdown("""### Audit of the three disputed preprocessing claims

1. **Common-mode removal — not a published ambient-interferometry step.** Section 4.1 reports time differentiation, running-absolute-mean normalization, correlation windows, nearby-pair stacking, travel-time shifting, and final bandpass filtering. It does not report subtracting the across-channel median. That operation existed in a legacy local pipeline, so it is tested as a sensitivity branch rather than inserted into the reproduction baseline.

2. **Equation 6 source-spectrum division — theoretically present but operationally underspecified.** The paper writes the Green-function proportionality with division by the ambient-field power spectrum. It does not state how that spectrum was estimated, averaged, water-level stabilized, or applied window by window. The baseline therefore preserves the explicitly described cross-correlation stack; a separately labelled stabilized source-power branch tests whether this ambiguity controls the result.

3. **R±10 nearby-receiver stack — genuinely required and previously missing from the legacy null path.** The paper explicitly sums 21 correlations for each output receiver and states that this smoothing is required to extract the signal. The corrected operator performs that literal sum in the fixed-top Figure 7c geometry without travel-time shifts. Twenty-one traces do not imply a 21-fold SNR increase: under independent-noise assumptions the ideal amplitude-SNR gain is at most the square root of 21, approximately 4.6, and is smaller for correlated neighboring noise.

Thus the quoted opinion was partly correct: item 3 identified a decisive omission in the old single-pair test; item 1 incorrectly promoted a legacy project choice to a published requirement; and item 2 identified a real fidelity ambiguity but overstated how completely the paper specifies its implementation.
""")

python("""from pathlib import Path
import json, numpy as np, matplotlib.pyplot as plt
ROOT = Path.cwd()
if ROOT.name != 'awd_clean': ROOT = ROOT / 'awd_clean'
OUT = ROOT / 'ambient_transfer'
plt.rcParams.update({'figure.dpi': 120, 'axes.grid': True, 'grid.alpha': .2})
products = {
 'lellouch': OUT/'lellouch2019_reproduction_v1/lellouch2019_2024-12-20_start0_requestedall_used1440.npz',
 'lellouch_meta': OUT/'lellouch2019_reproduction_v1/lellouch2019_2024-12-20_start0_requestedall_used1440.json',
 'unfiltered_meta': OUT/'seasonal_unfiltered_aggregate.json',
 'fk_grid': OUT/'fk_mask_sensitivity_v2/ambient_fk_mask_sensitivity_v2.npz',
 'fk_grid_meta': OUT/'fk_mask_sensitivity_v2/ambient_fk_mask_sensitivity_v2.json',
 'fk_nulls': OUT/'fk_full_pipeline_null_v2_n300_r20/fk_full_pipeline_null_v2_aggregate.npz',
 'fk_null_example': OUT/'fk_full_pipeline_null_v2_n300_r20/fk_full_pipeline_null_v2_2024-12-20_start0_n300_r0-0.npz',
 'fk_null_example_meta': OUT/'fk_full_pipeline_null_v2_n300_r20/fk_full_pipeline_null_v2_2024-12-20_start0_n300_r0-0.json',
 'directional_meta': OUT/'fk_directional_audit_v1/ambient_fk_directional_audit_v1.json',
 'white_noise': OUT/'fk_qc_notebook_v2/ambient_fk_white_noise_v1.npz',
 'white_noise_meta': OUT/'fk_qc_notebook_v2/ambient_fk_white_noise_v1.json',
 'convergence': OUT/'fk_qc_notebook_v2/ambient_fk_convergence_v1.npz',
 'convergence_meta': OUT/'fk_qc_notebook_v2/ambient_fk_convergence_v1.json',
 'alias': OUT/'alias_sensitivity_2024-12-20_start0_n300.npz',
 'alias_meta': OUT/'alias_sensitivity_2024-12-20_start0_n300.json',
 'sign_figure': ROOT/'fk_sign_synthetic_test.png',
 'sign_meta': ROOT/'fk_sign_synthetic_test.json',
 'awd_sign_figure': OUT/'awd_production_operator_validation_v1/awd_fk_production_operator_validation_v1.png',
 'awd_sign_meta': OUT/'awd_production_operator_validation_v1/awd_fk_production_operator_validation_v1.json',
 'matched': OUT/'lellouch_fk_matched_v2/lellouch_fk_matched_2024-12-20_start0_n300.npz',
 'matched_meta': OUT/'lellouch_fk_matched_v2/lellouch_fk_matched_2024-12-20_start0_n300.json',
}
for name, path in products.items(): print(f'{name:16s}', 'OK' if path.exists() else 'MISSING', path)
assert all(path.exists() for path in products.values())
""")

markdown("""## 0. Data provenance, raw-file audit, and checkpoint map

This notebook uses three input classes that must not be conflated:

1. **Measured ambient DAS:** one-minute HDF5 records from the 2024–2025 continuous archive on the cemented SAFOD main-hole fiber.
2. **Reproducible checkpoints:** NPZ arrays plus JSON metadata written by the named analysis scripts after reading those HDF5 records. They prevent a two-week stack from being recomputed during every notebook review.
3. **Controls:** deterministic Gaussian white noise, exact analytic plane waves, and a separate June 2026 AWD product used only for the propagation-direction check.

The cell below opens the exact raw file used in the five-hour analysis, reads its acquisition header, and prints a figure-by-figure ledger. `RawData` is stored on disk as **time sample × channel/locus** and each loader transposes it to **channel × time**. Distance is along-fiber coordinate, not independently surveyed true vertical or measured depth. The exact baseline uses channel 23 as the provisional wellhead from the independent G0 noise transition and reruns channel 0 as a predeclared coordinate sensitivity; neither is silently treated as surveyed absolute depth.

The expensive checkpoints are not unaudited black boxes. The notebook prints their producing script/command, dates, duration, and geometry. Section 0a then recomputes a small raw subset end to end so that the path from HDF5 samples to a correlation gather can be inspected without waiting for all long-duration jobs.
""")

python("""import sys, h5py, pandas as pd
from IPython.display import display
sys.path.insert(0, str(ROOT))
from ambient_transfer_test import CSV, corrected_path

archive_index = pd.read_csv(CSV, sep=r'\\s+')
archive_index = archive_index[archive_index.nSamples > 0].copy()
archive_index['time_utc'] = pd.to_datetime(archive_index.startTime, utc=True, errors='coerce')
development_rows = archive_index[archive_index.time_utc.dt.strftime('%Y-%m-%d') == '2024-12-20'].sort_values('time_utc').reset_index(drop=True)
null_example_manifest = np.load(products['fk_null_example'])
raw_example_path = Path(str(null_example_manifest['used_files'][0]))

def decoded(value):
    return value.decode() if isinstance(value, (bytes, np.bytes_)) else value

with h5py.File(raw_example_path, 'r') as handle:
    acquisition = handle['Acquisition']
    raw_group = handle['Acquisition/Raw[0]']
    raw_dataset = raw_group['RawData']
    sample_rate = float(raw_group.attrs['OutputDataRate'])
    header = {
        'HDF5 dataset key': raw_dataset.name,
        'On-disk array order': 'time sample × channel/locus',
        'On-disk shape': str(tuple(raw_dataset.shape)),
        'Stored dtype': str(raw_dataset.dtype),
        'Sample rate': f'{sample_rate:.1f} Hz',
        'Samples and duration per file': f'{raw_dataset.shape[0]:,}; {raw_dataset.shape[0] / sample_rate:.1f} s',
        'Number of loci/channels': f"{int(acquisition.attrs['NumberOfLoci'])}",
        'Along-fiber spacing': f"{float(acquisition.attrs['SpatialSamplingInterval']):.9f} m",
        'Gauge length': f"{float(acquisition.attrs['GaugeLength']):.6f} m",
        'Raw unit in header': str(decoded(raw_group.attrs.get('RawDataUnit', 'not recorded'))),
        'StartLocusIndex': f"{int(acquisition.attrs.get('StartLocusIndex', 0))}",
    }

display(pd.DataFrame({'Acquisition quantity': list(header), 'Value read from exact HDF5 header': list(header.values())}))
print('Archive index:', CSV)
print('Exact measured example:', raw_example_path)
print('December 20 usable one-minute files:', len(development_rows))
print('Full-day first UTC/file:', development_rows.iloc[0].time_utc, corrected_path(development_rows.iloc[0].file))
print('Five-hour last UTC/file:', development_rows.iloc[299].time_utc, corrected_path(development_rows.iloc[299].file))
print('Full-day last UTC/file:', development_rows.iloc[-1].time_utc, corrected_path(development_rows.iloc[-1].file))

figure_inputs = pd.DataFrame([
    ['1', 'measured HDF5 → checkpoint', '2024-12-20; all 1,440 one-minute files', 'ch 0; targets ch 49–784 with R±10', 'ambient_lellouch2019_reproduction_v1.py --date 2024-12-20 --nfiles all', products['lellouch']],
    ['2', 'measured raw HDF5', 'first 12 s of exact 2024-12-20 00:00:06 file printed above', 'ch 0–799', 'opened directly with h5py', raw_example_path],
    ['3–4', 'measured HDF5 → checkpoints', 'seven held-out complete days in 2024–2025', 'channel-0 geometry; equal day weights', 'ambient_fk_mask_sensitivity_v2.py + aggregate_ambient_fk_mask_sensitivity_v2.py', products['fk_grid']],
    ['5', 'controls only', 'analytic waves + separate June 2026 AWD product', 'no ambient stack contribution', 'sign-test scripts', products['awd_sign_figure']],
    ['6', 'measured HDF5 → checkpoint', '2024-12-20; first 300 files = 5 h', 'same geometry for every branch', 'ambient_lellouch_fk_matched_v2.py --date 2024-12-20 --start 0 --nfiles 300', products['matched']],
    ['7', 'measured HDF5 → checkpoint', 'same first 300 files', 'resampling method alone changes', 'ambient_alias_sensitivity.py --date 2024-12-20 --start 0 --nfiles 300', products['alias']],
    ['8a–d', 'generated control + measured checkpoints', 'seed 20260813; 800 ch; 500 Hz; 12 s records', 'independent N(0,1) samples before preprocessing', 'ambient_fk_white_noise_v1.py and in-notebook raw audit', products['white_noise']],
    ['9a–c', 'measured HDF5 + surrogates', 'first 300 December 20 files; 20 realizations', 'pre-F–K channel permutation and time shifts', 'ambient_fk_full_pipeline_null_v2.py --date 2024-12-20 --nfiles 300 --nulls 20', products['fk_nulls']],
    ['10', 'measured HDF5 → daily checkpoints', 'eight complete dates spanning 2024-05-11 to 2025-03-04', 'all combinations of daily stacks', 'ambient_fk_convergence_v1.py', products['convergence']],
], columns=['Figure', 'Input class', 'Dates/files/duration', 'Geometry/control', 'Producing code or command', 'Checkpoint loaded here'])
display(figure_inputs)

selection_path = OUT/'seasonal_day_selection.json'
selection = json.loads(selection_path.read_text())
seasonal_rows=[]
for item in selection['days']:
    chosen = archive_index[archive_index.time_utc.dt.strftime('%Y-%m-%d') == item['date']].sort_values('time_utc').reset_index(drop=True)
    assert len(chosen) == int(item['nfiles']), (item['date'],len(chosen),item['nfiles'])
    seasonal_rows.append([
        item['date'],item['season'],len(chosen),chosen.iloc[0].time_utc,chosen.iloc[-1].time_utc,
        corrected_path(chosen.iloc[0].file),corrected_path(chosen.iloc[-1].file),
    ])
print('Frozen seasonal selection:',selection_path)
print('Selection rule:',selection['criteria'])
display(pd.DataFrame(seasonal_rows,columns=['UTC date','season','raw HDF5 count','first record UTC','last record UTC','first raw HDF5','last raw HDF5']))

lellouch_check=json.loads(products['lellouch_meta'].read_text())
matched_check=json.loads(products['matched_meta'].read_text())
alias_check=json.loads(products['alias_meta'].read_text())
fk_grid_check=json.loads(products['fk_grid_meta'].read_text())
null_check=json.loads(products['fk_nulls'].with_suffix('.json').read_text())
convergence_check=json.loads(products['convergence_meta'].read_text())
checkpoint_checks=pd.DataFrame([
    ['Figure 1 full day',1440,lellouch_check['used_files'],lellouch_check['date']=='2024-12-20'],
    ['Figure 6 matched five hours',300,matched_check['used_files'],matched_check['date']=='2024-12-20' and matched_check['start']==0],
    ['Figure 7 alias test',300,alias_check['used_files'],alias_check['date']=='2024-12-20' and alias_check['start']==0],
    ['Figures 3–4 seasonal masks',sum(item['nfiles'] for item in selection['days']),sum(v['used_files_by_mask']['production_2p5_4p5'] for v in fk_grid_check['day_completeness'].values()),all(not v['missing_chunks'] for v in fk_grid_check['day_completeness'].values())],
    ['Figures 9b–9c measured nulls',300,null_check['used_files'],null_check['date']=='2024-12-20' and null_check['start']==0 and null_check['null_realizations']==20],
    ['Figure 10 daily convergence',8,len(convergence_check['dates']),set(convergence_check['dates'])==set(item['date'] for item in selection['days'])],
],columns=['Checkpoint','Expected raw files or days','Recorded files or days','Metadata/date/complete check'])
checkpoint_checks['PASS']=checkpoint_checks['Expected raw files or days'].eq(checkpoint_checks['Recorded files or days']) & checkpoint_checks['Metadata/date/complete check']
display(checkpoint_checks)
assert bool(checkpoint_checks.PASS.all()), checkpoint_checks
print('June 2026 AWD files included in an ambient correlation stack: NO')
""")

markdown("""**Data-provenance reading guide.** The header table is read from the exact first HDF5 record used by the five-hour matched and measured-surrogate analyses, rather than copied from a generic instrument specification. The ledger identifies whether each panel comes from raw measured data, a reproducible checkpoint, or a synthetic/independent control. The June 2026 AWD record appears only in Figure 5 and cannot create an ambient-noise correlation ridge elsewhere in this notebook.""")

markdown("""### 0a. Legacy two-file rebuild — withdrawn from the exact decision

**Withdrawal reason.** This historical cell uses `ambient_lellouch2019_reproduction_v1.py`, which linearly detrends and resets differentiation/RAM at each one-minute file. It is retained only to expose the earlier data path. It does not implement continuous preprocessing and its figure must not be used for the Figure 7c decision. The authoritative code and synthetic validation are in the opening v6 section.
""")

python("""from ambient_lellouch2019_reproduction_v1 import (
    acquisition_metadata, geometry, load_required_channels,
    strain_rate_and_ram, correlation_window, align_and_stack, final_bandpass,
    MAX_LAG_SECONDS, WINDOW_SECONDS, STEP_SECONDS,
)
RAW_AUDIT_FILES = 2
audit_paths = [corrected_path(path) for path in development_rows.file.iloc[:RAW_AUDIT_FILES]]
audit_fs, audit_dx, audit_n_channels, _ = acquisition_metadata(audit_paths[0])
audit_targets, audit_required = geometry(audit_dx, audit_n_channels)
audit_source_index = int(np.flatnonzero(audit_required == 0)[0])
audit_sum = None; audit_windows = 0; audit_lags = None
for audit_path in audit_paths:
    audit_raw, current_fs = load_required_channels(audit_path, audit_required)
    assert current_fs == audit_fs
    audit_processed = strain_rate_and_ram(audit_raw, audit_fs)
    window_samples = int(round(WINDOW_SECONDS * audit_fs))
    step_samples = int(round(STEP_SECONDS * audit_fs))
    for start_sample in range(0, audit_processed.shape[1] - window_samples + 1, step_samples):
        current_lags, current_corr = correlation_window(
            audit_processed[:, start_sample:start_sample + window_samples],
            audit_source_index, audit_fs, MAX_LAG_SECONDS,
        )
        audit_sum = current_corr if audit_sum is None else audit_sum + current_corr
        audit_lags = current_lags
        audit_windows += 1
audit_mean = audit_sum / audit_windows
_, audit_simple, audit_aligned = align_and_stack(
    audit_mean, audit_lags, audit_required, audit_targets, audit_dx,
)
audit_distances = audit_targets.astype(float) * audit_dx
audit_simple = final_bandpass(audit_simple, audit_fs)
audit_aligned = final_bandpass(audit_aligned, audit_fs)
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True, sharey=True, constrained_layout=True)
audit_limit = np.nanpercentile(np.abs(np.concatenate([audit_simple.ravel(), audit_aligned.ravel()])), 98.5)
for ax, section, title in zip(axes, [audit_simple, audit_aligned], ['Raw audit: simple R±10 stack', 'Raw audit: locally aligned at 3.2 km/s']):
    ax.imshow(section, extent=[audit_lags[0], audit_lags[-1], audit_distances[-1], audit_distances[0]], aspect='auto', cmap='RdBu_r', vmin=-audit_limit, vmax=audit_limit)
    ax.plot(audit_distances/3200, audit_distances, 'k--', lw=1)
    ax.set(title=title, xlabel='Correlation lag (s)', ylabel='Receiver offset from channel 0 (m)')
plt.show()
print('Raw audit files:')
for item in audit_paths: print(' ', item)
print('Raw audit 30-s windows:', audit_windows)
print('This cell read HDF5 and recomputed these sections; it did not load an NPZ stack.')
""")

markdown("""**Legacy raw-audit figure — not the exact operator. End-to-end HDF5-to-correlation reconstruction on two contiguous one-minute records.** The left panel is the simple R±10 receiver stack and the right panel applies the separately declared 3.2-km/s local alignment. Both are rebuilt in this notebook from the exact raw files printed beneath the figure, using the same strain-rate proxy, 5-s running-absolute-mean normalization, 30-s windows with 15-s overlap, channel-0 virtual source, receiver geometry, and final 5–20 Hz bandpass as Figure 1. Two minutes are intentionally insufficient for a scientific ambient-noise conclusion; this is a transparent executable provenance check showing that the cached long-duration workflow starts from the stated HDF5 dataset and produces the expected array shapes and coordinates.""")

markdown("""## 1. Legacy unfiltered attempt — withdrawn from the exact decision

This was the prior attempted baseline. It has no F–K filter, but it is **not** the exact control because preprocessing resets at file boundaries, includes an unreported detrend, uses the former 5-s RAM choice, and shows an imported 3.2-km/s alignment. It remains visible only for provenance. The dashed lines are comparison trajectories; they are not fitted velocities. The geometry, windowing, nearby-receiver stack, final 5–20 Hz filter, and post-stack 3.2 km/s alignment follow the ambient-interferometry description of [Lellouch et al. (2019)](https://doi.org/10.1029/2019JB017533). The data epoch and interrogator are different, and the present channel-0 registration is assumed to correspond to the same top coordinate on the same cemented main-hole fiber; those are experiment-specific caveats, not facts supplied by the 2019 paper.
""")
python("""meta = json.loads(products['lellouch_meta'].read_text())
z = np.load(products['lellouch']); lags = z['lags_s']; dist = z['distances_m']
simple=z['local_21_channel_simple']; aligned=z['local_21_channel_aligned_3200']
fig, axes = plt.subplots(1,2,figsize=(13,5),sharex=True,sharey=True,constrained_layout=True)
lim = np.nanpercentile(np.abs(np.concatenate((simple.ravel(),aligned.ravel()))), 98.5)
for ax,section,title in zip(axes,[simple,aligned],['Simple R±10 receiver stack','R±10 aligned locally at 3.2 km/s']):
    ax.imshow(section, extent=[lags[0], lags[-1], dist[-1], dist[0]], aspect='auto', cmap='RdBu_r', vmin=-lim, vmax=lim)
    for v, style, color in [(2500, '--', '.4'), (3200, '-', 'k'), (4000, '--', '.4')]:
        ax.plot(dist/v, dist, style, color=color, lw=1.2, label=f'{v/1000:g} km/s')
    ax.set(xlabel='Correlation lag (s)', ylabel='Receiver offset from channel 0 (m)', title=title)
axes[0].legend(); plt.show()
print('Used files:', meta['used_files'], '30-s windows:', meta['used_30_s_windows'])
print('F-K status:', meta['important_departures_or_boundaries'][-1])
""")

markdown("""**Figure 1 (legacy/withdrawn). Unfiltered fixed-top ambient-noise correlations and a published-velocity control.** Correlations between channel 0 and receivers at approximately 50-m increments are computed from 1,440 one-minute records using a strain-rate proxy, running-absolute-mean temporal normalization, 30-s windows with 15-s overlap, and a final 5–20 Hz bandpass. The left panel sums each target correlation with correlations to its R±10 neighboring receivers without time shifts; this is the first published stacking pass. The right panel shifts those neighboring traces using 3.2 km/s and then stacks them. In this frozen product, 3.2 km/s is imported from Lellouch et al. rather than estimated from the present left panel. The right panel is therefore a velocity-conditioned control, not the second step of a completed data-driven reproduction. The separate 50-m constant-offset geometry and its quadratic travel-time picks are not contained in this figure.""")

markdown("""## 2. Raw → strain rate → temporal normalization → F–K → correlation

This cell follows one real record through the operations in their actual order. [Lellouch et al. (2019)](https://doi.org/10.1029/2019JB017533) convert strain to a strain-rate proxy, use running-absolute-mean normalization, divide the continuous record into 30-s windows with 15-s overlap, correlate channel 0 with receivers, stack nearby receivers, and finally filter the correlations from 5–20 Hz. Temporal normalization suppresses the leverage of earthquakes and other transients; it does not remove those waveforms ([Bensen et al., 2007](https://doi.org/10.1111/j.1365-246X.2007.03374.x)). The F–K extension is inserted after normalization and before correlation, so its effect must be tested against the unfiltered path with all other operations held fixed.
""")
python("""try:
    import sys, h5py, pandas as pd
    sys.path.insert(0, str(ROOT))
    from scipy.ndimage import uniform_filter1d
    from scipy.signal import detrend
    path = raw_example_path
    if not path.exists(): raise FileNotFoundError(f'Exact shared HDF5 record not mounted: {path}')
    with h5py.File(path,'r') as handle:
        raw_group=handle['Acquisition/Raw[0]']; dataset=raw_group['RawData']
        fs=float(raw_group.attrs.get('OutputDataRate',dataset.attrs.get('OutputDataRate',500.0))); dx=float(handle['Acquisition'].attrs.get('SpatialSamplingInterval',1.0))
        excerpt_samples=min(int(12*fs),dataset.shape[0])
        raw=np.asarray(dataset[:excerpt_samples,:800],dtype='float32').T
    strain = detrend(raw, axis=1, type='linear').astype('float32')
    rate = np.empty_like(strain); rate[:,0]=0; rate[:,1:]=np.diff(strain,axis=1)*fs
    ram = uniform_filter1d(np.abs(rate),size=int(5*fs),axis=1,mode='nearest')
    normalized = rate/np.maximum(ram,np.percentile(ram,5,axis=1,keepdims=True)*.1+1e-12)
    dec = normalized[::2,::2]; f=np.fft.fftfreq(dec.shape[1],1/(fs/2)); k=np.fft.fftfreq(dec.shape[0],dx*2)
    K,F=np.meshgrid(k,f,indexing='ij'); velocity=np.abs(F)/np.maximum(np.abs(K),1e-12)
    mask=(np.abs(F)>=5)&(np.abs(F)<=20)&(velocity>=2500)&(velocity<=4500)&(np.abs(K)>0)&(F*K<0)
    filtered=np.fft.ifft2(np.fft.fft2(dec)*mask).real
    samples=slice(0,min(int(12*fs),raw.shape[1])); dsamples=slice(0,min(int(12*fs/2),filtered.shape[1]))
    fig, ax = plt.subplots(1, 5, figsize=(18, 4), constrained_layout=True)
    panels=[(strain[:400,samples],'Raw strain/phase'),(rate[:400,samples],'Time derivative: strain rate'),(normalized[:400,samples],'5-s RAM normalization')]
    for a,(data,title) in zip(ax[:3],panels):
        limit=np.nanpercentile(np.abs(data),99)
        a.imshow(data,aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit,interpolation='nearest'); a.set_title(title); a.set(xlabel='Sample',ylabel='Channel')
    power=np.abs(np.fft.fft2(dec))**2; order=np.argsort(k)
    ax[3].pcolormesh(f,k[order]*1000,10*np.log10(power[order]/power.max()+1e-12),shading='auto',cmap='magma',vmin=-40,vmax=0)
    ax[3].set(xlim=(0,30),ylim=(-10,10),xlabel='Frequency (Hz)',ylabel='Wavenumber (cycles/km)',title='Pre-filter F–K power')
    fshown=filtered[:400,dsamples]; flimit=np.nanpercentile(np.abs(fshown),99)
    ax[4].imshow(fshown,aspect='auto',cmap='RdBu_r',vmin=-flimit,vmax=flimit,interpolation='nearest')
    ax[4].set_title('Downgoing 2.5–4.5 km/s fan'); ax[4].set(xlabel='Decimated sample',ylabel='Channel')
    plt.show(); print('Raw file:', path, 'shape:', raw.shape, 'fs:', fs, 'dx:', dx)
except Exception as exc:
    print('Raw display skipped:', repr(exc))
""")

markdown("""**Figure 2. Processing progression for a 12-s excerpt of one measured DAS record.** From left to right, the panels show the raw phase/strain-like record, its time derivative used as a strain-rate proxy, division by a 5-s running mean of absolute amplitude, the pre-filter two-dimensional frequency–wavenumber power spectrum, and the wavefield retained by the 5–20 Hz, 2.5–4.5 km/s, F×K<0 fan. The first three panels use the original temporal and spatial sampling; the final two use 2× temporal and spatial decimation. The F–K panel is displayed before the filtered wavefield to make the selection order explicit. The shorter excerpt is used only to document the operator without loading redundant samples into the notebook kernel; all quantitative products use their durations stated in the corresponding sections. This figure documents what the operator does to measured data; it does not demonstrate that the retained energy is a physical arrival.""")

markdown(r"""### Coordinate, correlation, and Fourier conventions

Let $x$ be distance along the fiber, with $x=0$ at channel 0 and increasing $x$ down the cemented main-hole fiber. The code computes the source–receiver correlation

\[
C_{0r}(\tau)=\int u_0(t)\,u_r(t+\tau)\,dt .
\]

With this convention, a receiver waveform delayed relative to channel 0 appears at positive lag. The two-dimensional FFT uses the NumPy sign convention in time and channel coordinate. An increasing-coordinate plane wave

\[
u(x,t)=\cos\left[2\pi f_0\left(t-\frac{x}{v}\right)\right]
\]

has opposite signs of temporal frequency and spatial wavenumber, so it is retained by $F K<0$. Reversing the propagation term to $t+x/v$ gives $F K>0$. Apparent velocity is

\[
v_{\mathrm{app}}=\frac{|f|}{|k|},
\]

where $f$ is in hertz and $k$ is in cycles per meter. A line drawn on a correlation gather follows

\[
\tau(\Delta x)=\frac{\Delta x}{v_{\mathrm{app}}}.
\]

These equations define **fiber-coordinate direction and apparent velocity**. They do not prove vertical propagation, true depth, ray geometry, or wave type. F–K slowness estimation and its dependence on aperture and spacing are reviewed by [Rost and Thomas (2002)](https://doi.org/10.1029/2000RG000100).
""")

markdown("""## 3. Velocity bands and the verified propagation directions

For a plane wave, signed wavenumber carries propagation-direction information; F–K analysis is a standard array method for estimating slowness and direction ([Rost and Thomas, 2002](https://doi.org/10.1029/2000RG000100)). Each row below is a frozen velocity-band choice defined on the 20 December 2024 development day and evaluated here on seven separately selected complete seasonal days with equal day weights. The left column is the increasing-coordinate branch, interpreted as downgoing because channel 0 is the top channel and channel coordinate increases downhole. The right column is the decreasing-coordinate/upgoing branch. Each branch is evaluated on its correct lag side; neither panel is a mirrored copy of the other. “Upgoing” and “downgoing” here describe **fiber-coordinate direction**, not a separately proven ray geometry or wave type. Unequal or unexpected causal/acausal energy can arise under finite, irregular, or one-sided source illumination, for which correlation does not exactly reconstruct the Green's function ([Wapenaar, 2006](https://doi.org/10.1029/2006GL027747)).
""")
python("""fk = np.load(products['fk_grid']); fk_meta=json.loads(products['fk_grid_meta'].read_text()); lags = fk['lags']; dist = fk['distance']; vel = fk['velocities_m_s']
masks = [('production_2p5_4p5','2.5–4.5 km/s'), ('narrow_2p8_3p8','2.8–3.8 km/s'), ('broad_2p0_5p5','2.0–5.5 km/s'), ('direction_only','direction only')]
print('Development day:',fk_meta['development_date']); print('Held-out days:',fk_meta['available_heldout_dates']); print('Frozen split:',fk_meta['split_is_frozen'])
fig, axes = plt.subplots(4, 2, figsize=(12, 13), sharex=True, sharey=True, constrained_layout=True)
for i, (key, label) in enumerate(masks):
    for j, branch in enumerate(['negative','positive']):
        top = fk[f'held_out__{key}__{branch}__top']; sign=1 if branch=='negative' else -1
        keep=lags>=0 if sign==1 else lags<=0; travel=sign*lags[keep]; section=top[:,keep]
        if travel[0]>travel[-1]: travel=travel[::-1]; section=section[:,::-1]
        lim=np.nanpercentile(np.abs(section),98.5)
        axes[i,j].imshow(section,extent=[travel[0],travel[-1],dist[-1],dist[0]],aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim)
        axes[i,j].plot(dist/3200,dist,'k--',lw=1)
        direction='downgoing / F×K<0' if sign==1 else 'upgoing / F×K>0'
        axes[i,j].set_title(f'{label}; {direction}')
        axes[i,j].set(xlabel='Positive travel time (s)',ylabel='Receiver offset (m)')
plt.show()
""")

markdown("""**Figure 3. Held-out sensitivity of the correlation gather to apparent-velocity and propagation-direction selection.** Rows compare the predeclared 2.5–4.5 km/s production fan, a narrower 2.8–3.8 km/s fan, a broader 2.0–5.5 km/s fan, and direction-only selection with no velocity bounds. Mask choices were frozen using 20 December 2024; the displayed sections aggregate seven other complete seasonal days with equal day weights. The left column retains F×K<0 and is evaluated at positive lag, corresponding to propagation toward increasing fiber coordinate under the calibrated convention. The right column retains F×K>0 and is evaluated at negative lag, transformed to positive travel time for direct comparison. Dashed curves show the 3.2-km/s trajectory. Each panel is independently color-scaled to its 98.5th-percentile absolute amplitude, so relative amplitudes should not be compared between panels. The loss of the ordered ridge in the direction-only row shows that direction selection alone is insufficient; the velocity fan is the operation producing the visually coherent feature.""")

python("""fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True, constrained_layout=True)
for ax, (key, label) in zip(axes.flat, masks):
    for branch, color in [('negative','tab:blue'), ('positive','tab:orange')]:
        ax.plot(vel/1000, fk[f'held_out__{key}__{branch}__scores'], color=color, label=branch)
    ax.axvline(3.2, color='k', ls='--'); ax.set_title(label); ax.set(xlabel='Trial velocity (km/s)', ylabel='Median correlation')
    ax.legend(frameon=False)
plt.show()
""")

markdown("""**Figure 4. Held-out trial-velocity scores for the frozen F–K mask grid.** Median normalized correlation sampled along constant apparent-velocity trajectories is plotted for the two signed branches and four mask choices shown in Figure 3. Curves use the seven-day held-out seasonal aggregate, not the 20 December 2024 development day. The vertical dashed line marks 3.2 km/s as a reference rather than a fitted or independently known velocity. Peaks inside a velocity-restricted fan are conditional on that fan and cannot be used as independent velocity estimates. The direction-only curves provide the least conditioned comparison and remain weak near 3.2 km/s.""")

markdown("""### 3a. Sign-convention calibration

The Fourier sign is established before interpreting the ambient panels. Exact 10-Hz plane waves traveling in opposite coordinate directions are passed through the same signed mask; each is retained only by its expected branch. A surface AWD provides an empirical direction standard because energy must initially travel from the top toward increasing fiber coordinate. That empirical test is strongly directional at 25–60 Hz but nearly balanced at 5–20 Hz. Therefore the software sign is verified, while the physical origin of both ambient 5–20 Hz branches remains unresolved.
""")
python("""from IPython.display import Image, display
display(Image(filename=str(products['sign_figure']), width=900))
display(Image(filename=str(products['awd_sign_figure']), width=1000))
sign=json.loads(products['sign_meta'].read_text()); awd_sign=json.loads(products['awd_sign_meta'].read_text())
si=sign['results']['increasing_z']; sd=sign['results']['decreasing_z']
print('Synthetic increasing-coordinate correct/wrong RMS:',si['negative_output_rms']/si['positive_output_rms'])
print('Synthetic decreasing-coordinate correct/wrong RMS:',sd['positive_output_rms']/sd['negative_output_rms'])
print('AWD 5–20 Hz expected/opposite tube-energy ratio:',awd_sign['metrics']['5_20_hz']['expected_to_opposite_energy_in_fixed_downgoing_tube_ratio'])
print('AWD 25–60 Hz expected/opposite tube-energy ratio:',awd_sign['metrics']['25_60_hz']['expected_to_opposite_energy_in_fixed_downgoing_tube_ratio'])
""")

markdown("""**Figure 5. Synthetic and empirical calibration of the signed F–K branches.** The first figure passes exact 10-Hz, 3.2-km/s plane waves traveling toward increasing and decreasing fiber coordinate through the production signed masks; the expected branch exceeds the complementary output RMS by approximately an order of magnitude in both directions, validating the Fourier implementation. The second figure applies the production operator to the known surface-AWD Nano arrival. The increasing-coordinate branch is strongly preferred at 25–60 Hz but only weakly preferred at 5–20 Hz. Together these tests validate coordinate direction and rule out a simple sign or mirroring bug. They do not determine the source, path, or wave type of the two ambient 5–20 Hz branches.""")

markdown("""## 4. Duration-matched test: is F–K itself responsible for the enhancement?

This is the decisive apples-to-apples comparison. Every path uses the same 300 one-minute records, full 500-Hz/approximately 1.02-m sampling, strain-to-strain-rate conversion, 5-s temporal normalization, 30-s windows with 15-s overlap, channel-0 virtual source, 50-m receivers, nearby-receiver stack, and final 5–20 Hz correlation filter. No path is decimated. Only the F–K mask changes. This removes the preprocessing confounds present in earlier comparisons. The unfiltered path is the closest reproduction of [Lellouch et al. (2019)](https://doi.org/10.1029/2019JB017533); the four F–K paths are controlled extensions.
""")
python("""matched_npz=OUT/'lellouch_fk_matched_v2/lellouch_fk_matched_2024-12-20_start0_n300.npz'
matched_json=matched_npz.with_suffix('.json')
assert matched_npz.exists() and matched_json.exists(), 'Matched five-hour Sherlock product is required'
matched=np.load(matched_npz); matched_meta=json.loads(matched_json.read_text())
mlags=matched['lags_s']; mdist=matched['distances_m']
modes=['unfiltered','downgoing_2p5_4p5','upgoing_2p5_4p5','downgoing_direction_only','upgoing_direction_only']
fig,axes=plt.subplots(2,5,figsize=(18,8),constrained_layout=True)
common=np.concatenate([matched[f'{mode}__simple'].ravel() for mode in modes]); limit=np.nanpercentile(np.abs(common),98.5)
trial=np.arange(2000.,5000.1,25.)
def matched_curve(section,lag_sign):
    normalized=section/np.maximum(np.max(np.abs(section),axis=1,keepdims=True),np.finfo(float).eps)
    return np.asarray([np.median([trace[np.argmin(np.abs(mlags-lag_sign*d/v))] for trace,d in zip(normalized,mdist)]) for v in trial])
for col,mode in enumerate(modes):
    ax=axes[0,col]
    section=matched[f'{mode}__simple']; sign=-1 if mode.startswith('upgoing') else 1
    ax.imshow(section,extent=[mlags[0],mlags[-1],mdist[-1],mdist[0]],aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
    ax.plot(sign*mdist/3200,mdist,'k--',lw=1); ax.set(title=mode.replace('_',' '),xlabel='Lag (s)')
    axes[1,col].plot(trial/1000,matched_curve(section,sign),color='k',label='selected lag side')
    axes[1,col].plot(trial/1000,matched_curve(section,-sign),color='.55',ls='--',label='opposite lag side')
    axes[1,col].axvline(3.2,color='tab:red',ls=':'); axes[1,col].set(xlabel='Trial velocity (km/s)',ylabel='Median correlation')
axes[0,0].set_ylabel('Receiver offset (m)'); axes[1,0].legend(frameon=False,fontsize=8); plt.show()
for mode in modes:
    section=matched[f'{mode}__simple']; sign=-1 if mode.startswith('upgoing') else 1
    selected=matched_curve(section,sign); peak=int(np.argmax(selected)); at_3200=int(np.argmin(abs(trial-3200)))
    print(mode, 'selected-lag peak=',round(trial[peak]/1000,3),'km/s',
          'peak score=',round(selected[peak],4),'score@3.2=',round(selected[at_3200],4))
""")

markdown("""**Figure 6. Five-hour duration-matched comparison of unfiltered and F–K-filtered Lellouch-style correlations.** All five columns use the same 300 one-minute records from 20 December 2024 at full temporal and spatial sampling and identical strain-rate conversion, 5-s temporal normalization, 30-s windows with 15-s overlap, channel-0 virtual source, approximately 50-m receivers, local 21-channel receiver stack, and final 5–20 Hz bandpass. Only the pre-correlation F–K mask differs: none, the signed 2.5–4.5 km/s fan in each direction, or signed direction-only selection. The upper row shows correlation gathers with a common color scale and dashed 3.2-km/s trajectories on the lag side appropriate to each branch. The lower row shows median moveout score versus trial velocity for the selected lag side (solid black) and complementary lag side (gray dashed); red dotted lines mark 3.2 km/s. The two opposite signed velocity fans give nearly identical selected-lag curves: both peak at 2.675 km/s with scores of 0.504 and 0.503, respectively, whereas the unfiltered 3.2-km/s score is only 0.015. Opposite directions producing the same fan-bounded peak is evidence of mask conditioning, not evidence for two independently recovered arrivals. This is the clean test showing that the F–K operation, rather than a change in sampling, preprocessing, or geometry, is responsible for the visual and quantitative enhancement.""")

markdown("""### 4a. Decimation and anti-alias sensitivity

F–K resolution and aliasing depend on array spacing and aperture ([Rost and Thomas, 2002](https://doi.org/10.1029/2000RG000100)). The production code therefore has a separate five-hour sensitivity test: direct 2× slicing is compared with polyphase anti-aliased resampling in both time and channel coordinate and with processing at full temporal and spatial sampling. The F–K mask and score are otherwise frozen. This tests whether decimation creates the **conditional fan-filtered observable**; it is not an independent signal-detection test.
""")
python("""alias=np.load(products['alias']); alias_meta=json.loads(products['alias_meta'].read_text())
paths=[('direct','Direct 2× slicing'),('antialiased','Polyphase anti-aliased'),('full_resolution','Full resolution')]
fig,axes=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
for col,(key,label) in enumerate(paths):
    section=alias[f'{key}_top']; alags=alias[f'{key}_lags']; adist=alias[f'{key}_distances']; av=alias[f'{key}_velocities_m_s']; scores=alias[f'{key}_scores']; anull=alias[f'{key}_null']
    limit=np.nanpercentile(np.abs(section),98.5)
    axes[0,col].imshow(section,extent=[alags[0],alags[-1],adist[-1],adist[0]],aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
    axes[0,col].plot(adist/3200,adist,'k--',lw=1); axes[0,col].set(title=label,xlabel='Lag (s)',ylabel='Receiver offset (m)')
    axes[1,col].plot(av/1000,scores,color='tab:blue'); axes[1,col].axhline(np.quantile(anull,.95),color='tab:red',ls='--'); axes[1,col].axvline(3.2,color='k',ls=':')
    axes[1,col].set(xlabel='Trial velocity (km/s)',ylabel='Moveout score')
plt.show()
for key,label in paths: print(label,alias_meta['results'][key])
print('Decision: decimation implementation does not materially move the selected fan ridge; physical recovery remains unproven because Section 6 fails.')
""")

markdown("""**Figure 7. Five-hour anti-alias and decimation sensitivity of the selected fan-filtered observable.** Columns compare direct factor-of-two slicing, explicit polyphase anti-aliased resampling in time and channel coordinate, and full-resolution processing of the same 300 one-minute records. The upper panels show channel-0 correlation gathers and the lower panels show trial-velocity scores; red dashed lines are the path-specific 95th percentiles of receiver-order permutation scores and black dotted lines mark 3.2 km/s. Peak velocities are 3.075, 3.075, and 3.050 km/s, respectively, and peak scores differ by less than 0.005. Thus the conditional ridge is insensitive to the tested resampling implementation. Because all three paths impose the same velocity fan, this figure addresses aliasing stability only and is not evidence that the ridge exists before selection.""")

markdown("""## 5. White-noise input null

Ten independent white-noise ensembles, each comprising three 12-s arrays, were passed through detrending, 5–20 Hz filtering, 5-s running-absolute-mean normalization, 2× decimation, each frozen velocity/direction mask, channel-0 correlations, and the same moveout score. The synthetic and measured outputs are plotted below in the **same gather coordinates and normalization**. The synthetic duration is much shorter than the real stack and therefore does not provide a duration-matched p-value. This is a project-specific surrogate-data falsification test: surrogate methods test whether a statistic exceeds what the same analysis produces after destroying the property of interest ([Theiler et al., 1992](https://doi.org/10.1016/0167-2789(92)90102-S)). Passing white noise is necessary but not sufficient because white noise does not preserve the measured data's spectra, coherent artifacts, or nonstationarity.

The first realization is audited **before any filter**. Its samples must be independent in channel and time, its power must span the entire 0–250 Hz Nyquist band rather than only the target band, and cross-correlations between distinct channels must fluctuate around zero at the finite-record scale.
""")
python("""from ambient_transfer_test import preprocess
from ambient_fk_transfer_test import fk_filter as production_fk_filter
white=np.load(products['white_noise']); white_meta=json.loads(products['white_noise_meta'].read_text())
white_fs=float(white_meta['sampling']['fs_hz']); white_dx=float(white_meta['sampling']['dx_m'])
white_rng=np.random.default_rng(int(white_meta['seed']))
white_raw=white_rng.standard_normal((int(white_meta['sampling']['channels']),int(white_fs*white_meta['seconds_per_record'])),dtype=np.float32)
from scipy.signal import welch, correlate, correlation_lags
validation_channels=[0,49,196,392,784]
psd_frequency,psd_each=welch(white_raw[validation_channels],fs=white_fs,nperseg=1024,axis=1)
psd_relative=psd_each.mean(axis=0)/np.median(psd_each.mean(axis=0))
pair_a,pair_b=0,784
a=white_raw[pair_a]-white_raw[pair_a].mean(); b=white_raw[pair_b]-white_raw[pair_b].mean()
cross=correlate(a,b,mode='full',method='fft')/np.sqrt(np.sum(a*a)*np.sum(b*b))
cross_lag_samples=correlation_lags(a.size,b.size,mode='full'); cross_lags=cross_lag_samples/white_fs
lag_keep=np.abs(cross_lags)<=1.0; overlap=np.maximum(1,a.size-np.abs(cross_lag_samples)); cross_95=1.96/np.sqrt(overlap)
zero_lag_against_ch0=np.corrcoef(white_raw)[0,1:]
fig,axes=plt.subplots(1,3,figsize=(16,4.2),constrained_layout=True)
axes[0].plot(np.arange(1000)/white_fs,white_raw[pair_a,:1000],lw=.7,label=f'ch {pair_a}')
axes[0].plot(np.arange(1000)/white_fs,white_raw[pair_b,:1000],lw=.7,alpha=.75,label=f'ch {pair_b}')
axes[0].set(xlabel='Time (s)',ylabel='Generated amplitude',title='Independent input samples'); axes[0].legend(frameon=False)
axes[1].plot(psd_frequency,psd_relative,color='k',lw=1); axes[1].axvspan(5,20,color='tab:green',alpha=.15,label='later analysis band')
axes[1].axhline(1,color='.5',ls=':'); axes[1].set(xlim=(0,white_fs/2),xlabel='Frequency (Hz)',ylabel='PSD / median PSD',title='Broadband pre-filter spectrum'); axes[1].legend(frameon=False)
axes[2].plot(cross_lags[lag_keep],cross[lag_keep],color='k',lw=.8)
axes[2].fill_between(cross_lags[lag_keep],-cross_95[lag_keep],cross_95[lag_keep],color='tab:blue',alpha=.15,label='pointwise Gaussian 95% scale')
axes[2].axvline(0,color='.5',ls=':'); axes[2].set(xlabel='Lag (s)',ylabel='Normalized cross-correlation',title=f'Independent channels {pair_a} and {pair_b}'); axes[2].legend(frameon=False)
plt.show()
zero_index=np.flatnonzero(cross_lag_samples==0)[0]; outside_fraction=float(np.mean(np.abs(cross[lag_keep])>cross_95[lag_keep]))
print('Generator: numpy.random.default_rng(seed).standard_normal; seed =',white_meta['seed'])
print('Shape:',white_raw.shape,'sample rate:',white_fs,'Hz; Nyquist:',white_fs/2,'Hz; duration:',white_raw.shape[1]/white_fs,'s')
print('All-sample mean/std:',float(white_raw.mean()),float(white_raw.std()))
print('Mean normalized PSD below/inside/above 5–20 Hz:',float(psd_relative[psd_frequency<5].mean()),float(psd_relative[(psd_frequency>=5)&(psd_frequency<=20)].mean()),float(psd_relative[psd_frequency>20].mean()))
print('ch0–ch784 zero-lag correlation:',float(cross[zero_index]))
print('ch0 versus all 799 other channels: mean/std/max |zero-lag r|:',float(zero_lag_against_ch0.mean()),float(zero_lag_against_ch0.std()),float(np.max(np.abs(zero_lag_against_ch0))))
print('Fraction of |lag|≤1 s points outside pointwise 95% scale:',outside_fraction,'(about 0.05 is expected)')
""")

markdown("""**Figure 8a. Direct validation of the generated white-noise input before any filter.** The left panel overlays the first 2 s from channel 0 and channel 784, the farthest target receiver in the Lellouch-style geometry. The center panel is the mean Welch power spectral density of five widely separated channels, normalized by its median. Power spans the complete 0–250 Hz Nyquist band before preprocessing; the shaded 5–20 Hz interval only marks the band applied later. The right panel cross-correlates channels 0 and 784 over ±1 s and compares the result with the lag-dependent pointwise 95% scale expected for independent Gaussian samples with finite overlap. Printed diagnostics give the exact generator and seed, dimensions, mean and standard deviation, spectral power below/within/above 5–20 Hz, the selected pair's zero-lag correlation, channel 0's zero-lag correlations with all 799 other channels, and the fraction of lag samples outside the pointwise envelope. This validates that the input is broadband and mutually uncorrelated before the processing operator is applied.""")

python("""
white_processed=preprocess(white_raw,white_fs,norm_seconds=5.0)
white_filtered,white_filtered_fs,white_filtered_dx=production_fk_filter(white_processed,white_fs,white_dx,'negative')
white_dec=white_processed[::2,::2]; wf=np.fft.fftfreq(white_dec.shape[1],1/white_filtered_fs); wk=np.fft.fftfreq(white_dec.shape[0],white_filtered_dx)
wpower=np.abs(np.fft.fft2(white_dec))**2; worder=np.argsort(wk)
fig,axes=plt.subplots(1,4,figsize=(18,4.2),constrained_layout=True)
for ax,data,title,fs_plot,dx_plot in [
    (axes[0],white_raw,'Gaussian white-noise input',white_fs,white_dx),
    (axes[1],white_processed,'5–20 Hz + 5-s RAM',white_fs,white_dx),
    (axes[3],white_filtered,'F×K<0, 2.5–4.5 km/s output',white_filtered_fs,white_filtered_dx)]:
    limit=np.nanpercentile(np.abs(data),99)
    ax.imshow(data,extent=[0,data.shape[1]/fs_plot,(data.shape[0]-1)*dx_plot,0],aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit,interpolation='nearest')
    ax.set(title=title,xlabel='Time (s)',ylabel='Fiber coordinate (m)')
axes[2].pcolormesh(wf,wk[worder]*1000,10*np.log10(wpower[worder]/wpower.max()+1e-12),shading='auto',cmap='magma',vmin=-40,vmax=0)
ff=np.linspace(5,20,100)
for velocity,style in [(2500,'--'),(4500,':')]: axes[2].plot(ff,-1000*ff/velocity,color='cyan',ls=style,lw=1.2)
axes[2].set(xlim=(0,30),ylim=(-10,10),title='Pre-filter F–K power',xlabel='Frequency (Hz)',ylabel='Wavenumber (cycles/km)')
plt.show()
""")

markdown("""**Figure 8b. Validated white noise before and after the production F–K operator.** The broadband, independent 12-s, 800-channel Gaussian realization validated in Figure 8a is shown before preprocessing, after the same 5–20 Hz bandpass and 5-s running-absolute-mean normalization used by this null workflow, in the F–K domain, and after the F×K<0, 2.5–4.5 km/s mask. Cyan curves delimit the imposed velocity fan in the positive-frequency half-plane. Panel-specific amplitude limits are used because raw Gaussian samples, normalized samples, spectral power, and filtered amplitudes have different units. The final panel visibly contains sloping coherent texture even though the input contains no physical wave. That texture is the impulse response of the velocity-selective operator and is exactly why a filtered ridge cannot validate itself.""")

python("""def normalize_gather(section):
    scale=np.max(np.abs(section),axis=1,keepdims=True)
    return section/np.maximum(scale,np.finfo(float).eps)
w_lags=white['lags']; w_dist=white['distance']
fig,axes=plt.subplots(4,4,figsize=(16,13),sharex=True,sharey=True,constrained_layout=True)
columns=[('negative','Measured held-out'),('negative','White noise'),('positive','Measured held-out'),('positive','White noise')]
for row,(mask,label) in enumerate(masks):
    for col,(branch,source) in enumerate(columns):
        section=(fk[f'held_out__{mask}__{branch}__top'] if source.startswith('Measured') else white[f'{mask}__{branch}__example_top'])
        section=normalize_gather(section); sign=1 if branch=='negative' else -1
        axes[row,col].imshow(section,extent=[w_lags[0],w_lags[-1],w_dist[-1],w_dist[0]],aspect='auto',cmap='RdBu_r',vmin=-1,vmax=1,interpolation='nearest')
        axes[row,col].plot(sign*w_dist/3200,w_dist,'k--',lw=1)
        axes[row,col].set_title(f'{label}; {source}; {branch}')
        axes[row,col].set(xlabel='Correlation lag (s)',ylabel='Receiver offset (m)')
plt.show()
""")

markdown("""**Figure 8c. Direct gather-by-gather comparison of measured data and white noise after identical signed F–K masks.** Rows apply the production 2.5–4.5 km/s fan, narrow 2.8–3.8 km/s fan, broad 2.0–5.5 km/s fan, and direction-only mask. Columns alternate the held-out measured aggregate and the first predeclared white-noise ensemble for the F×K<0 and F×K>0 branches. Every panel uses the same lag and receiver-offset axes, each trace is divided by its own maximum absolute amplitude, and all panels use the common range −1 to 1. Dashed curves mark 3.2 km/s on the branch-appropriate lag side. The velocity-bounded white-noise panels develop moveout-like ridges comparable in geometry to the measured panels, especially for the narrow fan, whereas direction-only white noise does not. This is direct visual evidence that the fan can manufacture the expected geometry from noise; differences in absolute amplitude are intentionally not assessed because the measured and synthetic stack durations differ.""")

python("""fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True,constrained_layout=True)
for ax,(mask,label) in zip(axes.flat,masks):
    for branch,color in [('negative','tab:blue'),('positive','tab:orange')]:
        ensemble=white[f'{mask}__{branch}__ensemble_scores']; real_curve=fk[f'held_out__{mask}__{branch}__scores']
        lo,median,hi=np.quantile(ensemble,[.05,.5,.95],axis=0)
        ax.fill_between(white['velocities_m_s']/1000,lo,hi,color=color,alpha=.16)
        ax.plot(white['velocities_m_s']/1000,median,color=color,ls='--',lw=1,label=f'{branch} white median')
        ax.plot(white['velocities_m_s']/1000,real_curve,color=color,lw=2,label=f'{branch} measured')
    ax.axvline(3.2,color='k',ls=':'); ax.axhline(0,color='.6',lw=.8)
    ax.set(title=label,xlabel='Trial apparent velocity (km/s)',ylabel='Signed moveout score'); ax.legend(frameon=False,fontsize=8)
plt.show()
wi=np.argmin(abs(white['velocities_m_s']-3200))
for mask,label in masks:
    for branch in ['negative','positive']:
        result=white_meta['masks'][mask][branch]
        print(label,branch,'real exceeds white95:',result['real_exceeds_white_95_at_3200'])
""")

markdown("""**Figure 8d. Measured and white-noise moveout curves for the same masks.** Solid curves are the held-out measured-data scores; dashed curves and shaded intervals are the median and 5th–95th percentile across ten identically processed white-noise ensembles. Blue and orange denote the two signed branches, and the vertical dotted line marks 3.2 km/s. The narrow fan produces large scores over its retained band even for white noise—the median absolute score at 3.2 km/s is approximately 0.64—so the magnitude and location of that peak are mask-conditioned. The production-fan measured curve exceeds this short-duration white-noise ensemble at 3.2 km/s, but the measured-data scrambling test below is the stronger control because it preserves the real spectra and nonstationarity.""")

markdown("""## 6. Channel-scrambling test — input and resulting gathers

These are stronger than shuffling the finished correlation panel. For every real input record, channel 0 is fixed and all other preprocessed traces are reassigned to random fiber coordinates before F–K filtering. A second null independently circularly shifts each non-source trace before filtering. Both preserve important attributes of the measured traces while destroying ordered interchannel propagation. The nulls are tailored to this analysis; [Theiler et al. (1992)](https://doi.org/10.1016/0167-2789(92)90102-S) provides the general surrogate-test logic, not this exact geophysical implementation.

**Scope of this null product.** It uses the same 5–20 Hz preprocessing, 5-s normalization, 2× decimation, signed fan, channel-0 virtual source, 50-m receiver spacing, and moveout statistic as the production F–K analysis. It correlates each one-minute record directly and does not reproduce the Lellouch 30-s/15-s-overlap and nearby-receiver-stack details used in Section 4. Because this product fails rather than passes, it remains a valid warning that the fan can generate the selected statistic after spatial order is destroyed; it is not presented as an exact matched-pipeline p-value.
""")
python("""from ambient_fk_full_pipeline_null_v2 import stable_rng, channel_permutation
null=np.load(products['fk_nulls']); null_example=np.load(products['fk_null_example']); null_example_meta=json.loads(products['fk_null_example_meta'].read_text())
example_path=Path(str(null_example['used_files'][0]))
with h5py.File(example_path,'r') as handle:
    raw_group=handle['Acquisition/Raw[0]']; dataset=raw_group['RawData']; example_fs=float(raw_group.attrs.get('OutputDataRate',dataset.attrs.get('OutputDataRate',500.0))); example_dx=float(handle['Acquisition'].attrs.get('SpatialSamplingInterval',1.0))
    example_raw=np.asarray(dataset[:min(int(12*example_fs),dataset.shape[0]),:800],dtype='float32').T
ordered_input=preprocess(example_raw,example_fs,norm_seconds=5.0)
scramble_rng=stable_rng(int(null_example_meta['random_seed']),'channel_permutation',0,str(example_path))
scrambled_input=channel_permutation(ordered_input,scramble_rng)
ordered_filtered,ordered_fs,ordered_dx=production_fk_filter(ordered_input,example_fs,example_dx,'negative')
scrambled_filtered,_,_=production_fk_filter(scrambled_input,example_fs,example_dx,'negative')
fig,axes=plt.subplots(1,4,figsize=(18,4.2),constrained_layout=True)
input_limit=np.nanpercentile(np.abs(np.concatenate([ordered_input.ravel(),scrambled_input.ravel()])),99)
filtered_limit=np.nanpercentile(np.abs(np.concatenate([ordered_filtered.ravel(),scrambled_filtered.ravel()])),99)
for ax,data,title,fs_plot,dx_plot,limit in [
    (axes[0],ordered_input,'Measured input: ordered channels',example_fs,example_dx,input_limit),
    (axes[1],scrambled_input,'Same input: channels scrambled',example_fs,example_dx,input_limit),
    (axes[2],ordered_filtered,'Ordered → production fan',ordered_fs,ordered_dx,filtered_limit),
    (axes[3],scrambled_filtered,'Scrambled → production fan',ordered_fs,ordered_dx,filtered_limit)]:
    ax.imshow(data,extent=[0,data.shape[1]/fs_plot,(data.shape[0]-1)*dx_plot,0],aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit,interpolation='nearest')
    ax.set(title=title,xlabel='Time (s)',ylabel='Assigned fiber coordinate (m)')
plt.show()
""")

markdown("""**Figure 9a. What the channel-scrambling null does before correlation.** A 12-s excerpt from the first file in null realization 0 is shown after the production 5–20 Hz and 5-s RAM preprocessing with its measured channel order, after all non-source channels are assigned to random fiber coordinates while channel 0 remains fixed, and after the same F×K<0, 2.5–4.5 km/s filter is applied to each input. The random generator is fixed by seed 20260805, method name, realization ID 0, and file path; this example was not selected by appearance. Ordered and scrambled inputs share one amplitude scale, and their filtered outputs share another. The fan imposes sloping coherent texture even after the physical channel order has been destroyed, demonstrating the mechanism tested by the 300-file correlation gathers below.""")

python("""n_lags=null_example['lags_s']; n_dist=null_example['receiver_offsets_m']
fig,axes=plt.subplots(2,3,figsize=(14,8),sharex=True,sharey=True,constrained_layout=True)
for row,branch in enumerate(['negative','positive']):
    sign=1 if branch=='negative' else -1
    panels=[
        (null_example[f'observed_{branch}_top'],'Measured, ordered'),
        (null_example[f'null_channel_permutation_{branch}_top_first'],'Channel-scrambled, r0'),
        (null_example[f'null_circular_time_shift_{branch}_top_first'],'Time-shifted, r0')]
    for col,(section,title) in enumerate(panels):
        section=normalize_gather(section)
        axes[row,col].imshow(section,extent=[n_lags[0],n_lags[-1],n_dist[-1],n_dist[0]],aspect='auto',cmap='RdBu_r',vmin=-1,vmax=1,interpolation='nearest')
        axes[row,col].plot(sign*n_dist/3200,n_dist,'k--',lw=1)
        axes[row,col].set(title=f'{title}; {branch}',xlabel='Correlation lag (s)',ylabel='Receiver offset (m)')
plt.show()
""")

markdown("""**Figure 9b. Legacy single-pair F–K diagnostic; not a null test of the published stack.** Rows show the F×K<0 and F×K>0 branches; columns show the ordered 300-file measured stack, fixed channel-scramble realization 0, and fixed independent-time-shift realization 0. Every panel uses the same lag and receiver-offset axes, trace-wise normalization, color limits, and branch-appropriate 3.2-km/s trajectory. Both coherence-destroying surrogates retain a fan-aligned ridge that is at least as visually prominent as the ordered-data ridge. Within this legacy operator, the apparent moveout survives removal of spatial ordering. Because this operator omits the published 30-s/15-s nearby-pair and two-pass stack, the result cannot decide whether the Lellouch-style reproduction passes or fails.""")

python("""vv=null['velocities_m_s']; i=np.argmin(abs(vv-3200)); observed=abs(float(null['observed_negative_scores'][i]))
cp=np.abs(null['null_channel_permutation_negative_scores'][:,i]); ts=np.abs(null['null_circular_time_shift_negative_scores'][:,i])
fig,axes=plt.subplots(1,2,figsize=(13,4.8),sharey=True,constrained_layout=True)
for ax,branch in zip(axes,['negative','positive']):
    observed_curve=null[f'observed_{branch}_scores']; cp_curves=null[f'null_channel_permutation_{branch}_scores']; ts_curves=null[f'null_circular_time_shift_{branch}_scores']
    cp_lo,cp_hi=np.quantile(cp_curves,[.1,.9],axis=0); ts_lo,ts_hi=np.quantile(ts_curves,[.1,.9],axis=0)
    ax.fill_between(vv/1000,cp_lo,cp_hi,color='tab:blue',alpha=.15,label='scramble 10–90%')
    ax.fill_between(vv/1000,ts_lo,ts_hi,color='tab:orange',alpha=.15,label='time-shift 10–90%')
    ax.plot(vv/1000,observed_curve,color='k',lw=2,label='measured ordered')
    ax.plot(vv/1000,null_example[f'null_channel_permutation_{branch}_scores'][0],color='tab:blue',lw=1.2,label='scramble r0')
    ax.plot(vv/1000,null_example[f'null_circular_time_shift_{branch}_scores'][0],color='tab:orange',lw=1.2,label='time-shift r0')
    ax.axvspan(2.5,4.5,color='.8',alpha=.15); ax.axvline(3.2,color='k',ls=':'); ax.axhline(0,color='.6',lw=.8)
    ax.set(xlim=(2,5),title=branch,xlabel='Trial apparent velocity (km/s)',ylabel='Signed moveout score'); ax.legend(frameon=False,fontsize=8)
plt.show()
print('Real score:',observed,'channel-scramble 95%:',np.quantile(cp,.95),'time-shift 95%:',np.quantile(ts,.95))
print('Decision for legacy single-pair operator only: FAIL. This is not a matched test of the published nearby-pair stack.')
""")

markdown("""**Figure 9c. Legacy single-pair velocity-score nulls.** Black curves are the ordered 300-file scores. Blue and orange curves show the fixed realization-0 channel-scramble and time-shift scores; shaded intervals span the 10th–90th percentiles of all 20 realizations. The gray band is the imposed 2.5–4.5 km/s velocity fan and the dotted line marks 3.2 km/s. At 3.2 km/s, the ordered F×K<0 absolute score is 0.266, compared with 95th-percentile values of 0.359 for channel scrambling and 0.326 for time shifting; every stored surrogate realization equals or exceeds the ordered score. Thus the legacy F–K/single-pair operator fails its own spatial-order and phase-coherence controls. It does not test the complete published nearby-pair stack and must not be used to reject that workflow.""")

markdown("""## 7. Is F–K required, and is the recovered feature defensible?

The phrase “F–K is required” has two meanings. The visual statement is true: a velocity-bounded fan makes an ordered ridge appear while the unfiltered and direction-only paths remain weak near 3.2 km/s. The scientific statement remains unresolved. The five-hour F–K comparison uses the nearby-receiver display stack, but the available measured-data surrogate test comes from the older single-pair pipeline. Those products cannot be combined into one acceptance decision. A decisive test must send ordered data and every surrogate through the identical 30-s/15-s, R±10, first-pass velocity, and shifted-restack sequence. This conservative decision is consistent with Green's-function theory: an attractive cross-correlation is not automatically a Green's function when illumination and source-boundary assumptions are unverified ([Wapenaar and Fokkema, 2006](https://doi.org/10.1190/1.2213955)).
""")
python("""un_meta=json.loads(products['unfiltered_meta'].read_text()); directional=json.loads(products['directional_meta'].read_text())
rows=[
 ('Unfiltered eight-day stack',False,un_meta['score_3200'],un_meta['p_peak'],'No independent 3.2 km/s recovery'),
 ('Direction-only F–K, downgoing',False,directional['key_values']['held_out_direction_only_negative_score_3200'],None,'No recovery without velocity restriction'),
 ('2.5–4.5 km/s fan, downgoing',True,observed,1.0,'Visible ridge, but fails pre-filter null'),
]
for row in rows: print(row)
assert observed < np.quantile(cp,.95) and observed < np.quantile(ts,.95)
""")

markdown("""## 8. Stack-duration convergence: 5 hours through 8 days

Every combination of the eight selected complete days is compared with the final eight-day section; the independently processed five-hour checkpoint is plotted separately. Eight days, rather than the suggested two-week upper bound, is the maximum supported by the frozen complete-day product set used here; no 14-day result is implied. Shading is the 10th–90th percentile across combinations. Comparing short stacks with a longer reference follows the convergence logic tested by [Seats et al. (2012)](https://doi.org/10.1111/j.1365-246X.2011.05263.x), while [Bensen et al. (2007)](https://doi.org/10.1111/j.1365-246X.2007.03374.x) emphasize temporal repeatability as a reliability criterion. Because each day subset is contained in the eight-day reference, the subset-to-reference correlations are positively dependent and the threshold is descriptive rather than an independent error estimate. The selected fan is internally stable after one day; the unfiltered section reaches the same descriptive threshold only at seven days. Because the fan fails the input-level null, its rapid “convergence” is interpreted as operator stability, not physical Green’s-function convergence.
""")
python("""conv=np.load(products['convergence']); conv_meta=json.loads(products['convergence_meta'].read_text()); days=np.arange(1,9)
fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
for key,color,label in [('unfiltered','0.25','unfiltered'),('filtered','tab:blue','2.5–4.5 km/s fan')]:
    corr=conv[f'{key}_corr_quantiles']; scoreq=conv[f'{key}_score_quantiles']
    axes[0].plot(days,corr[:,1],'o-',color=color,label=label); axes[0].fill_between(days,corr[:,0],corr[:,2],color=color,alpha=.18)
    axes[1].plot(days,scoreq[:,1],'o-',color=color,label=label); axes[1].fill_between(days,scoreq[:,0],scoreq[:,2],color=color,alpha=.18)
axes[0].axhline(.9,color='tab:red',ls='--'); axes[0].set(xlabel='Stacked complete days',ylabel='Correlation with eight-day section',title='Section stability')
axes[1].set(xlabel='Stacked complete days',ylabel='3.2 km/s score',title='Moveout-score stability')
five=conv_meta['five_hour']
axes[0].scatter([5/24],[five['unfiltered_section_correlation']],marker='s',color='.25',zorder=4)
axes[0].scatter([5/24],[five['filtered_section_correlation']],marker='s',color='tab:blue',zorder=4)
axes[1].scatter([5/24],[five['unfiltered_score_3200']],marker='s',color='.25',zorder=4)
axes[1].scatter([5/24],[five['filtered_score_3200']],marker='s',color='tab:blue',zorder=4)
for ax in axes: ax.set_xlim(0,8.25)
for ax in axes: ax.legend(frameon=False)
plt.show(); print(conv_meta)
""")

markdown("""**Figure 10. Stack-duration stability relative to the eight-day reference.** For each stack length from one to eight complete selected days, every available day combination is averaged and compared with the final eight-day section. Square symbols show the separately processed five-hour checkpoint. Lines show medians and shaded intervals show the 10th–90th percentiles across day combinations. The left panel measures full-section waveform correlation with the eight-day reference; the red dashed line marks the predeclared descriptive threshold of 0.9. The right panel shows the 3.2-km/s moveout score. Because each partial day stack contributes to the eight-day reference, these correlations are positively dependent and cannot be read as independent prediction accuracy. The fan-filtered section crosses the internal threshold after approximately one day, whereas the unfiltered section does so only at seven days. Because the fan fails the input-level surrogate tests in Figure 9, its fast stability is interpreted as convergence of the selected operator output, not convergence of a physical Green's function. A defensible physical-signal minimum is not established between five hours and the tested eight days.""")

markdown("""## 9. Decision and scope

The corrected one-day Figure 7c test is complete. All required gates passed operationally: 1,440 files; 5,759 unique windows; the fixed-top unshifted R±10 gather; no F–K filter; a complete velocity scan; receiver-order, white-noise, and pre-correlation channel-scramble controls; source-coordinate and RAM sensitivities; and separately labelled common-mode and Equation-6 sensitivities.

The scientific outcome is negative for **20 December 2024**:

- the measured scan maximum occurs at 5.85 km s⁻¹, near the flat-moveout edge of the 1.5–6.0 km s⁻¹ search, rather than at the published approximately 3.2 km s⁻¹;
- the measured scan maximum is below the receiver-order 95th percentile (6.13 versus 6.32; familywise p = 0.147);
- the causal and acausal scores at 3.2 km s⁻¹ are comparable (2.75 and 2.83), not a clear causal-dominant packet;
- the baseline and pre-correlation channel-scramble gathers have flattened correlation 0.9976, showing that the dominant waveform is insensitive to receiver order; and
- neither 900-channel common-mode subtraction nor the declared post-average Equation-6 stabilization reveals a significant 3.2 km s⁻¹ mode.

Therefore the exact published one-day operator does **not** reproduce Figure 7c on this matched 2024 day. This conclusion is narrower than “ambient interferometry fails” and narrower than “F–K filtering is invalid.” The next clean scientific question is whether the same exact operator converges on independently chosen days or a predeclared multi-day stack; any F–K-assisted result remains a separately labelled extension whose evidence must include matched input-level controls.
""")

if EXACT_SUMMARY.is_file():
    markdown("""## 10. Paper-faithful Figure 7c result — v8

The dependency-gated full-day calculation and v4 aggregation completed. The table and figures below are generated from the corrected operator and matched controls, not from the withdrawn legacy checkpoints. Hourly baseline spectra retain v2 provenance because the v3 sensitivity and v4 lag-filter corrections leave their cross spectra unchanged; v4 reconstructs, filters, tests, and plots every aggregate.""")
    python("""import json
import numpy as np
import pandas as pd
from IPython.display import display, Image

exact_dir = ROOT/'ambient_transfer'/'lellouch2019_exact_stack'
branch_stems = {
    'paper baseline': 'aggregate_2024-12-20_src23_ram0p1_cross_correlation_ordered_r0',
    'source channel 0': 'aggregate_2024-12-20_src0_ram0p1_cross_correlation_ordered_r0',
    'RAM 5 s': 'aggregate_2024-12-20_src23_ram5_cross_correlation_ordered_r0',
    'common mode (900 ch)': 'aggregate_2024-12-20_src23_ram0p1_cross_correlation_ordered_r0_cm',
    'Equation 6 stabilized': 'aggregate_2024-12-20_src23_ram0p1_source_power_stabilized_ordered_r0',
    'white noise': 'aggregate_2024-12-20_src23_ram0p1_cross_correlation_white_noise_r0',
    'channel permutation': 'aggregate_2024-12-20_src23_ram0p1_cross_correlation_channel_permutation_r0',
}
branch_json = {}
branch_npz = {}
for label, stem in branch_stems.items():
    json_path = exact_dir/(stem+'.json')
    npz_path = exact_dir/(stem+'.npz')
    assert json_path.is_file() and npz_path.is_file(), (json_path, npz_path)
    item = json.loads(json_path.read_text())
    data = np.load(npz_path, allow_pickle=False)
    assert item['workflow_version'] == 'lellouch2019_exact_stack_v4'
    assert item['files'] == 1440 and item['windows_30_s_15_s_step'] == 5759
    assert item['fk_filter'] is False
    assert item['bandpass_order_relative_to_lag_crop'] == 'full correlation first; crop second'
    branch_json[label], branch_npz[label] = item, data

baseline = branch_npz['paper baseline']['r_plus_minus_10_correlation']
def gather_similarity(test):
    row_correlation = [np.corrcoef(reference, candidate)[0, 1]
                       for reference, candidate in zip(baseline, test)]
    return float(np.median(row_correlation)), float(np.corrcoef(
        baseline.ravel(), test.ravel())[0, 1])

records = []
for label in branch_stems:
    item, data = branch_json[label], branch_npz[label]
    median_row, flattened = gather_similarity(data['r_plus_minus_10_correlation'])
    records.append({
        'branch': label,
        'best v (m/s)': item['best_causal_velocity_m_s'],
        'best score': item['best_causal_score'],
        'score at 3200': item['causal_score_at_3200'],
        'acausal at 3200': item['acausal_score_at_3200'],
        'null 95% scan max': float(np.percentile(
            data['receiver_order_null_maxima'], 95)),
        'familywise p': item['receiver_order_scan_max_null_p'],
        'median trace corr. to baseline': median_row,
        'flattened corr. to baseline': flattened,
    })
display(pd.DataFrame(records).round(4))

for label in ('paper baseline', 'common mode (900 ch)',
              'white noise', 'channel permutation'):
    display(Image(filename=str(exact_dir/(branch_stems[label]+'.png')), width=1100))
""")
    markdown("""**Figure 11. Paper-faithful full-day Figure 7c test and matched falsification controls.** The first three-panel figure is the measured-data baseline: (a) individual center-receiver correlations, (b) the literal unshifted R−10:R+10 sum reported for Figure 7c, and (c) causal and acausal apparent-velocity scans. The dotted horizontal line in panel (c) is the 95th percentile of 10,000 receiver-order null scan maxima; the vertical 3.2 km s⁻¹ line and red moveout trajectory are comparison references, not fitted or imposed filters. The baseline maximum occurs at 5.85 km s⁻¹ with score 6.13, below the null threshold 6.32 (familywise p = 0.147), while the causal and acausal scores at 3.2 km s⁻¹ are 2.75 and 2.83. The second figure subtracts the instantaneous median estimated from all 900 acquisition channels before RAM; it removes the dominant shared waveform but does not expose significant ordered moveout (p = 0.922). The third figure passes deterministic-seed iid broadband Gaussian input through the complete operator and remains nonsignificant (p = 0.337). The fourth permutes measured receiver channels before differentiation and RAM, preserving temporal content while destroying spatial order; it is also nonsignificant (p = 0.901) yet remains nearly identical to the measured baseline (flattened gather correlation 0.9976). Every branch uses 1,440 files, 5,759 windows, full-correlation 5–20 Hz filtering before lag cropping, and no F–K filter. Together these results show that the exact one-day 2024 calculation does not recover the ordered approximately 3.2 km s⁻¹ packet shown by Lellouch et al. (2019). They do not test whether a separately controlled F–K extension or a longer, independently selected stack can recover a physical signal.""")
else:
    markdown("""## 10. Paper-faithful Figure 7c result — pending

The authoritative implementation is `ambient_lellouch2019_exact_stack.py` at commit `bb48942`. Full-day matrix job `38988141` is queued/running on Sherlock. Dependency-gated aggregation job `38988173` will continue automatically only after all 168 hourly branch tasks succeed. No user notification is required. Until the aggregate JSON appears, all legacy Figure 7c and F–K verdicts remain withdrawn.""")
markdown("""## 11. Beyond one day — five more days, a coherent stack, and why a small p is still not a reproduction — v10 (census claims withdrawn)

Section 10 tests 20 December 2024 alone, matched to the paper's own one-day design.
Three questions remained: whether another day does better, whether stacking more
data recovers the arrival, and what to conclude if some subset does return a small
p value. All three are now answered.

**Every complete day, not just the tested one.** One day is a thin basis for a
statement about an archive, so configuration 0 was run on every other complete day.
(This was originally motivated by a raw energy census that appeared to rank the days;
that census is **withdrawn** — see the note at the end of this section — so the days
are effectively an unselected set, which is a stronger basis for the result than a
selected one.) Configuration 0 was therefore run on every other
complete day. Days whose manifest is not exactly continuous at 60.000 s cannot be
streamed end to end; 2024-11-30, 2024-10-28 and 2025-03-04 each carry two timing
anomalies, and 2024-11-30 was analysed over its
leading 21.3-hour contiguous block rather than excluded.

| date | peak | at (m/s) | causal/acausal at 3,200 | p |
|---|---:|---:|---:|---:|---:|
| 2024-11-30 (21.3 h) | 6.536 | 5850 | 0.96 | 0.1345 |
| 2024-12-20 | 6.131 | 5850 | 0.97 | 0.1470 |
| 2024-06-17 | 1.901 | 5925 | 1.12 | 0.7413 |
| 2024-06-26 | 1.297 | 5925 | 0.86 | 0.4206 |
| 2025-02-24 | 2.532 | 5875 | 1.05 | 0.7517 |
| 2024-05-11 | 1.026 | 1675 | 1.02 | 0.3091 |

No day reaches significance; the minimum over six days is 0.1345, on the richest
day. Fisher's combination across the independent days gives χ² = 9.08 on 10 degrees
of freedom, p = 0.524 — the χ² ≈ df expected from noise. Note also that the observed
score never clears the **per-velocity** null at any of 181 velocities, so the negative
does not depend on the familywise correction at all.

**Peak velocities in the table are grid ceilings, not velocity estimates.** The score
samples each trace in a gate at t = offset/v, so as v rises every gate slides toward
zero lag where these gathers carry a dominant broad lobe; hence
corr(trial velocity, score) = +0.976 and the argmax lands wherever the scan stops.
Capping the grid at 3,500 / 4,000 / 8,000 / 20,000 m/s moves the "peak" to
3,475 / 3,925 / 7,700 / 18,775 m/s. A moveout-free control — every trace replaced by
the across-receiver median — reproduces the observed curve to within 3.5 % at every
velocity, including the exact peak location. The p-values remain valid because the
receiver-order permutation preserves the same gate geometry (97 % of null curves also
peak at ≥ 5,500 m/s).

**Withdrawn (2026-08-14).** The raw energy census that originally motivated the day
selection, and that was offered as the mechanism for the non-reproduction, is
retracted. It had no geometric baseline: 98.4 % of in-band (f,k) cells lie below
1,500 m/s by construction, so its "81–99 % of energy below 1,500 m/s" is at or below
the white-noise expectation and cannot distinguish the SAFOD wavefield from noise; per
cell the data are in fact 2.6–2.8× enriched at body-wave velocities. Its "downgoing
share never exceeds 49.9 %" is false — the tree's own outputs include 51.6 %. A
cross-epoch replacement using Lellouch's released 2017 earthquake records is likewise
withdrawn: the two arms were not processed alike, because `extract_all.py` takes
DASutils `median=True` by default and stripped the common mode from the 2024–25 arm
only. **The non-reproduction therefore has no validated mechanism at present.** The
per-day and multi-day numbers above are unaffected and were verified by an independent
chunk re-sum reproducing the stored aggregates to max abs difference 0.000e+00.

**A coherent stack, not a pooled p value.** Combining p values is not the same
experiment as stacking data: a coherent arrival grows with stacked windows while
noise averages down. Because chunk files store summed cross spectra, days are
exactly additive. Stacking the four full days that share an acquisition rate —
5,760 files, 23,036 windows, 96.0 hours — gives peak 1.912 at 5,925 m/s against a
null 95th percentile of 2.649, p = 0.9184, and causal/acausal of 0.99 at 3,200 m/s.
Four times the data leaves the statistic further from threshold than the single best
day.

**The archive is not homogeneous in acquisition rate.** 2024-05-11 was recorded at
5000 Hz against 500 Hz on every other day, so its cross spectra lie on a different
frequency grid and cannot be coherently pooled. Any future analysis combining days
must group by rate.

**Why a p below 0.05 here would still not be a reproduction.** Stacking the two
best-scoring days returns p = 0.039, and the two days were selected precisely
because they scored lowest. Running all ten pairs of the five usable 500 Hz days
gives two below 0.05 against 0.5 expected by chance; binomial P(≥2) ≈ 0.086, and the
pairs share days so are not independent. The decisive objection, however, is not
multiplicity but velocity: **every pair reaching p < 0.05 peaks at 5,850–5,900 m/s
with causal/acausal of 0.97–1.01.** A maximum at the top edge of the declared scan
is flat moveout — energy arriving at all receivers within a few milliseconds across
700 m — which is the signature of common-mode or instrumental structure, not a
propagating body wave. Figure 7c is a specific claim: a packet near 3,200 m/s with
the causal side dominant. At 3,200 m/s the score is consistently about half the peak
and the causal side never dominates.

This archive therefore does contain a statistically detectable coherent component,
and it is not the arrival Lellouch et al. (2019) report. Accepting it on its p value
alone would repeat, in a new form, the error that Sections 5–7 document for the F–K
fan. The multi-day tool now requires three conditions jointly — p < 0.05, a peak
inside the 2,500–4,000 m/s fan, and causal dominance — before it will report a
recovery.
""")

python("""from pathlib import Path
import numpy as np
from IPython.display import Image, display

base = ROOT if 'ROOT' in globals() else Path.cwd()
days_dir = base/'ambient_transfer'/'lellouch2019_exact_stack_days'
auth = base/'ambient_transfer'/'lellouch2019_exact_stack'/'_authoritative_38993456'

rows = []
for f in sorted(days_dir.glob('aggregate_*_src23_ram0p1_cross_correlation_ordered_r0.npz')) + \
         sorted(auth.glob('aggregate_2024-12-20_src23_ram0p1_cross_correlation_ordered_r0.npz')):
    z = np.load(f, allow_pickle=True)
    vg = z['velocity_grid_m_s']; cs = z['causal_moveout_scores']; ac = z['acausal_moveout_scores']
    nl = z['receiver_order_null_maxima']; k = int(np.nanargmax(cs))
    rows.append((f.name.split('_')[1], int(z['n_windows']), float(cs[k]), float(vg[k]),
                 float(np.interp(3200, vg, cs)/np.interp(3200, vg, ac)),
                 float((np.sum(nl >= cs[k])+1)/(len(nl)+1))))
rows.sort(key=lambda r: r[5])
print(f"{'date':<12}{'windows':>9}{'peak':>9}{'at m/s':>9}{'c/a@3200':>10}{'p':>9}")
for r in rows:
    print(f"{r[0]:<12}{r[1]:>9}{r[2]:>9.3f}{r[3]:>9.0f}{r[4]:>10.2f}{r[5]:>9.4f}")
print()
print('No day clears alpha = 0.05. Every peak but one sits at the top edge of the scan,')
print('i.e. flat moveout, rather than at the ~3200 m/s the paper reports.')

for png in (base/'ambient_lellouch2019_multiday_stack.png',
            base/'ambient_lellouch2019_multiday_stack_best2.png'):
    if png.is_file():
        display(Image(filename=str(png), width=1100))
""")

markdown("""**Figure 12. Coherent multi-day stacks.** The first panel pair is the 96.0-hour stack of the four full days sharing a 500 Hz acquisition rate; the second is the two best-scoring days only. Left panels show the R−10:R+10 gather with the 3.2 km s⁻¹ trajectory overlaid as a reference, not a fit; right panels show the causal velocity scan against the 95th percentile of the receiver-order null. Neither stack places its maximum near 3.2 km s⁻¹. The 96-hour stack is nonsignificant (p = 0.918); the hand-picked two-day stack reaches p = 0.039 but peaks at 5,850 m s⁻¹ with causal/acausal 0.97, which is flat moveout rather than a propagating arrival, and was selected as the best of ten pairs. Neither is a reproduction of Figure 7c.""")

markdown("""## References

- Bensen, G. D., et al. (2007), Processing seismic ambient noise data to obtain reliable broad-band surface wave dispersion measurements, *Geophysical Journal International*, 169, 1239–1260. [https://doi.org/10.1111/j.1365-246X.2007.03374.x](https://doi.org/10.1111/j.1365-246X.2007.03374.x)
- Ehsaninezhad, L., C. Wollin, V. Rodríguez Tribaldos, B. Schwarz, and C. M. Krawczyk (2024), Urban subsurface exploration improved by denoising of virtual shot gathers from distributed acoustic sensing ambient noise, *Geophysical Journal International*, 237, 1751–1764. [https://doi.org/10.1093/gji/ggae134](https://doi.org/10.1093/gji/ggae134)
- Isken, M. P., H. Vasyura-Bathke, T. Dahm, and S. Heimann (2022), De-noising distributed acoustic sensing data using an adaptive frequency–wavenumber filter, *Geophysical Journal International*, 231, 944–949. [https://doi.org/10.1093/gji/ggac229](https://doi.org/10.1093/gji/ggac229)
- Lellouch, A., et al. (2019), Seismic velocity estimation using passive downhole distributed acoustic sensing records: Examples from the San Andreas Fault Observatory at Depth, *Journal of Geophysical Research: Solid Earth*, 124, 6931–6948. [https://doi.org/10.1029/2019JB017533](https://doi.org/10.1029/2019JB017533)
- Rost, S., and C. Thomas (2002), Array seismology: Methods and applications, *Reviews of Geophysics*, 40, 1008. [https://doi.org/10.1029/2000RG000100](https://doi.org/10.1029/2000RG000100)
- Seats, K. J., J. F. Lawrence, and G. A. Prieto (2012), Improved ambient noise correlation functions using Welch's method, *Geophysical Journal International*, 188, 513–523. [https://doi.org/10.1111/j.1365-246X.2011.05263.x](https://doi.org/10.1111/j.1365-246X.2011.05263.x)
- Shapiro, N. M., and M. Campillo (2004), Emergence of broadband Rayleigh waves from correlations of the ambient seismic noise, *Geophysical Research Letters*, 31, L07614. [https://doi.org/10.1029/2004GL019491](https://doi.org/10.1029/2004GL019491)
- Suprajitno, M., and S. A. Greenhalgh (1985), Separation of upgoing and downgoing waves in vertical seismic profiling by contour-slice filtering, *Geophysics*, 50, 950–962. [https://doi.org/10.1190/1.1441973](https://doi.org/10.1190/1.1441973)
- Theiler, J., et al. (1992), Testing for nonlinearity in time series: The method of surrogate data, *Physica D*, 58, 77–94. [https://doi.org/10.1016/0167-2789(92)90102-S](https://doi.org/10.1016/0167-2789(92)90102-S)
- Wapenaar, K., and J. Fokkema (2006), Green's function representations for seismic interferometry, *Geophysics*, 71, SI33–SI46. [https://doi.org/10.1190/1.2213955](https://doi.org/10.1190/1.2213955)
- Wapenaar, K. (2006), Green's function retrieval by cross-correlation in case of one-sided illumination, *Geophysical Research Letters*, 33, L19304. [https://doi.org/10.1029/2006GL027747](https://doi.org/10.1029/2006GL027747)
""")

nbf.write(nb, OUT)
print(OUT)
