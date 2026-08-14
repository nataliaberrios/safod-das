"""Build the audited advisor-facing ambient F-K QC notebook."""
from pathlib import Path
import nbformat as nbf

HERE = Path(__file__).resolve().parent
OUT = HERE / "Ambient_FK_QC_workflow.ipynb"
nb = nbf.v4.new_notebook()
nb.metadata = {"kernelspec": {"display_name": "Python [conda env:das]", "language": "python", "name": "das"}}

def markdown(text):
    nb.cells.append(nbf.v4.new_markdown_cell(text))

def python(text):
    nb.cells.append(nbf.v4.new_code_cell(text))

markdown("""# Ambient-noise F–K QC workflow — v2

This is the single advisor-facing notebook for the ambient-noise analysis. It follows one linear chain:

1. reproduce the Lellouch-style unfiltered calculation;
2. show the raw → normalized → correlation → F–K progression;
3. compare velocity bands and both signed coordinate-direction branches;
4. run white-noise and pre-filter channel-scramble nulls through the same operator;
5. show full-pipeline null statistics; and
6. place the available 5-hour, one-day, and multi-day products on a stack-duration axis.

The notebook is a decision document, not a gallery. A processing choice passes only if the real-data result survives the corresponding input-level null. The current 2.5–4.5 km/s fan produces an apparent ridge but fails the pre-filter channel-scramble gate, so it is not accepted as an independently recovered physical arrival.
""")

markdown("""v2""")

markdown("""## Frozen analysis contract

The Lellouch baseline uses channel 0 as the virtual source, receivers at approximately 50 m increments, 5–20 Hz filtering, 30 s windows with 15 s overlap, a time-derivative/strain-rate proxy, running absolute-mean temporal normalization, and 3.2 km/s alignment only where explicitly shown. Earthquake-containing records remain in the archive; temporal normalization suppresses their amplitude leverage but does not remove their waveforms.

Lellouch et al. (2019) used F–K filtering for earthquake S-wave isolation, not in the ambient-interferometry subsection. Our ambient F–K work is therefore an explicitly labelled extension. Channel coordinate increases downhole. A synthetic check maps F×K<0, evaluated at positive lag, to increasing-coordinate propagation and F×K>0, evaluated at negative lag, to decreasing-coordinate propagation. The known surface-AWD arrival independently validates that physical mapping at 25–60 Hz; its 5–20 Hz energy is nearly directionally balanced, so the ambient-band physical labels retain that caveat.
""")

markdown("""## Published basis and project-specific extensions

The workflow deliberately separates what is reported in the literature from what is being tested here.

| Choice | Basis | Status in this notebook |
|---|---|---|
| Differentiate DAS strain in time | [Lellouch et al. (2019)](https://doi.org/10.1029/2019JB017533) convert the 2017 SAFOD strain records to strain rate. | Reproduced. |
| Correlate short, overlapping windows and stack them | Lellouch et al. use 30-s windows with 15-s overlap at SAFOD; [Seats et al. (2012)](https://doi.org/10.1111/j.1365-246X.2011.05263.x) analyze why overlapping-window correlation can improve ambient-noise convergence. | Reproduced, then convergence-tested. |
| Reduce transient amplitude leverage before correlation | Lellouch et al. use running-absolute-mean normalization; [Bensen et al. (2007)](https://doi.org/10.1111/j.1365-246X.2007.03374.x) describe temporal normalization as a way to suppress earthquake and instrumental transients. | Reproduced in form. The **5-s normalization window is a declared project parameter**, not a value reported by Lellouch et al. |
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
| Physical array | Cemented fiber in the SAFOD main hole | Same physical fiber, assumed unchanged channel origin | Reasonable continuity assumption, documented rather than independently re-surveyed. |
| Recording epoch | June–July 2017 | 2024–2025 archive; matched example is 20 December 2024 | The ambient source field and instrument response can differ. |
| Interrogator sampling | 2,500 Hz, approximately 1-m channel spacing, 10-m gauge length | 500 Hz and 1.020952-m header spacing; current gauge/interrogator metadata must be taken from the 2024 acquisition records | The method is reproduced, not the acquisition hardware. |
| DAS observable | Strain differentiated in time to strain rate | Same | Direct reproduction. |
| Temporal normalization | Running absolute mean | Same form; 5-s window frozen here | The 5-s duration is a project choice because the paper does not report that window length. |
| Correlation segmentation | 30-s windows, 15-s overlap | Same | Direct reproduction. |
| Virtual-source geometries | Fixed top receiver with approximately 50-m receiver increments; separate constant-offset analysis | Channel 0 with approximately 50-m receiver increments is used here | Fixed-top geometry reproduced; channel 0 is treated as the top coordinate under the shared-fiber assumption. |
| Nearby-receiver stack | Target receiver plus R±10 channels; simple stack followed by 3.2-km/s local alignment | Both simple and locally aligned versions shown in Figure 1 | Direct reproduction; the alignment is explicitly labelled as velocity-conditioned. |
| Display band | 5–20 Hz after correlation stack | Same | Direct reproduction. |
| F–K use | Reported for earthquake/VSP processing, not the ambient-interferometry subsection | Applied before ambient correlation as the tested extension | All F–K conclusions require matched controls and input-level nulls. |

The source for the 2017 acquisition and processing entries is [Lellouch et al. (2019)](https://doi.org/10.1029/2019JB017533). Values for the present archive are read from the current HDF5 headers and frozen products.
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

markdown("""## 1. Lellouch-style reproduction — no F–K filter

This is the required baseline and control. The F–K filter is not applied here. The dashed lines are comparison trajectories; they are not fitted velocities. The geometry, windowing, nearby-receiver stack, final 5–20 Hz filter, and post-stack 3.2 km/s alignment follow the ambient-interferometry description of [Lellouch et al. (2019)](https://doi.org/10.1029/2019JB017533). The data epoch and interrogator are different, and the present channel-0 registration is assumed to correspond to the same top coordinate on the same cemented main-hole fiber; those are experiment-specific caveats, not facts supplied by the 2019 paper.
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

markdown("""**Figure 1. Unfiltered Lellouch-style ambient-noise correlation gathers before and after local velocity alignment.** Correlations between the channel-0 virtual source and receivers at approximately 50-m increments are computed from 1,440 one-minute records using a strain-rate proxy, running-absolute-mean temporal normalization, 30-s windows with 15-s overlap, and a final 5–20 Hz bandpass. The left panel averages each target receiver with its R±10 neighboring channels without time shifts. The right panel first shifts those neighboring traces by their local offset divided by 3.2 km/s and then averages them, reproducing the paper's alignment step. Both panels share one color scale. Receiver offset along the fiber is plotted vertically and correlation lag horizontally; black and gray trajectories show constant apparent velocities of 3.2, 2.5, and 4.0 km/s for comparison, not fitted travel times. These are no-F–K controls. The right panel is already conditioned locally on 3.2 km/s, so any visual enhancement there cannot serve as an independent velocity measurement.""")

markdown("""## 2. Raw → strain rate → temporal normalization → F–K → correlation

This cell follows one real record through the operations in their actual order. [Lellouch et al. (2019)](https://doi.org/10.1029/2019JB017533) convert strain to a strain-rate proxy, use running-absolute-mean normalization, divide the continuous record into 30-s windows with 15-s overlap, correlate channel 0 with receivers, stack nearby receivers, and finally filter the correlations from 5–20 Hz. Temporal normalization suppresses the leverage of earthquakes and other transients; it does not remove those waveforms ([Bensen et al., 2007](https://doi.org/10.1111/j.1365-246X.2007.03374.x)). The F–K extension is inserted after normalization and before correlation, so its effect must be tested against the unfiltered path with all other operations held fixed.
""")
python("""try:
    import sys, h5py, pandas as pd
    sys.path.insert(0, str(ROOT))
    from scipy.ndimage import uniform_filter1d
    from scipy.signal import detrend
    from ambient_transfer_test import CSV, corrected_path
    db = pd.read_csv(CSV, sep=r'\\s+'); db = db[db.nSamples > 0]
    path = next((corrected_path(p) for p in db.file if corrected_path(p).exists()), None)
    if path is None: raise FileNotFoundError('No shared HDF5 record mounted')
    with h5py.File(path,'r') as handle:
        dataset=handle['Acquisition/Raw[0]/RawData']
        fs=float(dataset.attrs.get('OutputDataRate',500.0)); dx=float(handle['Acquisition'].attrs.get('SpatialSamplingInterval',1.0))
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

Ten independent white-noise ensembles, each comprising three 12-s arrays, were passed through detrending, 5–20 Hz filtering, 5-s running-absolute-mean normalization, 2× decimation, each frozen velocity/direction mask, channel-0 correlations, and the same moveout score. The question is quantitative: does the held-out real-data score exceed the 95th percentile of identically processed white noise? The synthetic duration is much shorter than the real stack and therefore gives a conservative, relatively broad finite-noise distribution rather than a duration-matched p-value. This is a project-specific surrogate-data falsification test: surrogate methods test whether a statistic exceeds what the same analysis produces after destroying the property of interest ([Theiler et al., 1992](https://doi.org/10.1016/0167-2789(92)90102-S)). Passing white noise is necessary but not sufficient because white noise does not preserve the real data's single-channel spectra or nonstationarity.
""")
python("""white=np.load(products['white_noise']); white_meta=json.loads(products['white_noise_meta'].read_text()); wi=np.argmin(abs(white['velocities_m_s']-3200))
fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
for ax,(mask,label) in zip(axes.flat,masks):
    for branch,color in [('negative','tab:blue'),('positive','tab:orange')]:
        values=np.abs(white[f'{mask}__{branch}__ensemble_scores'][:,wi]); real_score=abs(fk[f'held_out__{mask}__{branch}__scores'][wi])
        ax.hist(values,bins=7,alpha=.45,color=color,label=f'{branch} white'); ax.axvline(real_score,color=color,lw=2)
    ax.set(title=label,xlabel='Absolute score at 3.2 km/s',ylabel='White-noise ensembles'); ax.legend(frameon=False)
plt.show()
for mask,label in masks:
    for branch in ['negative','positive']:
        result=white_meta['masks'][mask][branch]
        print(label,branch,'real exceeds white95:',result['real_exceeds_white_95_at_3200'])
""")

markdown("""**Figure 8. White-noise full-operator null at 3.2 km/s.** Histograms show absolute moveout scores from 10 independent Gaussian white-noise ensembles, each built from three 12-s records and passed through the production preprocessing, decimation, signed F–K mask, channel-0 correlation geometry, and scoring operator. Vertical lines show the corresponding held-out real-data score for each signed branch and mask. The production and narrow fans exceed their own white-noise 95th percentiles, whereas the broad and direction-only filters do not recover the held-out real score. The narrow fan nevertheless gives white noise a median absolute 3.2-km/s score of approximately 0.64, demonstrating that the absolute score is strongly conditioned by the imposed mask; only a comparison with the same mask is meaningful. The synthetic stacks are shorter than the real stack, so this is a conservative implementation gate rather than a duration-matched significance estimate. It shows that finite white-noise input does not routinely reproduce the real production-fan score, but it is not sufficient validation because white noise lacks the spectra, coherent artifacts, and nonstationarity of the measured data.""")

markdown("""## 6. Pre-filter channel-scramble and time-shift nulls

These are stronger than shuffling the finished correlation panel. For every real input record, channel 0 is fixed and all other preprocessed traces are reassigned to random fiber coordinates before F–K filtering. A second null independently circularly shifts each non-source trace before filtering. Both preserve important attributes of the measured traces while destroying ordered interchannel propagation. The nulls are tailored to this analysis; [Theiler et al. (1992)](https://doi.org/10.1016/0167-2789(92)90102-S) provides the general surrogate-test logic, not this exact geophysical implementation.

**Scope of this null product.** It uses the same 5–20 Hz preprocessing, 5-s normalization, 2× decimation, signed fan, channel-0 virtual source, 50-m receiver spacing, and moveout statistic as the production F–K analysis. It correlates each one-minute record directly and does not reproduce the Lellouch 30-s/15-s-overlap and nearby-receiver-stack details used in Section 4. Because this product fails rather than passes, it remains a valid warning that the fan can generate the selected statistic after spatial order is destroyed; it is not presented as an exact matched-pipeline p-value.
""")
python("""null=np.load(products['fk_nulls']); vv=null['velocities_m_s']; i=np.argmin(abs(vv-3200)); observed=abs(float(null['observed_negative_scores'][i]))
cp=np.abs(null['null_channel_permutation_negative_scores'][:,i]); ts=np.abs(null['null_circular_time_shift_negative_scores'][:,i])
fig,ax=plt.subplots(figsize=(8,4)); ax.hist(cp,bins=10,alpha=.6,label='pre-filter channel scramble'); ax.hist(ts,bins=10,alpha=.6,label='pre-filter time shifts'); ax.axvline(observed,color='k',lw=2,label=f'real={observed:.3f}'); ax.set(xlabel='Absolute 3.2 km/s score',ylabel='Realizations',title='Production-fan full-pipeline nulls'); ax.legend(frameon=False); plt.show()
print('Real score:',observed,'channel-scramble 95%:',np.quantile(cp,.95),'time-shift 95%:',np.quantile(ts,.95))
print('Decision: FAIL — the real production-fan score does not exceed either input-level null.')
""")

markdown("""**Figure 9. Pre-filter spatial-order and phase-coherence surrogate tests.** The black line is the absolute 3.2-km/s score of the measured production-fan stack. Histograms show scores obtained when, before F–K filtering, non-source channels are randomly reassigned to fiber coordinates or independently circularly shifted in time while channel 0 is fixed. Both nulls preserve measured single-channel content while destroying the ordered interchannel relationship needed for physical moveout. The observed score lies below both 95th-percentile null thresholds and has an exceedance probability of 1.0 in the stored ensembles. The selected fan therefore fails this QC gate: it can produce an equal or stronger moveout statistic after the physical spatial or phase ordering is destroyed.""")

markdown("""## 7. Is F–K required, and is the recovered feature defensible?

The phrase “F–K is required” has two meanings. The visual statement is true: a velocity-bounded fan makes an ordered ridge appear while the unfiltered and direction-only paths remain weak near 3.2 km/s. The scientific statement is not yet true: in the matched five-hour calculation, opposite signed fans produce nearly identical selected-lag curves and the same 2.675-km/s peak, and the production fan gives an equal or stronger ridge after channel scrambling and time shifting. Therefore the fan is sufficient to impose an ordered kernel, and the current result cannot be accepted as independent signal recovery. This conservative decision is consistent with Green's-function theory: an attractive cross-correlation is not automatically a Green's function when illumination and source-boundary assumptions are unverified ([Wapenaar and Fokkema, 2006](https://doi.org/10.1190/1.2213955)).
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

markdown("""## 9. Requirement-by-requirement QC decision

1. **Lellouch reproduction:** implemented with the published top-source/50-m geometries, RAM normalization, 30-s windows, 15-s overlap, nearby-receiver stack, 3.2 km/s alignment only after simple stacking, and final 5–20 Hz filtering.
2. **Is F–K required?** It is required to make the present ridge visible, but the ridge is not accepted because the same mask fails input-level surrogate tests.
3. **Velocity and direction combinations:** production, narrow, broad, and direction-only masks are shown for verified downgoing and upgoing branches.
4. **White noise:** passed through the complete production-fan operator and correlation score; the fan clears this necessary gate. This is not an exact Lellouch-window null.
5. **Channel scrambling:** applied before F–K filtering; the production fan fails this decisive gate. The pipeline departures are stated in Section 6.
6. **Processing progression:** raw strain, strain rate, RAM normalization, pre-filter F–K power, filtered wavefield, unfiltered correlation, matched filtered correlations, and decimation sensitivity are shown in order.
7. **Minimum stack duration:** the fan product appears stable after about one day, but this is operator stability. A defensible physical-signal convergence time is **not established within the eight tested days**.
8. **Advisor confidence:** every positive-looking result is placed beside the null that determines whether it is accepted.

**Current decision:** F–K filtering is implemented correctly, but the 2.5–4.5 km/s ambient ridge does not yet pass the full QC ladder. It must not be presented as a recovered P-wave Green’s function or formation velocity.
""")

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
