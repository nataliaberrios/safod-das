"""Burst-bootstrap Nano apparent-velocity and dispersion test.

No phase name is assigned. Times are relative to the ``UTC_Date`` column in
``p26.cc9.txt``; that value is used only as an alignment timestamp and is not
assumed to be an impact time, onset, or phase pick.
"""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt

HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
STEM = HERE / "nano_mode_identification"
PRE_S, APERTURE, STRIDE = 0.5, (80., 440.), 4
BANDS = np.array(((15.,25.),(25.,35.),(35.,45.),(45.,60.),(60.,80.)))
WIN_S = 0.120
TGRID = np.arange(-0.100, 0.3201, 0.004)
VGRID = np.arange(1500., 5000.1, 50.)
PPOS = 1 / VGRID[::-1]
PGRID = np.r_[-PPOS[::-1], PPOS]
NBOOT, SEED, BOOT_CHUNK = 999, 20260802, 32

def bp(x, fs, band):
    return sosfiltfilt(butter(4, band, btype="bandpass", fs=fs, output="sos"), x, axis=-1)

def moving_sum(x, n):
    c = np.r_[0., np.cumsum(x, dtype=np.float64)]
    return c[n:] - c[:-n]

def scan(section, x, fs):
    """Signed-slowness/intercept semblance; t0 is window start, not onset."""
    nw = int(round(WIN_S * fs)); tstart = TGRID[0]
    samples = np.arange(int(round((PRE_S+tstart)*fs)),
                        int(round((PRE_S+TGRID[-1]+WIN_S)*fs))+1)
    it = np.rint((TGRID-tstart)*fs).astype(int)
    out = np.full((len(PGRID), len(TGRID)), np.nan, np.float32)
    ich = np.arange(len(x))[:, None]
    for ip, p in enumerate(PGRID):
        indices = samples[None, :] + np.rint(x*p*fs).astype(int)[:, None]
        if indices.min() < 0 or indices.max() >= section.shape[1]: continue
        y = section[ich, indices]
        num = moving_sum(np.sum(y, axis=0)**2, nw)
        den = len(x) * moving_sum(np.sum(y*y, axis=0), nw)
        valid = it < len(num)
        out[ip, valid] = np.divide(num[it[valid]], den[it[valid]],
                                   out=np.zeros(valid.sum()), where=den[it[valid]]>0)
    return out

def signed_peak(surface, positive=True):
    use = PGRID > 0 if positive else PGRID < 0
    ids = np.flatnonzero(use); local = surface[ids]
    i, j = np.unravel_index(np.nanargmax(local), local.shape)
    return ids[i], j, float(local[i,j])

def boot_peaks(scores, weights):
    n = len(weights); p=np.full(n,np.nan); t=np.full(n,np.nan); margin=np.full(n,np.nan)
    pos, neg = PGRID>0, PGRID<0
    for first in range(0,n,BOOT_CHUNK):
        last=min(first+BOOT_CHUNK,n)
        means=np.tensordot(weights[first:last],scores,axes=(1,0))
        for k,s in enumerate(means,first):
            q=int(np.nanargmax(s[pos])); i,j=np.unravel_index(q,(pos.sum(),s.shape[1]))
            p[k]=PGRID[np.flatnonzero(pos)[i]]; t[k]=TGRID[j]
            margin[k]=np.nanmax(s[pos])-np.nanmax(s[neg])
    return p,t,margin

def pct(x): return np.nanpercentile(x,(2.5,50,97.5))

def main():
    with np.load(STACKS) as z:
        counts=z["n_common"]; active=np.flatnonzero(counts>0)
        stacks=z["nano_stacks"][active]; fs=float(z["fs"]); dx=float(z["dx_nano"])
    c0=int(np.ceil(APERTURE[0]/dx)); c1=min(stacks.shape[1],int(np.floor(APERTURE[1]/dx))+1)
    channels=np.arange(c0,c1,STRIDE)
    finite=np.all(np.isfinite(stacks[:,channels]),axis=(0,2))
    rms=np.sqrt(np.mean(stacks[:,channels]**2,axis=(0,2))); good=rms[np.isfinite(rms)&(rms>0)]
    channels=channels[finite & np.isfinite(rms) & (rms>np.median(good)*1e-6)]
    if len(channels)<32: raise RuntimeError(f"only {len(channels)} common channels pass QC")
    data=stacks[:,channels]; xabs=channels*dx; x=xabs-xabs[0]; nb=len(active)
    rng=np.random.default_rng(SEED)
    bcounts=np.vstack([rng.multinomial(nb,np.full(nb,1/nb)) for _ in range(NBOOT)])
    weights=bcounts/nb
    scores=np.empty((len(BANDS),nb,len(PGRID),len(TGRID)),np.float32)
    boots_p=[]; boots_t=[]; boots_m=[]; summary=[]; burst_rows=[]
    for ib,band in enumerate(BANDS):
        filtered=bp(data,fs,band)
        for j in range(nb):
            scores[ib,j]=scan(filtered[j],x,fs)
            ip,it,sp=signed_peak(scores[ib,j],True); im,jm,sm=signed_peak(scores[ib,j],False)
            burst_rows.append(dict(band_low_hz=band[0],band_high_hz=band[1],
                canonical_burst_index=int(active[j]),n_common_drops=int(counts[active[j]]),
                positive_slowness_ms_per_m=PGRID[ip]*1e3,positive_t0_s=TGRID[it],
                positive_semblance=sp,negative_slowness_ms_per_m=PGRID[im]*1e3,
                negative_t0_s=TGRID[jm],negative_semblance=sm,direction_margin=sp-sm))
        mean=np.nanmean(scores[ib],axis=0); ip,it,sp=signed_peak(mean,True); im,_,sm=signed_peak(mean,False)
        pb,tb,mb=boot_peaks(scores[ib],weights); boots_p.append(pb);boots_t.append(tb);boots_m.append(mb)
        plo,pmed,phi=pct(pb*1e3); vlo,vmed,vhi=pct(1/pb);tlo,tmed,thi=pct(tb)
        rows=burst_rows[-nb:]
        summary.append(dict(band_low_hz=band[0],band_high_hz=band[1],band_center_hz=band.mean(),
            population_slowness_ms_per_m=PGRID[ip]*1e3,population_speed_mps=1/PGRID[ip],
            population_t0_s=TGRID[it],population_positive_semblance=sp,
            population_negative_semblance=sm,population_direction_margin=sp-sm,
            fraction_bursts_positive_preferred=np.mean([r["direction_margin"]>0 for r in rows]),
            bootstrap_probability_positive_preferred=np.mean(mb>0),
            bootstrap_slowness_p2p5_ms_per_m=plo,bootstrap_slowness_median_ms_per_m=pmed,
            bootstrap_slowness_p97p5_ms_per_m=phi,bootstrap_speed_p2p5_mps=vlo,
            bootstrap_speed_median_mps=vmed,bootstrap_speed_p97p5_mps=vhi,
            bootstrap_t0_p2p5_s=tlo,bootstrap_t0_median_s=tmed,bootstrap_t0_p97p5_s=thi))
        print(f"finished {band[0]:.0f}-{band[1]:.0f} Hz")
    pb=np.column_stack(boots_p);tb=np.column_stack(boots_t);mb=np.column_stack(boots_m)
    fc=BANDS.mean(axis=1); f0=fc-fc.mean(); slope=np.sum(pb*f0,axis=1)/np.sum(f0*f0); sci=pct(slope*1e6)
    with STEM.with_suffix(".csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    with (HERE/"nano_mode_burst_objectives.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(burst_rows[0]));w.writeheader();w.writerows(burst_rows)

    fig,ax=plt.subplots(2,2,figsize=(13.2,10),constrained_layout=True); colors=plt.cm.viridis(np.linspace(.08,.92,len(BANDS)))
    vb=1/pb; vm=np.median(vb,axis=0);vci=np.percentile(vb,(2.5,97.5),axis=0)
    ax[0,0].errorbar(fc,vm/1e3,yerr=np.vstack((vm-vci[0],vci[1]-vm))/1e3,fmt="o-",color="k",capsize=3)
    ax[0,0].set(xlabel="frequency-band center (Hz)",ylabel="positive apparent speed (km/s)",title="A  Burst-bootstrap speed (95% interval)");ax[0,0].grid(alpha=.25)
    ax[0,0].text(.03,.04,f"slowness trend 95% interval:\n{sci[0]:.3f} to {sci[2]:.3f} µs m⁻¹ Hz⁻¹",transform=ax[0,0].transAxes,fontsize=9)
    xp=np.linspace(x[0],x[-1],100)
    for i,(band,c) in enumerate(zip(BANDS,colors)):
        q=tb[:,i,None]+pb[:,i,None]*xp;lo,med,hi=np.percentile(q,(2.5,50,97.5),axis=0)
        ax[0,1].plot(med,xabs[0]+xp,color=c,label=f"{band[0]:.0f}–{band[1]:.0f} Hz");ax[0,1].fill_betweenx(xabs[0]+xp,lo,hi,color=c,alpha=.12)
    ax[0,1].invert_yaxis();ax[0,1].set(xlabel="time relative to p26.cc9.txt UTC_Date (s)",ylabel="distance along Nano fiber (m)",title="B  Coherent-window trajectories (95% interval)");ax[0,1].legend(frameon=False,fontsize=8)
    margins=[[r["direction_margin"] for r in burst_rows if r["band_low_hz"]==b[0]] for b in BANDS]
    ax[1,0].boxplot(margins,positions=fc,widths=6,showfliers=False);ax[1,0].axhline(0,color=".25",ls="--")
    ax[1,0].set(xlabel="frequency-band center (Hz)",ylabel="positive − negative peak semblance",title="C  Signed ordering evidence by burst")
    pos=PGRID>0;va=1/PGRID[pos]/1e3
    for i,(b,c) in enumerate(zip(BANDS,colors)):
        ax[1,1].plot(va,np.nanmax(np.nanmean(scores[i],axis=0)[pos],axis=1),color=c,label=f"{b[0]:.0f}–{b[1]:.0f} Hz")
    ax[1,1].set(xlabel="positive apparent speed (km/s)",ylabel="equal-burst mean peak semblance",title="D  Frequency-dependent moveout objectives");ax[1,1].legend(frameon=False,fontsize=8)
    fig.suptitle("SAFOD AWD Nano: frequency-dependent moveout and burst-bootstrap uncertainty");fig.savefig(STEM.with_suffix(".png"),dpi=220)
    np.savez_compressed(STEM.with_suffix(".npz"),bands_hz=BANDS,p_grid_s_per_m=PGRID,t0_grid_s=TGRID,
        burst_scores=scores,bootstrap_slowness_s_per_m=pb,bootstrap_t0_s=tb,
        bootstrap_direction_margin=mb,bootstrap_slowness_slope_s_per_m_per_hz=slope,
        bootstrap_counts=bcounts,active_burst_indices=active,selected_channels=channels,
        coordinate_m=xabs,window_s=WIN_S,seed=SEED)
    resolved=not(sci[0]<=0<=sci[2]); direction=all(r["bootstrap_probability_positive_preferred"]>=.95 for r in summary)
    lines=["SAFOD AWD Nano mode-identification summary",
      "Times are relative to p26.cc9.txt UTC_Date, used only as an alignment timestamp.",
      "t0 is the coherence-window start at the first aperture channel, not an onset or phase pick.",
      f"Bursts={nb}; represented common drops={counts[active].sum()}; bootstrap resamples={NBOOT}; seed={SEED}.",
      f"Slowness-trend 2.5/50/97.5 percentiles: {sci[0]:.6f}, {sci[1]:.6f}, {sci[2]:.6f} microseconds m^-1 Hz^-1.",""]
    for r in summary: lines.append(f"{r['band_low_hz']:.0f}-{r['band_high_hz']:.0f} Hz: speed {r['bootstrap_speed_median_mps']:.0f} [{r['bootstrap_speed_p2p5_mps']:.0f},{r['bootstrap_speed_p97p5_mps']:.0f}] m/s; P(+ preferred)={r['bootstrap_probability_positive_preferred']:.3f}.")
    lines += ["",f"Positive ordering supported in every band: {direction}.",
      f"Monotonic frequency-dependent slowness resolved: {resolved}.",
      "Resolved dispersion would be compatible with, but not unique to, guided propagation.",
      "Unresolved dispersion would be necessary but insufficient evidence for direct body-wave propagation.",
      "This test alone cannot distinguish direct P energy from a weakly dispersive guided/coupled mode."]
    STEM.with_suffix(".txt").write_text("\n".join(lines)+"\n");print("\n".join(lines))

if __name__ == "__main__": main()
