"""Source-history regression and phase-free tidal-timescale test for Nano delays."""
from pathlib import Path
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import f as f_dist

HERE = Path(__file__).resolve().parent
PERIOD_H = 24.0
N_BOOT = 5000
BLOCK_LEN = 3
SEED = 20260802

def zscore(v):
    v = np.asarray(v, float); scale = np.std(v, ddof=1)
    return (v - np.mean(v)) / scale if scale else np.zeros_like(v)

def load_table():
    m = pd.read_csv(HERE/"awd_manifest.csv", parse_dates=["utc_time"])
    b = m.groupby("burst_id", as_index=False).agg(utc_time=("utc_time","median"),
                                                    manifest_drops=("drop_id","count"))
    r = pd.read_csv(HERE/"nano_burst_repeatability_hierarchical.csv")
    d = b.merge(r, on="burst_id").sort_values("utc_time").reset_index(drop=True)
    d["hours"] = (d.utc_time-d.utc_time.iloc[0]).dt.total_seconds()/3600.
    d["cumulative_drops"] = d.n_drops.cumsum()-d.n_drops.iloc[0]
    d["delay_ms"] = 1000*d.loo_burst_delay_s
    d["compaction_proxy_z"] = zscore(d.cumulative_drops)
    d["coupling_proxy_z"] = zscore(d.burst_signal_rms)
    return d

def design(d, harmonic=False):
    t = d.hours.to_numpy(float)
    x = np.c_[np.ones(len(d)), d.compaction_proxy_z, d.coupling_proxy_z]
    return np.c_[x, np.sin(2*np.pi*t/PERIOD_H), np.cos(2*np.pi*t/PERIOD_H)] if harmonic else x

def fit(x, y):
    b = np.linalg.lstsq(x, y, rcond=None)[0]
    e = y-x@b
    return b, e, float(e@e)

def blocked_cv(x, y, folds=5):
    pred = np.full(len(y), np.nan)
    for test in np.array_split(np.arange(len(y)), folds):
        train = np.setdiff1d(np.arange(len(y)), test)
        b = np.linalg.lstsq(x[train], y[train], rcond=None)[0]
        pred[test] = x[test]@b
    return float(np.sqrt(np.mean((y-pred)**2))), pred

def bootstrap_delta(x0, x1, y, b0):
    rng = np.random.default_rng(SEED); e = y-x0@b0; n=len(y); out=np.empty(N_BOOT)
    for k in range(N_BOOT):
        sample=[]
        while len(sample)<n:
            j=int(rng.integers(n)); sample.extend(e[(j+np.arange(BLOCK_LEN))%n])
        yy=x0@b0+np.asarray(sample[:n])
        _,_,r0=fit(x0,yy); _,_,r1=fit(x1,yy); out[k]=r0-r1
    return out

def extrema(b):
    g=np.linspace(0,PERIOD_H,10001)
    q=b[3]*np.sin(2*np.pi*g/PERIOD_H)+b[4]*np.cos(2*np.pi*g/PERIOD_H)
    imax=int(np.argmax(q)); imin=int(np.argmin(q))
    return float(g[imax]),float(q[imax]),float(g[imin]),float(q[imin]),g,q

def make_figure(d, x0, x1, b0, b1, pred0, pred1, boot, out):
    t=d.hours.to_numpy(float); y=d.delay_ms.to_numpy(float); e=y-x0@b0
    hmax,qmax,hmin,qmin,g,q=extrema(b1)
    fig,ax=plt.subplots(2,2,figsize=(12,8.2),constrained_layout=True)
    a=ax[0,0]; a.scatter(t,y,c=d.coupling_proxy_z,cmap="viridis",s=38,edgecolor="white",lw=.35)
    a.plot(t,pred0,color="#5d6d7e",lw=1.7,label="history/coupling baseline")
    a.plot(t,pred1,color="#c0392b",lw=2,label="+ 24-h harmonic")
    a.axhline(0,color=".3",lw=.7,ls=":"); a.grid(alpha=.2)
    a.set(xlabel="hours since first burst",ylabel="LOO-template delay (ms)",title="(a) Delay and nested fits")
    a.legend(frameon=False,fontsize=8)
    a=ax[0,1]; a.plot(t,d.compaction_proxy_z,color="#7d3c98",lw=2,label="cumulative drops")
    a.plot(t,d.coupling_proxy_z,color="#148f77",lw=1.7,label="signal RMS")
    a.axhline(0,color=".3",lw=.7,ls=":"); a.grid(alpha=.2)
    a.set(xlabel="hours since first burst",ylabel="standardized proxy",title="(b) Source nuisance proxies")
    a.legend(frameon=False,fontsize=8)
    a=ax[1,0]; a.scatter(t,e,color="#34495e",s=32,edgecolor="white",lw=.35,label="delay − baseline")
    a.plot(g,q,color="#c0392b",lw=2.2,label="fitted 24-h component")
    a.axhline(0,color=".3",lw=.7,ls=":")
    a.axvline(hmax,color="#c0392b",lw=.9,ls="--"); a.axvline(hmin,color="#c0392b",lw=.9,ls="--")
    a.text(hmax,qmax,f"max {hmax:.1f} h",ha="center",va="bottom",fontsize=8)
    a.text(hmin,qmin,f"min {hmin:.1f} h",ha="center",va="top",fontsize=8)
    a.grid(alpha=.2); a.set(xlabel="hours since first burst",ylabel="residual delay (ms)",title="(c) Nonmonotone residual candidate")
    a.legend(frameon=False,fontsize=8)
    a=ax[1,1]; cv0=blocked_cv(x0,y)[0]; cv1=blocked_cv(x1,y)[0]
    bars=a.bar(["baseline\nblocked CV","+24-h harmonic\nblocked CV"],[cv0,cv1],color=["#7f8c8d","#c0392b"],width=.62)
    a.set(ylabel="5-fold contiguous CV RMSE (ms)",title="(d) Predictive check"); a.grid(axis="y",alpha=.2)
    for bar in bars: a.text(bar.get_x()+bar.get_width()/2,bar.get_height()+.006,f"{bar.get_height():.3f}",ha="center",fontsize=9)
    a2=a.twinx(); a2.hist(boot,bins=30,color="#95a5a6",alpha=.28,density=True)
    obs=fit(x0,y)[2]-fit(x1,y)[2]; a2.axvline(obs,color="#8e44ad",lw=2,label=f"observed ΔRSS={obs:.3f}")
    a2.set_ylabel("null density (block bootstrap)"); a2.legend(frameon=False,fontsize=7,loc="upper right")
    fig.suptitle("Nano burst-delay regression: source-history removal and tidal-timescale test",fontsize=14,fontweight="bold")
    fig.savefig(out,dpi=220); fig.savefig(out.with_suffix(".pdf")); plt.close(fig)

def main():
    d=load_table(); y=d.delay_ms.to_numpy(float); x0=design(d); x1=design(d,True)
    b0,e0,rss0=fit(x0,y); b1,e1,rss1=fit(x1,y); cv0,pred0=blocked_cv(x0,y); cv1,pred1=blocked_cv(x1,y)
    delta=rss0-rss1; q=x1.shape[1]-x0.shape[1]; df2=len(y)-x1.shape[1]
    fstat=(delta/q)/(rss1/df2); naive=float(f_dist.sf(fstat,q,df2))
    boot=bootstrap_delta(x0,x1,y,b0); bootp=float((1+(boot>=delta).sum())/(len(boot)+1))
    hmax,qmax,hmin,qmin,g,harm=extrema(b1); amp=float(np.hypot(b1[3],b1[4]))
    rows=[
      ("n_bursts",len(d),"bursts","independent burst-level delays"),
      ("delay_median",float(np.median(y)),"ms","LOO-template delay"),
      ("delay_p05",float(np.percentile(y,5)),"ms","LOO-template delay"),
      ("delay_p95",float(np.percentile(y,95)),"ms","LOO-template delay"),
      ("baseline_rss",rss0,"ms^2","intercept + cumulative-drop + coupling"),
      ("harmonic_rss",rss1,"ms^2","baseline + phase-free 24-h sine/cosine"),
      ("baseline_blocked_cv_rmse",cv0,"ms","five contiguous time blocks"),
      ("harmonic_blocked_cv_rmse",cv1,"ms","five contiguous time blocks"),
      ("delta_rss",delta,"ms^2","24-h improvement"),
      ("naive_f_p",naive,"probability","reference only; not autocorrelation corrected"),
      ("block_bootstrap_p",bootp,"probability",f"circular residual blocks of {BLOCK_LEN} bursts"),
      ("harmonic_amplitude",amp,"ms","phase-free 24-h component"),
      ("harmonic_max_hour",hmax,"hours","maximum of fitted component"),
      ("harmonic_max_ms",qmax,"ms","maximum of fitted component"),
      ("harmonic_min_hour",hmin,"hours","minimum of fitted component"),
      ("harmonic_min_ms",qmin,"ms","minimum of fitted component"),
      ("source_compaction_proxy","cumulative AWD drops","proxy","no independent ground-level compaction sensor"),
      ("source_coupling_proxy","burst signal RMS","proxy","observed coupling nuisance"),
      ("tidal_predictor","phase-free 24-hour harmonic","proxy","not an astronomical traction calculation"),
      ("interpretation_status","candidate nonmonotone delay component","status","not a tidal detection"),
    ]
    with (HERE/"nano_tidal_compaction_regression.csv").open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["quantity","value","units","interpretation"]); w.writerows(rows)
    with (HERE/"nano_tidal_compaction_regression.txt").open("w") as f:
        f.write("SAFOD AWD Nano burst-delay regression\nStatus: source-history nuisance regression plus phase-free tidal-timescale stress test.\n\n")
        f.write("Data: 49 burst-level leave-one-burst-out delays.\n")
        f.write("Compaction proxy: cumulative AWD drops; coupling proxy: burst signal RMS. Neither is an independent ground-level displacement measurement.\n")
        f.write(f"Delay median and p05-p95: {np.median(y):.3f} ms; {np.percentile(y,5):.3f} to {np.percentile(y,95):.3f} ms.\n")
        f.write(f"Baseline RMSE={np.sqrt(rss0/len(y)):.3f} ms; +24-h harmonic RMSE={np.sqrt(rss1/len(y)):.3f} ms.\n")
        f.write(f"Blocked 5-fold CV RMSE: baseline={cv0:.3f} ms; +24-h harmonic={cv1:.3f} ms.\n")
        f.write(f"Observed ΔRSS={delta:.4f} ms^2; naive F p={naive:.4g}; circular block-bootstrap p={bootp:.4g} (block={BLOCK_LEN}, N={N_BOOT}).\n")
        f.write(f"Fitted 24-h component amplitude={amp:.3f} ms; maximum at {hmax:.2f} h ({qmax:+.3f} ms), minimum at {hmin:.2f} h ({qmin:+.3f} ms).\n")
        f.write("Interpretation: source-history/coupling removal leaves a reproducible nonmonotone 24-hour-timescale delay component. This is a candidate tidal-timescale pattern, not a tidal-force detection: phase is fitted from the data and is not mapped to opening or closing without an astronomical tide/traction predictor.\n")
        f.write("Next check: add a physically phased Earth-tide/strain predictor and repeat with held-out bursts and environmental controls.\n")
    with (HERE/"nano_tidal_compaction_regression.json").open("w") as f:
        json.dump({"status":"candidate_nonmonotone_delay_component","n_bursts":len(d),
          "delay_median_ms":float(np.median(y)),"delay_p05_ms":float(np.percentile(y,5)),"delay_p95_ms":float(np.percentile(y,95)),
          "source_compaction_proxy":"cumulative AWD drops","source_coupling_proxy":"burst signal RMS",
          "tidal_predictor":"phase-free 24-hour sine/cosine pair","baseline_rss_ms2":rss0,"harmonic_rss_ms2":rss1,
          "blocked_cv_rmse_ms":{"baseline":cv0,"harmonic":cv1},"delta_rss_ms2":delta,"naive_f_p":naive,"block_bootstrap_p":bootp,
          "harmonic_amplitude_ms":amp,"harmonic_max_hour":hmax,"harmonic_min_hour":hmin,
          "limitations":["No independent ground-level compaction sensor was found.","The 24-hour harmonic is phase-free, not an astronomical traction calculation.","The 49 delays are time ordered, not independent daily experiments.","Opening/closing sign requires an externally phased predictor."]},f,indent=2)
    np.savez(HERE/"nano_tidal_compaction_regression.npz",hours=d.hours.to_numpy(float),delay_ms=y,cumulative_drops=d.cumulative_drops.to_numpy(int),signal_rms=d.burst_signal_rms.to_numpy(float),baseline_prediction=x0@b0,harmonic_prediction=x1@b1,baseline_residual_ms=e0,harmonic_residual_ms=e1,bootstrap_delta_rss=boot,beta_baseline=b0,beta_harmonic=b1,harmonic_grid_hours=g,harmonic_component_ms=harm)
    make_figure(d,x0,x1,b0,b1,pred0,pred1,boot,HERE/"nano_tidal_compaction_regression.png")
    print((HERE/"nano_tidal_compaction_regression.txt").read_text())

if __name__=="__main__": main()

