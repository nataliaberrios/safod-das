#!/usr/bin/env python3
"""Memory-safe Lellouch-style ambient interferometry transfer test."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
from scipy.signal import butter, sosfiltfilt, detrend
from scipy.ndimage import uniform_filter1d
import matplotlib.pyplot as plt

CSV = Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv")
OUT = Path(__file__).resolve().parent / "ambient_transfer"

def corrected_path(p):
    return Path(str(p).replace("/data/SAFODAS1-harddrive-transfer", "/data/SAFOD/SAFODAS1-harddrive-transfer"))

def load_segment(path):
    with h5py.File(path, "r") as h:
        raw_group = h["Acquisition/Raw[0]"]
        dataset = raw_group["RawData"]
        d = dataset[:].astype(np.float32)
        a = h["Acquisition"].attrs
        fs = float(raw_group.attrs.get("OutputDataRate", dataset.attrs.get("OutputDataRate", 500.0)))
        dx = float(a.get("SpatialSamplingInterval", 1.0))
    return d.T, fs, dx

def preprocess(x, fs, fmin=5.0, fmax=20.0, norm_seconds=5.0):
    x = detrend(x, axis=1, type="linear")
    sos = butter(4, [fmin, fmax], btype="bandpass", fs=fs, output="sos")
    x = sosfiltfilt(sos, x, axis=1)
    nwin = max(3, int(norm_seconds * fs))
    kernel = np.ones(nwin, dtype=np.float32) / nwin
    # Vectorized running absolute mean; equivalent to the prior loop but much faster.
    m = uniform_filter1d(np.abs(x), size=nwin, axis=1, mode="nearest")
    floor = np.percentile(m, 5, axis=1, keepdims=True) * 0.1 + 1e-12
    return x / np.maximum(m, floor)

def normalized_corr_pairs(x, pairs, fs, max_lag=0.35, batch=64):
    n=x.shape[1]; nfft=1<<int(np.ceil(np.log2(2*n-1)))
    ml=int(round(max_lag*fs))
    lags=np.arange(-ml,ml+1)/fs; out=np.zeros((len(pairs),len(lags)),float)
    # Keep FFT work bounded: the previous all-channel FFT was memory intensive.
    for first in range(0,len(pairs),batch):
        pp=pairs[first:first+batch]
        si=np.asarray([i for i,j in pp],dtype=int); sj=np.asarray([j for i,j in pp],dtype=int)
        fi=np.fft.rfft(x[si],n=nfft,axis=1); fj=np.fft.rfft(x[sj],n=nfft,axis=1)
        for k,((i,j),a,b) in enumerate(zip(pp,fi,fj)):
            # Non-negative lags begin at index zero; negative lags wrap to
            # the end of the padded circular-correlation array.
            c=np.fft.irfft(np.conj(a)*b,n=nfft)
            c=np.concatenate((c[-ml:],c[:ml+1]))
            den=np.sqrt(np.sum(x[i]**2)*np.sum(x[j]**2))
            if den>0: out[first+k]=c/den
        del fi,fj
    return lags,out

def stack_day(day,nfiles,norm_seconds,start=0):
    rows=day.iloc[start:start+nfiles]
    source_targets=[int(round(50*k/1.0209523)) for k in range(1,15)]
    source_targets=[c for c in source_targets if c<900]
    top_pairs=[(0,c) for c in source_targets]
    offset=int(round(50/1.0209523))
    fixed_pairs=[(i,i+offset) for i in range(0,900-offset,5)]
    top_stack=fixed_stack=None; lags=None; used=[]
    for ii,row in enumerate(rows.itertuples(index=False),1):
        f=corrected_path(row.file)
        if not f.exists(): continue
        try:
            x,fs,dx=load_segment(f); x=preprocess(x,fs,norm_seconds=norm_seconds)
            lags,top=normalized_corr_pairs(x,top_pairs,fs); _,fixed=normalized_corr_pairs(x,fixed_pairs,fs)
            top_stack=top if top_stack is None else top_stack+top
            fixed_stack=fixed if fixed_stack is None else fixed_stack+fixed
            used.append(str(f))
            if ii%10==0: print(f"processed {ii}/{nfiles}; usable={len(used)}",flush=True)
        except Exception as e: print("skip",f,repr(e),flush=True)
    if not used: raise RuntimeError("No usable files")
    return lags,np.asarray(source_targets)*dx,top_stack/len(used),fixed_stack/len(used),used,dx,fs

def velocity_score(top,lags,dist,v):
    vals=[]
    for row,d in zip(top,dist):
        t=d/v
        k=int(np.argmin(np.abs(lags-t))) if abs(t)<=lags.max() else None
        vals.append(np.nan if k is None else row[k])
    return float(np.nanmedian(vals)),float(np.nanmean(vals)),float(np.nanstd(vals))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--date",default="2025-02-20"); ap.add_argument("--nfiles",type=int,default=10); ap.add_argument("--start",type=int,default=0); ap.add_argument("--norm-seconds",type=float,default=5.0); args=ap.parse_args()
    db=pd.read_csv(CSV,sep=r"\s+"); db=db[db.nSamples>0].copy(); db["t"]=pd.to_datetime(db.startTime,utc=True,errors="coerce")
    day=db[db.t.dt.strftime("%Y-%m-%d")==args.date].sort_values("t").reset_index(drop=True)
    OUT.mkdir(exist_ok=True); lags,dist,top,fixed,used,dx,fs=stack_day(day,args.nfiles,args.norm_seconds,args.start)
    scores={str(v):velocity_score(top,lags,dist,v) for v in [1500,2000,2500,3000,3200,3400,4000,5000]}
    stem=f"transfer_{args.date}_start{args.start}_n{len(used)}"; np.savez(OUT/(stem+".npz"),lags=lags,distances=dist,top_stack=top,fixed_stack=fixed,used=np.array(used),dx=dx,fs=fs)
    fig,axs=plt.subplots(1,3,figsize=(16,5),constrained_layout=True); extent=[lags[0],lags[-1],dist[-1],dist[0]]
    axs[0].imshow(top,extent=extent,aspect="auto",cmap="RdBu_r",vmin=np.percentile(top,2),vmax=np.percentile(top,98)); axs[0].set_title(f"Top source to receivers ({len(used)} files)"); axs[0].set_xlabel("Lag (s)"); axs[0].set_ylabel("Receiver position relative to top source (m)")
    for v in [2000,2500,3000,3200,3400,4000]: axs[0].plot(dist/v,dist,"k--",lw=.7,alpha=.6)
    vv=np.linspace(1200,6000,97); axs[1].plot(vv,[velocity_score(top,lags,dist,v)[0] for v in vv]); axs[1].axvline(3200,color="r",ls="--"); axs[1].set_xlabel("Trial velocity (m/s)"); axs[1].set_ylabel("Median normalized correlation")
    offset=int(round(50/dx))
    pair_midpoints=np.asarray([(i+(i+offset))/2 for i in range(0,900-offset,5)])*dx
    axs[2].imshow(fixed,extent=[lags[0],lags[-1],pair_midpoints[-1],pair_midpoints[0]],aspect="auto",cmap="RdBu_r",vmin=np.percentile(fixed,2),vmax=np.percentile(fixed,98)); axs[2].set_title(f"Sliding 50-m pairs ({fixed.shape[0]} pairs)"); axs[2].set_xlabel("Lag (s)"); axs[2].set_ylabel("Pair midpoint along fiber (m)")
    fig.savefig(OUT/(stem+".png"),dpi=220)
    report={"date":args.date,"requested_files":args.nfiles,"used_files":len(used),"fs_hz":fs,"dx_m":dx,"norm_seconds":args.norm_seconds,"scores":scores}; (OUT/(stem+".json")).write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=="__main__": main()
