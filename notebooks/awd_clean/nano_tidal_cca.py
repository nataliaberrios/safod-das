"""Canonical-correlation test for the source-corrected Nano delay pattern."""
from pathlib import Path
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PERIODS_H = (24.0, 12.4206)
BLOCK_LEN = 3
N_PERM = 5000
SEED = 20260802

def z(v):
    v=np.asarray(v,float); return (v-v.mean())/v.std(ddof=1)

def load():
    m=pd.read_csv(HERE/"awd_manifest.csv",parse_dates=["utc_time"])
    b=m.groupby("burst_id",as_index=False).agg(utc_time=("utc_time","median"))
    r=pd.read_csv(HERE/"nano_burst_repeatability_hierarchical.csv")
    d=b.merge(r,on="burst_id").sort_values("utc_time").reset_index(drop=True)
    d["hours"]=(d.utc_time-d.utc_time.iloc[0]).dt.total_seconds()/3600.
    d["cumulative_drops"]=d.n_drops.cumsum()-d.n_drops.iloc[0]
    d["delay_ms"]=1000*d.loo_burst_delay_s
    B=np.c_[np.ones(len(d)),z(d.cumulative_drops),z(d.burst_signal_rms)]
    d["delay_source_corrected_ms"]=d.delay_ms-(B@np.linalg.lstsq(B,d.delay_ms,rcond=None)[0])
    return d

def matrices(d):
    t=d.hours.to_numpy(float)
    X=np.c_[np.sin(2*np.pi*t/PERIODS_H[0]),np.cos(2*np.pi*t/PERIODS_H[0]),
            np.sin(2*np.pi*t/PERIODS_H[1]),np.cos(2*np.pi*t/PERIODS_H[1])]
    Y=np.c_[z(d.delay_source_corrected_ms),z(d.burst_signal_rms)]
    return X,Y

def invsqrt(S):
    w,v=np.linalg.eigh(S)
    return (v*(1/np.sqrt(np.maximum(w,1e-8))))@v.T

def cca(X,Y):
    xm=X.mean(0); ym=Y.mean(0); xx=X-xm; yy=Y-ym; n=len(X)
    sxx=xx.T@xx/(n-1); syy=yy.T@yy/(n-1); sxy=xx.T@yy/(n-1)
    A=invsqrt(sxx)@sxy@invsqrt(syy)
    u,s,vt=np.linalg.svd(A,full_matrices=False)
    wx=invsqrt(sxx)@u[:,0]; wy=invsqrt(syy)@vt.T[:,0]
    sx=xx@wx; sy=yy@wy; corr=float(np.corrcoef(sx,sy)[0,1])
    if corr<0: wy=-wy; sy=-sy; corr=-corr
    return {"corr":corr,"wx":wx,"wy":wy,"xscore":sx,"yscore":sy,"singular_values":s,
            "xmean":xm,"ymean":ym}

def block_order(rng,n):
    starts=np.arange(0,n,BLOCK_LEN); out=[]
    for j in rng.permutation(len(starts)):
        out.extend((starts[j]+np.arange(BLOCK_LEN))%n)
    return np.asarray(out[:n])

def permutation_null(X,Y):
    rng=np.random.default_rng(SEED); n=len(Y); row=np.empty(N_PERM); block=np.empty(N_PERM)
    for k in range(N_PERM):
        row[k]=cca(X,Y[rng.permutation(n)])["corr"]
        block[k]=cca(X,Y[block_order(rng,n)])["corr"]
    return row,block

def blocked_cv(X,Y,nfold=5):
    values=[]
    for test in np.array_split(np.arange(len(X)),nfold):
        train=np.setdiff1d(np.arange(len(X)),test)
        fit=cca(X[train],Y[train])
        sx=(X[test]-fit["xmean"])@fit["wx"]; sy=(Y[test]-fit["ymean"])@fit["wy"]
        values.append(float(abs(np.corrcoef(sx,sy)[0,1])))
    return np.asarray(values)

def make_figure(d,fit,row,block,cv,out):
    t=d.hours.to_numpy(float)
    fig,ax=plt.subplots(2,2,figsize=(11.8,8.2),constrained_layout=True)
    a=ax[0,0]; a.plot(t,fit["xscore"],color="#c0392b",lw=2,label="tidal-basis canonical score")
    a.plot(t,fit["yscore"],color="#34495e",lw=1.8,label="delay/coupling canonical score")
    a.axhline(0,color=".3",lw=.7,ls=":"); a.grid(alpha=.2)
    a.set(xlabel="hours since first burst",ylabel="canonical score",title="(a) CCA scores versus time"); a.legend(frameon=False,fontsize=8)
    a=ax[0,1]; a.scatter(fit["xscore"],fit["yscore"],c=t,cmap="plasma",s=40,edgecolor="white",lw=.35)
    a.set(xlabel="tidal-basis score",ylabel="delay/coupling score",title=f"(b) Canonical variates (r={fit['corr']:.3f})"); a.grid(alpha=.2)
    a.text(.04,.94,"color = elapsed time",transform=a.transAxes,fontsize=8,va="top")
    a=ax[1,0]; labels=["24h sin","24h cos","12.42h sin","12.42h cos"]
    x=np.arange(4); a.bar(x-.18,fit["wx"],.36,color="#c0392b",label="tidal basis")
    a.bar(x+.18,[fit["wy"][0],fit["wy"][1],0,0],.36,color="#34495e",label="delay/coupling (delay, RMS)")
    a.set_xticks(x,labels,rotation=25,ha="right"); a.axhline(0,color=".3",lw=.7); a.grid(axis="y",alpha=.2)
    a.set(ylabel="CCA weight",title="(c) Canonical weights"); a.legend(frameon=False,fontsize=8)
    a=ax[1,1]; a.hist(row,bins=30,color="#95a5a6",alpha=.45,density=True,label="row permutation")
    a.hist(block,bins=30,color="#5d6d7e",alpha=.40,density=True,label="3-burst block permutation")
    a.axvline(fit["corr"],color="#8e44ad",lw=2,label=f"observed r={fit['corr']:.3f}")
    a.set(xlabel="first canonical correlation",ylabel="density",title="(d) Null tests"); a.grid(alpha=.2); a.legend(frameon=False,fontsize=7)
    fig.suptitle("Canonical correlation of source-corrected delay with tidal-timescale bases",fontsize=14,fontweight="bold")
    fig.savefig(out,dpi=220); fig.savefig(out.with_suffix(".pdf")); plt.close(fig)

def main():
    d=load(); X,Y=matrices(d); fit=cca(X,Y); row,block=permutation_null(X,Y); cv=blocked_cv(X,Y)
    p_row=float((1+(row>=fit["corr"]).sum())/(len(row)+1))
    p_block=float((1+(block>=fit["corr"]).sum())/(len(block)+1))
    with (HERE/"nano_tidal_cca.csv").open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["quantity","value","units","interpretation"])
        w.writerows([
          ("n_bursts",len(d),"bursts","burst-level observations"),
          ("canonical_correlation_1",fit["corr"],"r","full-sample CCA"),
          ("row_permutation_p",p_row,"probability","rows independently permuted"),
          ("block_permutation_p",p_block,"probability","3-burst blocks permuted"),
          ("blocked_cv_corr_mean",float(np.mean(cv)),"r","five contiguous test folds"),
          ("blocked_cv_corr_median",float(np.median(cv)),"r","five contiguous test folds"),
          ("blocked_cv_corr_values",";".join(f"{v:.3f}" for v in cv),"r","fold-by-fold values"),
          ("x_set","24 h and 12.4206 h sine/cosine bases","predictor","phase-free tidal-timescale basis"),
          ("y_set","source-corrected delay and burst signal RMS","response","delay corrected using history/coupling baseline"),
          ("status","exploratory CCA; not tidal detection","status","block null and CV are controlling checks"),
        ])
    with (HERE/"nano_tidal_cca.txt").open("w") as f:
        f.write("SAFOD AWD Nano canonical-correlation test\n")
        f.write("Status: exploratory pattern-discovery test; not a causal or tidal-force detection.\n\n")
        f.write("X set: phase-free 24-hour and 12.4206-hour sine/cosine bases.\n")
        f.write("Y set: source-corrected burst delay and burst signal RMS.\n")
        f.write(f"Full-sample first canonical correlation: {fit['corr']:.3f}.\n")
        f.write(f"Independent row-permutation p={p_row:.4f}; 3-burst block-permutation p={p_block:.4f}.\n")
        f.write(f"Five contiguous-fold test correlations: {', '.join(f'{v:.3f}' for v in cv)}; mean={np.mean(cv):.3f}, median={np.median(cv):.3f}.\n")
        f.write("Interpretation: CCA finds a full-sample association with the phase-free tidal-timescale basis, but the serial-dependence-aware block null and unstable blocked-CV correlations do not establish a robust tidal signal. Use CCA to prioritize a physically phased tide/traction predictor, not to assign opening or closing.\n")
    with (HERE/"nano_tidal_cca.json").open("w") as f:
        json.dump({"status":"exploratory_cca_not_tidal_detection","n_bursts":len(d),"canonical_correlation":fit["corr"],"row_permutation_p":p_row,"block_permutation_p":p_block,"blocked_cv_correlations":cv.tolist(),"x_weights":fit["wx"].tolist(),"y_weights":fit["wy"].tolist(),"x_set":"phase-free 24-hour and 12.4206-hour sine/cosine bases","y_set":"source-corrected delay and burst signal RMS","limitations":["CCA phase is free, not an astronomical tide phase.","No independent ground-compaction measurement is available.","Block null and blocked CV are not robustly significant.","Opening/closing sign requires an external fault-traction convention."]},f,indent=2)
    np.savez(HERE/"nano_tidal_cca.npz",hours=d.hours.to_numpy(float),delay_ms=d.delay_ms.to_numpy(float),source_corrected_delay_ms=d.delay_source_corrected_ms.to_numpy(float),x_score=fit["xscore"],y_score=fit["yscore"],x_weights=fit["wx"],y_weights=fit["wy"],row_null=row,block_null=block,blocked_cv=cv)
    make_figure(d,fit,row,block,cv,HERE/"nano_tidal_cca.png")
    print((HERE/"nano_tidal_cca.txt").read_text())

if __name__=="__main__": main()

