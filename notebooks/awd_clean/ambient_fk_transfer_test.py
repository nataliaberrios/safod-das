import sys,json,argparse
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from ambient_transfer_test import corrected_path,load_segment,preprocess,normalized_corr_pairs
import matplotlib.pyplot as plt
CSV=Path("/oak/stanford/groups/ettore88/data/SAFOD/SAFODAS1-harddrive-transfer/SAFOD_2024_2025.csv")
OUT=Path(__file__).resolve().parent/"ambient_transfer"

def fk_filter(x,fs,dx,mode):
    xd=x[::2,::2]; dxf=2*dx; fsf=fs/2
    f=np.fft.fftfreq(xd.shape[1],1/fsf); k=np.fft.fftfreq(xd.shape[0],dxf)
    K,F=np.meshgrid(k,f,indexing="ij"); af=np.abs(F); ak=np.abs(K); v=af/np.maximum(ak,1e-12)
    if mode in ("positive","negative","both"):
        mask=(af>=5)&(af<=20)&(v>=2500)&(v<=4500)&(ak>0)
        if mode=="positive": mask &= (F*K)>0
        if mode=="negative": mask &= (F*K)<0
    elif mode=="complement":
        mask=(af>=5)&(af<=20)&(v>=4500)&(v<=6500)&(ak>0)&((F*K)>0)
    else:
        raise ValueError(mode)
    y=np.fft.ifft2(np.fft.fft2(xd)*mask).real
    return y,fsf,dxf

def velocity_scores(top,lags,dist):
    vv=np.linspace(1200,6000,193)
    scores=np.array([np.median([r[np.argmin(abs(lags-d/v))] for r,d in zip(top,dist)]) for v in vv])
    return vv,scores

def run(args):
    db=pd.read_csv(CSV,sep=r"\s+"); db=db[db.nSamples>0].copy()
    db["t"]=pd.to_datetime(db.startTime,utc=True,errors="coerce")
    day=db[db.t.dt.strftime("%Y-%m-%d")==args.date].sort_values("t").iloc[args.start:args.start+args.nfiles]
    modes=args.modes.split(","); stacks={m:None for m in modes}; used=0; lags=dist=None
    for row in day.itertuples(index=False):
        f=corrected_path(row.file)
        if not f.exists(): continue
        x0,fs0,dx0=load_segment(f); x0=preprocess(x0,fs0,norm_seconds=5.0)
        for mode in modes:
            if mode=="unfiltered": x=x0[::2,::2]; fs=fs0/2; dx=2*dx0
            else: x,fs,dx=fk_filter(x0,fs0,dx0,mode)
            offset=int(round(50/dx)); targets=[int(round(50*j/dx)) for j in range(1,15) if int(round(50*j/dx))<x.shape[0]]
            l,cc=normalized_corr_pairs(x,[(0,c) for c in targets],fs)
            if stacks[mode] is None:
                stacks[mode]=cc; lags=l; dist=np.asarray(targets)*dx
            else: stacks[mode]+=cc
        used+=1
    results={}
    rng=np.random.default_rng(20260803)
    for mode in modes:
        top=stacks[mode]/used; vv,scores=velocity_scores(top,lags,dist)
        null=np.array([np.nanmax(velocity_scores(top,lags,rng.permutation(dist))[1]) for _ in range(500)])
        peak=int(np.nanargmax(scores)); results[mode]={"used":used,"peak_v_mps":float(vv[peak]),"peak_score":float(scores[peak]),"score_3200":float(scores[np.argmin(abs(vv-3200))]),"null95":float(np.quantile(null,.95)),"p_peak":float((1+np.sum(null>=scores[peak]))/(len(null)+1)),"top":top,"scores":scores,"vv":vv}
    OUT.mkdir(exist_ok=True); stem=f"fk_controls_{args.date}_start{args.start}_n{args.nfiles}"
    first=next(iter(results.values())); np.savez(OUT/(stem+".npz"),**{f"{m}_top":r["top"] for m,r in results.items()},**{f"{m}_scores":r["scores"] for m,r in results.items()},vv=first["vv"],lags=lags,dist=dist)
    clean={m:{k:v for k,v in r.items() if k not in ("top","scores","vv")} for m,r in results.items()}; (OUT/(stem+".json")).write_text(json.dumps(clean,indent=2))
    fig,axs=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for ax,(m,r) in zip(axs.flat,results.items()):
        ax.imshow(r["top"],extent=[lags[0],lags[-1],dist[-1],dist[0]],aspect="auto",cmap="RdBu_r",vmin=np.percentile(r["top"],2),vmax=np.percentile(r["top"],98)); ax.set_title(f"{m}: p={r['p_peak']:.3f}"); ax.set_xlabel("Lag (s)"); ax.set_ylabel("Separation (m)")
    fig.savefig(OUT/(stem+".png"),dpi=250); print(json.dumps(clean,indent=2))

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--date",default="2025-02-20"); ap.add_argument("--start",type=int,default=0); ap.add_argument("--nfiles",type=int,default=10); ap.add_argument("--modes",default="unfiltered,positive,negative,both,complement"); run(ap.parse_args())
