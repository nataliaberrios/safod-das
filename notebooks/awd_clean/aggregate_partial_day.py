#!/usr/bin/env python3
"""Combine completed checkpoint chunks into an interim partial-day stack."""
from pathlib import Path
import json
import numpy as np

OUT=Path(__file__).resolve().parent/"ambient_transfer"
DATE="2024-12-20"; NFILES=300; CHUNK=10; RNG=np.random.default_rng(20260804)

def score(top,lags,dist):
    vv=np.linspace(1200,6000,193)
    s=np.array([np.nanmedian([r[np.argmin(abs(lags-d/v))] for r,d in zip(top,dist)]) for v in vv])
    return vv,s

def main():
    zsum={m:None for m in ("negative","positive","both")}; usum=None; total=0; lags=dist=None
    for start in range(0,NFILES,CHUNK):
        z=np.load(OUT/f"fk_controls_{DATE}_start{start}_n10.npz"); meta=json.loads((OUT/f"fk_controls_{DATE}_start{start}_n10.json").read_text()); used=int(meta["negative"]["used"])
        if lags is None: lags,dist=z["lags"],z["dist"]
        for m in zsum: zsum[m]=z[f"{m}_top"]*used if zsum[m] is None else zsum[m]+z[f"{m}_top"]*used
        u=np.load(OUT/f"transfer_{DATE}_start{start}_n10.npz")
        usum=u["top_stack"]*used if usum is None else usum+u["top_stack"]*used
        total+=used
    results={"date":DATE,"files":total,"duration_hours":total/60,"norm_seconds":5.0,"passband_hz":[5.0,20.0],"methods":{}}
    for m,s in zsum.items():
        top=s/total; vv,sc=score(top,lags,dist); null=np.array([np.nanmax(score(top,lags,RNG.permutation(dist))[1]) for _ in range(2000)]); ip=int(np.nanargmax(sc)); i32=int(np.argmin(abs(vv-3200)))
        results["methods"]["fk_"+m]={"peak_v_mps":float(vv[ip]),"peak_score":float(sc[ip]),"score_3200":float(sc[i32]),"null95":float(np.quantile(null,.95)),"p_peak":float((1+np.sum(null>=sc[ip]))/(len(null)+1))}
        np.savez(OUT/f"interim_{DATE}_n{total}_fk_{m}.npz",top=top,lags=lags,dist=dist,vv=vv,scores=sc)
    top=usum/total; z=np.load(OUT/f"transfer_{DATE}_start0_n10.npz"); # exact path has same lag geometry
    lags_u=z["lags"]; dist_u=z["distances"]; vv,sc=score(top,lags_u,dist_u); null=np.array([np.nanmax(score(top,lags_u,RNG.permutation(dist_u))[1]) for _ in range(2000)]); ip=int(np.nanargmax(sc)); i32=int(np.argmin(abs(vv-3200)))
    results["methods"]["unfiltered"]={"peak_v_mps":float(vv[ip]),"peak_score":float(sc[ip]),"score_3200":float(sc[i32]),"null95":float(np.quantile(null,.95)),"p_peak":float((1+np.sum(null>=sc[ip]))/(len(null)+1))}
    np.savez(OUT/f"interim_{DATE}_n{total}_unfiltered.npz",top=top,lags=lags_u,dist=dist_u,vv=vv,scores=sc)
    (OUT/f"interim_{DATE}_n{total}.json").write_text(json.dumps(results,indent=2)); print(json.dumps(results,indent=2))

if __name__=="__main__": main()
