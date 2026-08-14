import itertools, glob, re, sys, os
import numpy as np
sys.path.insert(0,".")
from ambient_lellouch2019_exact_stack import (VELOCITY_GRID_M_S, correlation_from_spectrum,
    moveout_scores, prepare_aggregate_spectra, receiver_order_null)
DIRS=["ambient_transfer/lellouch2019_exact_stack","ambient_transfer/lellouch2019_exact_stack_days"]
PAT="chunk_*_src23_ram0p1_cross_correlation_ordered_r0_start*_n0060.npz"
by={}
for d in DIRS:
    for f in sorted(glob.glob(os.path.join(d,PAT))):
        if "_cm_" in f: continue
        dt=re.search(r"chunk_(\d{4}-\d{2}-\d{2})_",f).group(1); by.setdefault(dt,[]).append(f)
ok=[]
for d in sorted(by):
    st=sorted(int(re.search(r"start(\d+)",f).group(1)) for f in by[d])
    if st==list(range(0,60*len(st),60)):
        with np.load(by[d][0]) as z:
            if float(z["fs_hz"])==500.0: ok.append(d)
print("days at 500 Hz, contiguous:",ok)
def stack(ds):
    c=n=None; nw=0; ref=None
    for d in ds:
        for f in sorted(by[d]):
            with np.load(f) as z:
                if c is None:
                    c=np.zeros_like(z["central_cross_spectrum_sum"]); n=np.zeros_like(z["neighbor_cross_spectrum_sum"])
                    ref=(float(z["fs_hz"]),int(z["n_fft"]),np.asarray(z["offsets_m"],float),str(z["spectral_mode"]))
                c+=z["central_cross_spectrum_sum"]; n+=z["neighbor_cross_spectrum_sum"]; nw+=int(z["n_windows"])
    fs,nfft,offs,mode=ref
    cs_,ns_,_,_,_=prepare_aggregate_spectra(c,n,None,nw,nfft,fs,mode,1e-3)
    lags,_=correlation_from_spectrum(cs_,1,nfft,fs,apply_bandpass=True)
    _,sec=correlation_from_spectrum(ns_,1,nfft,fs,apply_bandpass=True)
    sc=moveout_scores(sec,lags,offs,sign=1.0); ac=moveout_scores(sec,lags,offs,sign=-1.0)
    k=int(np.argmax(sc)); pk=float(sc[k])
    _,p=receiver_order_null(sec,lags,offs,pk,20260814,2000)
    return pk,float(VELOCITY_GRID_M_S[k]),p,float(np.interp(3200,VELOCITY_GRID_M_S,sc)),float(np.interp(3200,VELOCITY_GRID_M_S,ac)),nw
print("\nALL %d pairs of the %d usable days (2000-perm null):"%(len(list(itertools.combinations(ok,2))),len(ok)))
print("%-26s %8s %7s %8s %7s %s"%("pair","peak","@v","p","c/a@3200","<0.05?"))
res=[]
for a,b in itertools.combinations(ok,2):
    pk,v,p,c32,a32,nw=stack([a,b]); res.append(p)
    print("%-26s %8.3f %7.0f %8.4f %7.2f  %s"%(a[5:]+"+"+b[5:],pk,v,p,c32/a32,"YES" if p<0.05 else ""))
res=np.array(res)
print("\npairs with p<0.05: %d of %d   min p %.4f"%((res<0.05).sum(),len(res),res.min()))
print("expected by chance at alpha=0.05 over %d pairs: %.1f"%(len(res),0.05*len(res)))
