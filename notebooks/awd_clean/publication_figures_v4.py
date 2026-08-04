"""Publication-ready v4 figures for the SAFOD AWD-DAS dashboard.

The plots are intentionally separate from the diagnostic scripts: raw products
remain reproducible, while these outputs use the corrected 15-m source-offset
and vertical Nano-fiber geometry, cropped model coordinates, vector PDF output,
and a consistent journal-style visual system.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

HERE = Path(__file__).resolve().parent
MODEL = HERE / "vel_model"
X = np.array([-240,-6,-3,-1,0,.7,1.4,2,3,5,7,10,240.], float)
Y = np.array([-240,-8,-6,-4,-2,-1,0,1,2,4,6,8,240.], float)
Z = np.array([-150,-.5,0,.5,1,2,4,7,10,340.], float)
SOURCE = np.array([.015, 0., 0.])  # km; 15-m radial offset beside fiber
COLORS = {"blue":"#0072B2", "orange":"#D55E00", "teal":"#009E73", "gray":"#666666", "light":"#D9E6F2"}
plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":8.5,
    "axes.labelsize":8.5, "axes.titlesize":9, "xtick.labelsize":7.5,
    "ytick.labelsize":7.5, "legend.fontsize":7.2, "axes.linewidth":.7,
    "xtick.major.width":.6, "ytick.major.width":.6, "savefig.facecolor":"white"})

def panel(ax, letter):
    ax.text(-0.13, 1.04, letter, transform=ax.transAxes, fontweight="bold", fontsize=11, va="bottom")
    ax.tick_params(direction="out", length=3, width=.6)

def save(fig, stem):
    fig.savefig(HERE/(stem+".pdf"), bbox_inches="tight", pad_inches=.04)
    fig.savefig(HERE/(stem+".png"), dpi=400, bbox_inches="tight", pad_inches=.04)
    plt.close(fig)

def load_model():
    return np.loadtxt(MODEL/"Vp_model.dat").reshape(len(Z),len(Y),len(X))

def direct_paths(vp, depths):
    interp=RegularGridInterpolator((Z,Y,X),vp,bounds_error=False,fill_value=None)
    t=[]; v=[]
    for depth in depths:
        receiver=np.array([0.,0.,depth])
        q=np.linspace(SOURCE,receiver,401)
        speed=np.maximum(interp(q).astype(float),.25)
        ds=np.linalg.norm(np.diff(q,axis=0),axis=1)
        ti=np.sum(ds/((speed[:-1]+speed[1:])/2))
        t.append(ti); v.append(depth/ti)
    return np.asarray(t),np.asarray(v)

def wavelet(t, arrival, f0=40.):
    a=(np.pi*f0*(t-arrival))**2
    return (1-2*a)*np.exp(-a)

def forward():
    vp=load_model(); depths=np.linspace(.08,.50,71); tdir, vdir=direct_paths(vp,depths)
    guided=lambda f: 2.72+.004*(f-20.)
    time=np.arange(-.06,.24,.001)
    direct=np.array([wavelet(time,t) for t in tdir])
    tg=-.022+depths/(guided(40.)*1000.)
    guided_sec=np.array([wavelet(time,t) for t in tg])
    obs=pd.read_csv(HERE/"nano_mode_identification.csv")
    fig,ax=plt.subplots(2,2,figsize=(7.15,5.45),constrained_layout=True)
    # Crop padded model boundaries; show the actual SAFOD-scale region.
    xi=(X>=-3)&(X<=3); zi=(Z>=-.5)&(Z<=2)
    im=ax[0,0].pcolormesh(X[xi],Z[zi],vp[zi,6,:][:,xi],shading="auto",cmap="viridis",vmin=2,vmax=6)
    ax[0,0].plot([SOURCE[0],0],[SOURCE[2],depths[-1]],color="w",lw=1.2,ls="--")
    ax[0,0].set(xlabel="model X (km)",ylabel="model Z (km; positive down)",title="SAFOD prior: $V_P$ slice ($Y=0$)",xlim=(-3,3),ylim=(2,-.5)); panel(ax[0,0],"A")
    cb=fig.colorbar(im,ax=ax[0,0],pad=.02,fraction=.045); cb.set_label("$V_P$ (km s$^{-1}$)")
    extent=[time[0],time[-1],depths[0]*1000,depths[-1]*1000]
    for axis,sec,title,letter in [(ax[0,1],direct,"direct-P synthetic section","B"),(ax[1,0],guided_sec,"guided-mode sensitivity section","C")]:
        axis.imshow(sec,origin="lower",aspect="auto",extent=extent,cmap="RdBu_r",vmin=-1,vmax=1,interpolation="nearest")
        axis.set(xlabel="time relative to reference (s)",ylabel="fiber distance / depth (m)",title=title); panel(axis,letter)
    freq=obs.band_center_hz.to_numpy();
    ax[1,1].fill_between(freq,obs.bootstrap_speed_p2p5_mps,obs.bootstrap_speed_p97p5_mps,color=COLORS["blue"],alpha=.16,lw=0,label="95% burst-bootstrap interval")
    ax[1,1].plot(freq,obs.population_speed_mps,"o-",color=COLORS["blue"],lw=1.25,ms=4,label="observed Nano mode")
    ax[1,1].axhline(np.median(vdir)*1000,color=COLORS["orange"],lw=1.4,label=f"model direct-P median ({np.median(vdir)*1000:.0f} m s$^{{-1}}$)")
    ff=np.linspace(15,80,200); ax[1,1].plot(ff,guided(ff)*1000,"--",color=COLORS["gray"],lw=1.2,label="guided sensitivity curve")
    ax[1,1].set(xlabel="frequency (Hz)",ylabel="apparent speed (m s$^{-1}$)",title="observed mode and forward hypotheses",xlim=(12,82)); ax[1,1].grid(alpha=.18); ax[1,1].legend(frameon=False,loc="lower right"); panel(ax[1,1],"D")
    save(fig,"fig13_forward_publication")
    np.savez(HERE/"safod_forward_model_v4.npz",x_km=X,y_km=Y,z_km=Z,vp_km_s=vp,depths_km=depths,direct_arrival_s=tdir,direct_speed_km_s=vdir,guided_frequency_hz=ff,guided_speed_km_s=guided(ff))
    (HERE/"safod_forward_model_v4.txt").write_text(f"Corrected geometry: source offset 15 m; vertical Nano fiber; depths 80-500 m.\nMedian direct-P speed: {np.median(vdir)*1000:.1f} m/s; range {vdir.min()*1000:.1f}-{vdir.max()*1000:.1f} m/s.\nModel values are km/s on the Zhang-Thurber-Bedrosian grid. This is a forward sensitivity test, not an inversion.\n")

def sensitivity():
    d=np.load(HERE/"nano_dvv_injection_recovery.npz",allow_pickle=True); inj=d["injected_dvv"].astype(float); est=d["estimated_dvv"].astype(float); null=est[np.isclose(inj,0.)]; sig=float(np.std(null,ddof=1)); n=np.array([1,2,4,8,16,32,46]); null95=1.96*sig/np.sqrt(n); power95=3.24*sig/np.sqrt(n)
    c=pd.read_csv(HERE/"nano_stack_convergence.csv"); g=c.groupby("n_drops_per_substack")["independent_substack_ncc"]; ns=np.array(sorted(g.groups)); med=np.array([g.get_group(i).median() for i in ns]); lo=np.array([g.get_group(i).quantile(.16) for i in ns]); hi=np.array([g.get_group(i).quantile(.84) for i in ns])
    fig,ax=plt.subplots(1,2,figsize=(7.15,3.25),constrained_layout=True)
    ax[0].plot(n,null95,"o-",color=COLORS["blue"],lw=1.3,ms=4,label="95% null threshold"); ax[0].plot(n,power95,"s--",color=COLORS["orange"],lw=1.2,ms=3.5,label="approx. 95% power"); ax[0].axhline(5e-4,color=COLORS["gray"],ls=":",lw=1.2,label="illustrative $5\\times10^{-4}$ benchmark"); ax[0].set_xscale("log",base=2); ax[0].set_xticks(n); ax[0].set_xticklabels([str(i) for i in n]); ax[0].set(xlabel="independent stacks ($N$)",ylabel="$|\\Delta v_{app}/v_{app}|$",title="projected threshold"); ax[0].grid(alpha=.18); ax[0].legend(frameon=False,loc="upper right"); panel(ax[0],"A")
    ax[1].fill_between(ns,lo,hi,color=COLORS["blue"],alpha=.16,lw=0); ax[1].plot(ns,med,"o-",color=COLORS["blue"],lw=1.3,ms=4); ax[1].set_xscale("log",base=2); ax[1].set_xticks(ns); ax[1].set_xticklabels([str(i) for i in ns]); ax[1].set(xlabel="drops per independent substack",ylabel="NCC",title="measured repeatability convergence",ylim=(-.05,1.02)); ax[1].grid(alpha=.18); panel(ax[1],"B")
    save(fig,"fig14_monitoring_publication"); pd.DataFrame({"stack_size":n,"95pct_null_threshold":null95,"approx_95pct_power_threshold":power95}).to_csv(HERE/"monitoring_sensitivity_v4.csv",index=False)

def repeatability():
    d=pd.read_csv(HERE/"nano_drop_repeatability.csv"); b=pd.read_csv(HERE/"nano_burst_repeatability_hierarchical.csv"); c=pd.read_csv(HERE/"nano_stack_convergence.csv")
    fig,ax=plt.subplots(3,2,figsize=(7.15,8.1),constrained_layout=True); blue=COLORS["blue"]; orange=COLORS["orange"]; gray="#999999"
    x=d.drop_id.to_numpy(); burst=d.burst_id.to_numpy(); jitter=(x/np.maximum(1,d.n_drops_in_burst.to_numpy()-1)-.5)*.58
    ax[0,0].scatter(burst+jitter,d.loo_noise_ncc,s=3,color=gray,alpha=.25,label="noise window"); ax[0,0].scatter(burst+jitter,d.loo_signal_ncc,s=3,color=blue,alpha=.32,label="mode window"); ax[0,0].set(xlabel="burst index",ylabel="NCC",title="individual-drop repeatability",ylim=(-.5,1.02)); ax[0,0].legend(frameon=False,loc="lower right"); panel(ax[0,0],"A")
    ax[0,1].scatter(burst+jitter,1e3*d.loo_signal_delay_s,s=3,color=blue,alpha=.32); ax[0,1].axhline(0,color="k",lw=.6,ls="--"); ax[0,1].set(xlabel="burst index",ylabel="relative delay (ms)",title="drop timing stability"); panel(ax[0,1],"B")
    order=d.drop_id.to_numpy(); counts=d.groupby("drop_id").size(); keep=counts[counts>=10].index.to_numpy(); q=d[d.drop_id.isin(keep)].groupby("drop_id").relative_amplitude.agg(["median",lambda z:z.quantile(.16),lambda z:z.quantile(.84)]); ax[1,0].scatter(order,d.relative_amplitude,s=3,color=orange,alpha=.18); ax[1,0].plot(q.index,q["median"],color=orange,lw=1.3); ax[1,0].fill_between(q.index,q.iloc[:,1],q.iloc[:,2],color=orange,alpha=.16,lw=0); ax[1,0].axhline(1,color="k",lw=.6,ls="--"); ax[1,0].set(xlabel="drop order within burst",ylabel="relative RMS amplitude",title="amplitude stability"); panel(ax[1,0],"C")
    q=d[d.drop_id.isin(keep)].groupby("drop_id").beam_snr_db.agg(["median",lambda z:z.quantile(.16),lambda z:z.quantile(.84)]); ax[1,1].scatter(order,d.beam_snr_db,s=3,color=blue,alpha=.18); ax[1,1].plot(q.index,q["median"],color=blue,lw=1.3); ax[1,1].fill_between(q.index,q.iloc[:,1],q.iloc[:,2],color=blue,alpha=.16,lw=0); ax[1,1].axhline(0,color="k",lw=.6,ls="--"); ax[1,1].set(xlabel="drop order within burst",ylabel="beam SNR (dB)",title="individual-drop observability"); panel(ax[1,1],"D")
    ax[2,0].scatter(b.burst_id,b.loo_burst_signal_ncc,c=1e3*b.loo_burst_delay_s,cmap="coolwarm",vmin=-2,vmax=2,s=18,edgecolor="white",lw=.3); ax[2,0].set(xlabel="burst index",ylabel="leave-one-burst-out NCC",title="across-burst repeatability",ylim=(-.05,1.02)); panel(ax[2,0],"E")
    g=c.groupby("n_drops_per_substack").independent_substack_ncc; ns=np.array(sorted(g.groups)); med=np.array([g.get_group(i).median() for i in ns]); lo=np.array([g.get_group(i).quantile(.16) for i in ns]); hi=np.array([g.get_group(i).quantile(.84) for i in ns]); ax[2,1].fill_between(ns,lo,hi,color=blue,alpha=.16,lw=0); ax[2,1].plot(ns,med,"o-",color=blue,lw=1.3,ms=4); ax[2,1].set_xscale("log",base=2); ax[2,1].set_xticks(ns); ax[2,1].set_xticklabels([str(i) for i in ns]); ax[2,1].set(xlabel="drops per independent substack",ylabel="NCC",title="convergence with stacking",ylim=(-.05,1.02)); panel(ax[2,1],"F")
    save(fig,"fig12_repeatability_publication")

if __name__=="__main__":
    forward(); sensitivity(); repeatability()
