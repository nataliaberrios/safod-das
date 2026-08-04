"""Time-gated split-sample branch and intersection search for Deep AWD DAS."""

from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, sosfiltfilt

from fk_dispersion import weighted_stack


HERE = Path(__file__).resolve().parent
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
PRE_S = 0.5
TURNAROUND_CH = 1702
TIME_WINDOW = (-0.08, 2.85)
BANDS = ((3.0, 15.0), (15.0, 30.0))
APERTURE_M = 400.0
STEP_M = 200.0
P_POS = np.linspace(1 / 1800.0, 1 / 1300.0, 13)
P_GRID = np.r_[-P_POS[::-1], P_POS]
TIME_DECIMATION = 4
CHANNEL_STRIDE = 4
WINDOW_S = 0.12
N_PERM = 199
SEED = 20260802


def prepare(section, fs, band):
    sos = butter(4, band, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, np.nan_to_num(section), axis=-1)[:, ::TIME_DECIMATION]


def aligned(section, coordinate, fs, p):
    """Align a local linear trajectory to its arrival time at coordinate mean."""
    center = coordinate.mean()
    shifts = np.rint(p * (coordinate - center) * fs).astype(int)
    out = np.full_like(section, np.nan)
    for i, shift in enumerate(shifts):
        if shift > 0:
            out[i, :-shift] = section[i, shift:]
        elif shift < 0:
            out[i, -shift:] = section[i, :shift]
        else:
            out[i] = section[i]
    return out


def semblance_series(section, coordinate, fs, p):
    gather = aligned(section, coordinate, fs, p)
    valid = np.isfinite(gather)
    values = np.nan_to_num(gather)
    numerator = np.sum(values, axis=0) ** 2
    denominator = valid.sum(axis=0) * np.sum(values ** 2, axis=0)
    length = max(3, int(round(WINDOW_S * fs)))
    numerator = uniform_filter1d(numerator, length, mode="constant")
    denominator = uniform_filter1d(denominator, length, mode="constant")
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)


def scan(section, coordinate, fs, time):
    allowed = (time >= 0.0) & (time <= 2.75)
    best = {"positive": (-np.inf, np.nan, np.nan), "negative": (-np.inf, np.nan, np.nan)}
    for p in P_GRID:
        series = semblance_series(section, coordinate, fs, p)
        index = np.flatnonzero(allowed)[np.argmax(series[allowed])]
        key = "positive" if p > 0 else "negative"
        if series[index] > best[key][0]:
            best[key] = (float(series[index]), float(p), float(time[index]))
    return best


def fixed_validation(section, coordinate, fs, time, p, target_time):
    series = semblance_series(section, coordinate, fs, p)
    allowed = np.abs(time - target_time) <= 0.06
    return float(np.max(series[allowed])) if np.any(allowed) else np.nan


def trajectory_score(section, coordinate, fs, time, p, target_time):
    series = semblance_series(section, coordinate, fs, p)
    index = int(np.argmin(np.abs(time - target_time)))
    half = max(1, int(round(0.03 * fs)))
    return float(np.max(series[max(0, index-half):index+half+1]))


def main():
    rng = np.random.default_rng(SEED)
    with np.load(STACKS) as data:
        stacks = data["deep_stacks"]
        counts = data["n_common"]
        fs0 = float(data["fs"]); dx = float(data["dx_deep"])
        epochs = np.flatnonzero(counts > 0)
        odd = counts.copy(); odd[epochs[epochs % 2 == 0]] = 0
        even = counts.copy(); even[epochs[epochs % 2 == 1]] = 0
        discovery = weighted_stack(stacks, odd)
        validation = weighted_stack(stacks, even)

    i0 = int(round((PRE_S + TIME_WINDOW[0]) * fs0))
    i1 = int(round((PRE_S + TIME_WINDOW[1]) * fs0))
    fs = fs0 / TIME_DECIMATION
    time = np.arange(0, i1-i0, TIME_DECIMATION) / fs0 + TIME_WINDOW[0]
    legs = {
        "outbound": (discovery[:TURNAROUND_CH, i0:i1], validation[:TURNAROUND_CH, i0:i1]),
        "return": (discovery[TURNAROUND_CH:, i0:i1][::-1], validation[TURNAROUND_CH:, i0:i1][::-1]),
    }
    width = int(round(APERTURE_M / dx)); step = int(round(STEP_M / dx))
    rows = []; intersections = []
    products = {"time": time, "p_grid": P_GRID, "bands": np.asarray(BANDS)}

    for leg, (disc_leg, val_leg) in legs.items():
        starts = np.arange(0, disc_leg.shape[0]-width+1, step)
        centers = (starts + 0.5*(width-1))*dx
        products[f"{leg}_centers_m"] = centers
        for ib, band in enumerate(BANDS):
            disc_f = prepare(disc_leg, fs0, band)
            val_f = prepare(val_leg, fs0, band)
            for ia, start in enumerate(starts):
                channels = np.arange(start, start+width, CHANNEL_STRIDE)
                coordinate = channels * dx
                dsec = disc_f[channels]; vsec = val_f[channels]
                result = scan(dsec, coordinate, fs, time)
                stored = {}
                for direction in ("positive", "negative"):
                    dsem, p, tc = result[direction]
                    vsem = fixed_validation(vsec, coordinate, fs, time, p, tc)
                    row = {
                        "leg": leg, "band_low_hz": band[0], "band_high_hz": band[1],
                        "aperture_index": ia, "start_m": coordinate[0], "stop_m": coordinate[-1],
                        "center_m": centers[ia], "direction": direction,
                        "discovery_slowness_ms_per_m": p*1e3,
                        "discovery_center_time_s": tc,
                        "discovery_semblance": dsem, "validation_semblance": vsem,
                    }
                    rows.append(row); stored[direction] = row
                pos, neg = stored["positive"], stored["negative"]
                pp=pos["discovery_slowness_ms_per_m"]*1e-3
                pn=neg["discovery_slowness_ms_per_m"]*1e-3
                xc=centers[ia]
                xi=xc+(neg["discovery_center_time_s"]-pos["discovery_center_time_s"])/(pp-pn)
                ti=pos["discovery_center_time_s"]+pp*(xi-xc)
                if coordinate[0] <= xi <= coordinate[-1] and 0 <= ti <= 2.75:
                    observed=min(pos["validation_semblance"],neg["validation_semblance"])
                    null=np.empty(N_PERM)
                    for k in range(N_PERM):
                        perm=rng.permutation(vsec.shape[0])
                        sp=trajectory_score(vsec[perm],coordinate,fs,time,pp,pos["discovery_center_time_s"])
                        sn=trajectory_score(vsec[perm],coordinate,fs,time,pn,neg["discovery_center_time_s"])
                        null[k]=min(sp,sn)
                    probability=(1+np.sum(null>=observed))/(N_PERM+1)
                    intersections.append({
                        "leg":leg,"band_low_hz":band[0],"band_high_hz":band[1],
                        "aperture_index":ia,"intersection_fiber_m":xi,
                        "intersection_time_from_p26_UTC_Date_s":ti,
                        "positive_validation_semblance":pos["validation_semblance"],
                        "negative_validation_semblance":neg["validation_semblance"],
                        "joint_validation_statistic":observed,"permutation_p":probability,
                    })

    with (HERE/"deep_time_gated_branches.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with (HERE/"deep_candidate_intersections.csv").open("w",newline="") as f:
        if intersections:
            w=csv.DictWriter(f,fieldnames=list(intersections[0])); w.writeheader(); w.writerows(intersections)

    fig,axes=plt.subplots(2,2,figsize=(13,9),sharex=False,sharey=True)
    for row_i,leg in enumerate(("outbound","return")):
        for col_i,band in enumerate(BANDS):
            ax=axes[row_i,col_i]
            subset=[r for r in rows if r["leg"]==leg and r["band_low_hz"]==band[0]]
            for direction,marker,cmap in (("positive","o","viridis"),("negative","v","plasma")):
                q=[r for r in subset if r["direction"]==direction]
                sc=ax.scatter([r["discovery_center_time_s"] for r in q],[r["center_m"] for r in q],
                              c=[r["validation_semblance"] for r in q],marker=marker,s=55,
                              cmap=cmap,vmin=0,vmax=1,label=f"{direction} slowness")
            cand=[r for r in intersections if r["leg"]==leg and r["band_low_hz"]==band[0] and r["permutation_p"]<=0.05]
            ax.scatter([r["intersection_time_from_p26_UTC_Date_s"] for r in cand],
                       [r["intersection_fiber_m"] for r in cand],marker="*",s=130,
                       color="red",edgecolor="white",label="intersection p≤0.05")
            ax.invert_yaxis(); ax.grid(alpha=.2)
            ax.set_title(f"{leg.capitalize()}, {band[0]:.0f}–{band[1]:.0f} Hz")
            ax.set_xlabel("time relative to p26.cc9.txt UTC_Date (s)")
            ax.set_ylabel("distance from provisional leg start (m)")
            if row_i==0 and col_i==0: ax.legend(frameon=False,fontsize=8)
    fig.suptitle("Deep time-gated signed branches and candidate intersections\nodd-burst discovery; even-burst validation color")
    fig.tight_layout(); fig.savefig(HERE/"deep_time_gated_branches.png",dpi=200)
    np.savez(HERE/"deep_time_gated_branches.npz",**products)
    significant=[r for r in intersections if r["permutation_p"]<=0.05]
    report=(f"Time reference: p26.cc9.txt UTC_Date; physical meaning not independently documented.\n"
            f"Branch rows: {len(rows)}; geometric intersections: {len(intersections)}; "
            f"permutation p<=0.05: {len(significant)}.\n"
            "Intersection positions are provisional fiber coordinates, not true depth or accepted fracture locations.\n")
    (HERE/"deep_time_gated_branches.txt").write_text(report)
    print(report,end="")


if __name__ == "__main__":
    main()
