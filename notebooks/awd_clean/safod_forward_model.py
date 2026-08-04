"""Analytical SAFOD-model forward test for Nano mode identification.

This is deliberately a transparent sensitivity calculation, not a Vp inversion
or a full elastic-wave simulation.  It reads the Zhang--Thurber--Bedrosian
SAFOD gridded model, integrates slowness along representative source--receiver
paths, and compares that direct-P prediction with a weakly dispersive guided
mode.  The synthetic DAS panels use the same along-fiber moveout observable
as the dashboard and are intended to answer whether the existing data can
distinguish the two hypotheses.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

HERE = Path(__file__).resolve().parent
MODEL = HERE / "vel_model"
OUT = HERE
X = np.array([-240, -6, -3, -1, 0, .7, 1.4, 2, 3, 5, 7, 10, 240.], float)
Y = np.array([-240, -8, -6, -4, -2, -1, 0, 1, 2, 4, 6, 8, 240.], float)
Z = np.array([-150, -.5, 0, .5, 1, 2, 4, 7, 10, 340.], float)
SOURCE = np.array([-.20, 0.0, .50])  # km; representative shallow path
RECEIVER_OFFSET_KM = np.linspace(.08, .50, 71)
OBSERVED = HERE / "nano_mode_identification.csv"


def read_model(name):
    values = np.loadtxt(MODEL / name).reshape(len(Z), len(Y), len(X))
    return values


def integrated_direct_arrival(vp, offsets):
    """Integrate 1/Vp along straight shallow source--receiver paths."""
    interp = RegularGridInterpolator((Z, Y, X), vp, bounds_error=False,
                                     fill_value=None)
    times = []
    mean_v = []
    for offset in offsets:
        receiver = SOURCE + np.array([offset, 0., 0.])
        q = np.linspace(SOURCE, receiver, 401)
        speed = interp(q).astype(float)
        speed = np.maximum(speed, .25)  # model units are km/s
        distance = np.linalg.norm(np.diff(q, axis=0), axis=1)
        times.append(np.sum(distance / ((speed[:-1] + speed[1:]) / 2)))
        mean_v.append(offset / times[-1])
    return np.asarray(times), np.asarray(mean_v)


def ricker(t, f0):
    a = (np.pi * f0 * (t - .045)) ** 2
    return (1 - 2 * a) * np.exp(-a)


def synthetic_section(offsets_m, times, c_of_f=None, f0=40.):
    """Return a simple band-limited synthetic DAS moveout section."""
    q = np.arange(-.06, .24, .001)
    section = np.zeros((len(offsets_m), len(q)))
    if c_of_f is None:
        arrivals = times
    else:
        arrivals = -.022 + offsets_m / (c_of_f(f0) * 1000.)
    for i, arrival in enumerate(arrivals):
        section[i] = ricker(q - arrival + .045, f0)
    return q, section


def main():
    vp = read_model("Vp_model.dat")
    vs = read_model("Vs_model.dat")
    vpvs = read_model("Vpvs_model.dat")
    offsets = RECEIVER_OFFSET_KM
    direct_t, direct_v = integrated_direct_arrival(vp, offsets)
    direct_c = lambda f: np.interp(f, [15, 30, 45, 60, 80],
                                   [np.median(direct_v), np.median(direct_v),
                                    np.median(direct_v), np.median(direct_v),
                                    np.median(direct_v)])
    # A deliberately weakly dispersive alternative, anchored to the observed
    # Nano ridge rather than asserted as a physical tube-wave law.
    guided_c = lambda f: 2.72 + .004 * (f - 20.)
    t, direct_section = synthetic_section(offsets * 1000., direct_t, f0=40.)
    _, guided_section = synthetic_section(offsets * 1000., direct_t,
                                          c_of_f=guided_c, f0=40.)
    observed = pd.read_csv(OBSERVED)
    model_speed = np.median(direct_v) * 1000.
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    im = ax[0, 0].imshow(vp[:, :, 6], origin="upper", aspect="auto",
                         extent=[Y[0], Y[-1], Z[-1], Z[0]], cmap="viridis")
    ax[0, 0].axhline(SOURCE[2], color="w", ls="--", lw=1)
    ax[0, 0].set(xlabel="model Y (km)", ylabel="model Z (km)",
                 title="A  SAFOD prior: Vp slice at X=0")
    fig.colorbar(im, ax=ax[0, 0], label="Vp (km s$^{-1}$)")
    extent = [t[0], t[-1], offsets[0] * 1000, offsets[-1] * 1000]
    for axis, section, title in [(ax[0, 1], direct_section,
                                  "B  Direct-P synthetic moveout"),
                                 (ax[1, 0], guided_section,
                                  "C  Weakly dispersive guided-mode alternative")]:
        axis.imshow(section, origin="lower", aspect="auto", extent=extent,
                    cmap="RdBu_r", vmin=-1, vmax=1)
        axis.set(xlabel="relative time (s)", ylabel="along-fiber offset (m)",
                 title=title)
    freq = observed.band_center_hz.to_numpy()
    ax[1, 1].plot(freq, observed.population_speed_mps, "o-", color="#0b5d8b",
                   label="observed Nano mode")
    ax[1, 1].axhline(model_speed, color="#d95f02", lw=2,
                     label=f"SAFOD-model direct-P median ({model_speed:.0f} m/s)")
    ff = np.linspace(15, 80, 200)
    ax[1, 1].plot(ff, guided_c(ff) * 1000, "--", color="#555555",
                   label="guided-mode sensitivity curve")
    ax[1, 1].fill_between(freq, observed.bootstrap_speed_p2p5_mps,
                          observed.bootstrap_speed_p97p5_mps,
                          color="#0b5d8b", alpha=.15, label="bootstrap interval")
    ax[1, 1].set(xlabel="frequency (Hz)", ylabel="apparent speed (m s$^{-1}$)",
                 title="D  Observed versus forward hypotheses")
    ax[1, 1].legend(frameon=False, fontsize=8)
    fig.suptitle("SAFOD AWD Nano forward sensitivity test — not a Vp inversion",
                 fontsize=13)
    fig.savefig(OUT / "safod_forward_model.png", dpi=220)
    np.savez(OUT / "safod_forward_model.npz", x_km=X, y_km=Y, z_km=Z,
             vp_km_s=vp, vs_km_s=vs, vpvs=vpvs, offsets_km=offsets,
             direct_arrival_s=direct_t, direct_speed_km_s=direct_v,
             guided_frequency_hz=ff, guided_speed_km_s=guided_c(ff))
    (OUT / "safod_forward_model.txt").write_text(
        "SAFOD AWD Nano forward sensitivity test\n"
        "========================================\n"
        "Model: Vp_model.dat, Vs_model.dat, Vpvs_model.dat; coordinates parsed from MOD.head.\n"
        f"Grid: {len(X)} x {len(Y)} x {len(Z)} padded nodes; units km and km/s.\n"
        f"Representative straight path: source {SOURCE.tolist()} km, offsets {offsets[0]:.2f}-{offsets[-1]:.2f} km.\n"
        f"Median integrated direct-P speed: {model_speed:.1f} m/s; range {direct_v.min()*1000:.1f}-{direct_v.max()*1000:.1f} m/s.\n"
        "Guided alternative: c(f)=2.72+0.004(f-20) km/s, a sensitivity hypothesis anchored near the observed ridge.\n"
        "Interpretation: this calculation tests identifiability and geometry sensitivity; it does not invert Vp or prove a phase.\n"
    )
    print(f"Saved {OUT/'safod_forward_model.png'}; median direct-P {model_speed:.1f} m/s")


if __name__ == "__main__":
    main()
