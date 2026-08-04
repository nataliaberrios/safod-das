"""Synthetic sign check for the signed F-K masks used in the ambient branch.

The test uses the same array orientation and decimation as
ambient_fk_transfer_test.py.  A wave cos(2*pi*f*(t-z/v)) propagates toward
increasing fiber coordinate z; cos(2*pi*f*(t+z/v)) propagates toward
decreasing z.  The output identifies which F*K mask retains each wave.  This
does not assign physical upgoing/downgoing labels until fiber orientation is
known.
"""
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "fk_sign_synthetic_test.json"
OUT_FIG = HERE / "fk_sign_synthetic_test.png"


def signed_masks(nz, nt, fs, dx, fmin=5.0, fmax=20.0, vmin=2500.0, vmax=4500.0):
    f = np.fft.fftfreq(nt, 1.0 / fs)
    k = np.fft.fftfreq(nz, dx)
    K, F = np.meshgrid(k, f, indexing="ij")
    af, ak = np.abs(F), np.abs(K)
    v = af / np.maximum(ak, 1e-12)
    base = (af >= fmin) & (af <= fmax) & (v >= vmin) & (v <= vmax) & (ak > 0)
    return base & (F * K < 0), base & (F * K > 0), F, K


def apply_branch(x, mask):
    return np.fft.ifft2(np.fft.fft2(x) * mask).real


def main():
    # Use the same physical sampling convention as the current Nano archive.
    fs0, dx0 = 500.0, 1.020952
    nz0, nt0 = 512, 8192
    z = np.arange(nz0)[:, None] * dx0
    t = np.arange(nt0)[None, :] / fs0
    f0, v0 = 10.0, 3200.0
    # Avoid a zero-frequency/edge ambiguity while retaining the exact phase sign.
    plus = np.cos(2.0 * np.pi * f0 * (t - z / v0))
    minus = np.cos(2.0 * np.pi * f0 * (t + z / v0))

    # Exact decimation used by fk_filter in ambient_fk_transfer_test.py.
    fs, dx = fs0 / 2.0, dx0 * 2.0
    plus, minus = plus[::2, ::2], minus[::2, ::2]
    neg, pos, F, K = signed_masks(plus.shape[0], plus.shape[1], fs, dx)
    results = {}
    for name, wave in (("increasing_z", plus), ("decreasing_z", minus)):
        spec_power = np.abs(np.fft.fft2(wave)) ** 2
        total = float(np.sum(spec_power))
        neg_power = float(np.sum(spec_power[neg]))
        pos_power = float(np.sum(spec_power[pos]))
        neg_out = apply_branch(wave, neg)
        pos_out = apply_branch(wave, pos)
        results[name] = {
            "negative_F_times_K_power_fraction": neg_power / total,
            "positive_F_times_K_power_fraction": pos_power / total,
            "negative_output_rms": float(np.sqrt(np.mean(neg_out ** 2))),
            "positive_output_rms": float(np.sqrt(np.mean(pos_out ** 2))),
            "retained_branch": "negative (F*K<0)" if neg_power > pos_power else "positive (F*K>0)",
        }

    summary = {
        "synthetic_wave_plus": "cos(2*pi*f*(t-z/v)): propagation toward increasing fiber coordinate z",
        "synthetic_wave_minus": "cos(2*pi*f*(t+z/v)): propagation toward decreasing fiber coordinate z",
        "f_hz": f0,
        "velocity_m_s": v0,
        "original_fs_hz": fs0,
        "original_spacing_m": dx0,
        "post_decimation_fs_hz": fs,
        "post_decimation_spacing_m": dx,
        "mask": "5-20 Hz, 2.5-4.5 km/s, F*K sign, same decimation as production code",
        "results": results,
        "physical_label_boundary": "F*K branch labels are coordinate-direction labels only; upgoing/downgoing requires fiber orientation and time/depth convention.",
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    names = ["increasing_z", "decreasing_z"]
    neg_rms = [results[n]["negative_output_rms"] for n in names]
    pos_rms = [results[n]["positive_output_rms"] for n in names]
    x = np.arange(2)
    axes[0].bar(x - 0.18, neg_rms, 0.36, label="F*K < 0")
    axes[0].bar(x + 0.18, pos_rms, 0.36, label="F*K > 0")
    axes[0].set_xticks(x, ["increasing z", "decreasing z"])
    axes[0].set_ylabel("Filtered output RMS")
    axes[0].set_title("Exact production-mask sign test")
    axes[0].legend(frameon=False)
    axes[1].bar(x - 0.18, [results[n]["negative_F_times_K_power_fraction"] for n in names], 0.36, label="F*K < 0")
    axes[1].bar(x + 0.18, [results[n]["positive_F_times_K_power_fraction"] for n in names], 0.36, label="F*K > 0")
    axes[1].set_xticks(x, ["increasing z", "decreasing z"])
    axes[1].set_ylabel("Input spectral power fraction")
    axes[1].set_title("Power retained by signed wedges")
    axes[1].legend(frameon=False)
    fig.suptitle("Synthetic sign convention check: 10 Hz, 3.2 km/s")
    fig.savefig(OUT_FIG, dpi=250)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
