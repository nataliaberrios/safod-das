"""Audit p26.cc9 arrival-referenced times against PEG500 and Nano waveforms."""

import os
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from obspy import read, UTCDateTime
from scipy.signal import butter, hilbert, sosfiltfilt


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CATALOG = ROOT / "p26.cc9.txt"
NODE_DIR = Path("/oak/stanford/groups/ettore88/data/SAFOD/ActiveJune2026/Nodes/PEG500_nodes")
STACKS = HERE / "canonical_epoch_stacks_paired_deep_all.npz"
PRE, POST = 0.25, 0.50
NODE_BAND = (15.0, 80.0)
NANO_BAND = (30.0, 60.0)


def bp(x, fs, band):
    sos = butter(4, band, btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def catalog_rows():
    rows = []
    with CATALOG.open() as stream:
        next(stream)
        for line in stream:
            fields = line.split()
            if len(fields) >= 3:
                rows.append((fields[0], UTCDateTime(fields[1]), float(fields[2])))
    return rows


def node_windows(rows):
    grouped = defaultdict(list)
    for filename, time, cc in rows:
        grouped[filename].append((time, cc))
    windows, quality = [], []
    for filename, picks in grouped.items():
        path = NODE_DIR / filename
        if not path.exists():
            continue
        trace = read(str(path))[0]
        fs = trace.stats.sampling_rate
        n = int(round((PRE + POST) * fs))
        for time, cc in picks:
            start = int(round((time - trace.stats.starttime - PRE) * fs))
            if start >= 0 and start + n <= trace.stats.npts:
                windows.append(trace.data[start:start + n].astype(float))
                quality.append(cc)
    if not windows:
        raise RuntimeError("No reference-node windows could be extracted")
    return np.asarray(windows), np.asarray(quality), fs


def robust_onset(trace, time, pre=(-0.20, -0.05), search=(-0.05, 0.15)):
    env = np.abs(hilbert(trace))
    base = env[(time >= pre[0]) & (time <= pre[1])]
    threshold = np.median(base) + 8 * np.median(np.abs(base - np.median(base)))
    candidates = np.where((time >= search[0]) & (time <= search[1]) & (env > threshold))[0]
    return (time[candidates[0]] if candidates.size else np.nan), threshold, env


def main():
    rows = catalog_rows()
    node, cc, fs_node = node_windows(rows)
    node = bp(node, fs_node, NODE_BAND)
    # Weight high-quality catalog picks, but report an unweighted robustness check.
    weights = np.clip(cc, 0, 1) ** 2
    node_stack = np.average(node, axis=0, weights=weights)
    node_time = np.arange(node.shape[1]) / fs_node - PRE
    onset, threshold, node_env = robust_onset(node_stack, node_time)
    peak_window = (node_time >= -0.05) & (node_time <= 0.15)
    node_peak = node_time[peak_window][np.argmax(np.abs(node_stack[peak_window]))]

    product = np.load(STACKS)
    counts = product["n_common"]
    good = counts > 0
    nano = np.tensordot(
        counts[good].astype(float), product["nano_stacks"][good], axes=(0, 0)
    ) / counts[good].sum()
    fs_nano = float(product["fs"])
    dx = float(product["dx_nano"])
    nano = bp(nano, fs_nano, NANO_BAND)
    nano_time = np.arange(nano.shape[1]) / fs_nano - 0.5

    # Earliest robust energy by channel, without imposing a velocity trajectory.
    channels = np.arange(int(20 / dx), int(440 / dx) + 1, 4)
    nano_onsets = []
    for channel in channels:
        value, _, _ = robust_onset(
            nano[channel], nano_time, pre=(-0.20, -0.08), search=(-0.08, 0.25)
        )
        nano_onsets.append(value)
    nano_onsets = np.asarray(nano_onsets)
    distance = channels * dx

    fig, ax = plt.subplots(1, 3, figsize=(16, 5.5))
    # Reference seismometer: individual gray windows plus stack.
    scale = np.max(np.abs(node_stack)) + 1e-30
    for waveform in node[::max(1, len(node)//80)]:
        ax[0].plot(node_time, waveform / scale, color="0.75", lw=0.35, alpha=0.25)
    ax[0].plot(node_time, node_stack / scale, color="k", lw=1.5, label="weighted stack")
    ax[0].plot(node_time, node_env / scale, color="C1", lw=1, label="envelope")
    ax[0].axvline(0, color="crimson", ls="--", label="catalog UTC_Date")
    ax[0].axvline(onset, color="C2", ls=":", label=f"onset {onset*1e3:.1f} ms")
    ax[0].axvline(node_peak, color="C0", ls=":", label=f"peak {node_peak*1e3:.1f} ms")
    ax[0].set_xlim(-0.10, 0.20)
    ax[0].set_xlabel("time relative to catalog pick (s)")
    ax[0].set_ylabel("normalized reference-node amplitude")
    ax[0].set_title(f"A  PEG500 reference, {NODE_BAND[0]:.0f}-{NODE_BAND[1]:.0f} Hz\n{len(node)} drops")
    ax[0].legend(fontsize=7)

    clip = np.percentile(np.abs(nano[channels]), 99.5)
    ax[1].imshow(
        nano[channels], aspect="auto", origin="upper", cmap="RdBu_r",
        extent=[nano_time[0], nano_time[-1], distance[-1], distance[0]],
        vmin=-clip, vmax=clip,
    )
    ax[1].plot(nano_onsets, distance, "k.", ms=3, label="threshold onset")
    ax[1].axvline(0, color="gold", lw=1, label="catalog pick")
    ax[1].set_xlim(-0.08, 0.25)
    ax[1].set_xlabel("time relative to catalog pick (s)")
    ax[1].set_ylabel("distance along Nano fiber (m)")
    ax[1].set_title(f"B  Nano stack, {NANO_BAND[0]:.0f}-{NANO_BAND[1]:.0f} Hz")
    ax[1].legend(fontsize=7)

    valid = np.isfinite(nano_onsets)
    ax[2].plot(nano_onsets[valid] * 1e3, distance[valid], "k.-", lw=0.7, ms=4)
    ax[2].axvline(0, color="crimson", ls="--")
    ax[2].invert_yaxis()
    ax[2].set_xlabel("threshold-onset time relative to catalog (ms)")
    ax[2].set_ylabel("distance along Nano fiber (m)")
    ax[2].set_title("C  Unconstrained Nano onset picks\n(QC diagnostic, not a velocity model)")
    ax[2].grid(alpha=0.3)
    fig.suptitle("SAFOD AWD source-time audit: catalog times are seismometer CC picks")
    fig.tight_layout()
    fig.savefig(HERE / "source_time_audit.png", dpi=180)

    report = (
        "SAFOD AWD source-time audit\n"
        "Timing definition: p26.cc9 UTC_Date is a cross-correlation pick on nearby "
        "PEG500 node 453009664, not a physical impact trigger.\n"
        f"catalog rows: {len(rows)}\n"
        f"reference-node windows extracted: {len(node)}\n"
        f"reference-node sampling rate: {fs_node:.1f} Hz\n"
        f"reference-node analysis band: {NODE_BAND[0]:.1f}-{NODE_BAND[1]:.1f} Hz\n"
        f"stack envelope onset relative to UTC_Date: {onset*1e3:.3f} ms\n"
        f"stack absolute-amplitude peak relative to UTC_Date: {node_peak*1e3:.3f} ms\n"
        "Consequence: absolute DAS travel times require the reference node's source "
        "distance/coordinates and a physical impact-time calibration. Relative moveout "
        "and interval slowness do not require that constant correction.\n"
    )
    (HERE / "source_time_audit.txt").write_text(report)
    np.savez(
        HERE / "source_time_audit.npz", node_time=node_time,
        node_stack=node_stack, node_envelope=node_env, node_onset=onset,
        node_peak=node_peak, nano_distance=distance, nano_onsets=nano_onsets,
        catalog_cc=cc,
    )
    print(report)
    print("Saved source_time_audit.png/.npz/.txt")


if __name__ == "__main__":
    main()
