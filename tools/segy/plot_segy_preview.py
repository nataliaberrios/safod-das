#!/usr/bin/env python3
"""
Preview an OptaSense/PASSCAL-style SEGY file from SAFOD DAS.

This is intentionally conservative: it reads one file, a limited trace range,
and saves quick-look figures plus a header text file. It uses DASutils when
available and falls back to a minimal IEEE-float SEGY reader for previewing.
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path
from typing import BinaryIO

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


COMMON_DASUTILS_PATHS = (
    "/home/groups/ettore88/nberrios/safod_das_git/DAS-utilities/python",
    "/home/groups/ettore88/nberrios/safod-das/DAS-utilities/python",
)


def add_dasutils_paths(extra_path: str | None) -> None:
    paths = []
    if extra_path:
        paths.append(extra_path)
    env_path = os.environ.get("DASUTILS_PATH")
    if env_path:
        paths.append(env_path)
    paths.extend(COMMON_DASUTILS_PATHS)

    for path in paths:
        if path and Path(path).exists() and path not in sys.path:
            sys.path.insert(0, path)


def import_dasutils(extra_path: str | None):
    add_dasutils_paths(extra_path)
    try:
        import DASutils  # type: ignore

        return DASutils
    except Exception as exc:  # pragma: no cover - reported at runtime
        print(f"Warning: could not import DASutils ({exc}); using fallback reader.")
        return None


def open_binary(path: Path) -> BinaryIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def file_size_uncompressed(path: Path) -> int | None:
    if path.suffix.lower() == ".gz":
        return None
    return path.stat().st_size


def find_segy_file(input_path: Path, pattern: str, recursive: bool) -> Path:
    if input_path.is_file():
        return input_path
    if not input_path.is_dir():
        raise FileNotFoundError(input_path)

    globber = input_path.rglob if recursive else input_path.glob
    candidates = sorted(p for p in globber(pattern) if p.is_file())
    if not candidates and pattern == "*.sgy*":
        patterns = ("*.sgy", "*.segy", "*.sgy.gz", "*.segy.gz", "*.SGY", "*.SEGY")
        for pat in patterns:
            candidates.extend(sorted(p for p in globber(pat) if p.is_file()))
    if not candidates:
        raise FileNotFoundError(f"No SEGY files matching {pattern!r} under {input_path}")
    return candidates[0]


def read_basic_header(path: Path) -> dict[str, object]:
    with open_binary(path) as fid:
        fid.seek(3200)
        binary_header = fid.read(400)
        if len(binary_header) < 400:
            raise ValueError("File is too small to contain a SEGY binary header.")

        sample_interval_us = int.from_bytes(binary_header[16:18], byteorder="big", signed=False)
        n_samples = int.from_bytes(binary_header[20:22], byteorder="big", signed=False)
        n_samples_extended = int.from_bytes(binary_header[84:88], byteorder="big", signed=False)
        sample_format = int.from_bytes(binary_header[24:26], byteorder="big", signed=False)
        n_extended_headers = int.from_bytes(binary_header[304:306], byteorder="big", signed=True)

    # OptaSense SEGY v1.5 sets the standard 2-byte nsamples field to zero
    # when nsamples > 32767, then stores the true value in binary bytes 85-88.
    if n_samples == 0 and n_samples_extended > 0:
        n_samples = n_samples_extended

    fs = 1.0 / (sample_interval_us * 1e-6) if sample_interval_us else 0.0
    bytes_per_sample = 4
    data_start = 3600
    if n_extended_headers > 0:
        data_start += n_extended_headers * 3200
    size = file_size_uncompressed(path)
    n_traces = None
    if size is not None and n_samples:
        n_traces = (size - data_start) // (240 + bytes_per_sample * n_samples)

    return {
        "n_samples": n_samples,
        "n_samples_standard": int.from_bytes(binary_header[20:22], byteorder="big", signed=False),
        "n_samples_extended": n_samples_extended,
        "sample_interval_us": sample_interval_us,
        "fs_hz": fs,
        "sample_format_code": sample_format,
        "n_extended_headers": n_extended_headers,
        "data_start_byte": data_start,
        "n_traces": n_traces,
    }


def header_with_dasutils(path: Path, dasutils) -> dict[str, object]:
    header = read_basic_header(path)
    if dasutils is None:
        return header

    try:
        nt, fs, start_time, end_time, n_traces = dasutils.read_PASSCAL_SEGY_headers(str(path))
    except Exception as exc:
        print(f"Warning: DASutils header read failed ({exc}); using basic SEGY header.")
        return header

    nt = int(nt)
    if nt == 0 and int(header["n_samples"]) > 0:
        print(
            "Note: DASutils PASSCAL header reports n_samples=0; "
            "using OptaSense extended sample count from binary bytes 85-88."
        )
        nt = int(header["n_samples"])
        dt = 1.0 / float(fs) if fs else 0.0
        n_traces = header.get("n_traces", n_traces)
        end_time = start_time + __import__("datetime").timedelta(seconds=nt * dt)

    header.update(
        {
            "n_samples": nt,
            "fs_hz": float(fs),
            "start_time": str(start_time),
            "end_time": str(end_time),
            "n_traces": int(n_traces),
        }
    )
    return header


def read_preview_data(
    path: Path,
    header: dict[str, object],
    dasutils,
    trace_start: int,
    trace_count: int,
) -> np.ndarray:
    n_samples = int(header["n_samples"])

    if dasutils is not None:
        try:
            return dasutils.read_PASSCAL_segy(
                str(path),
                nTraces=trace_count,
                nSample=n_samples,
                TraceOff=trace_start,
            )
        except Exception as exc:
            print(f"Warning: DASutils data read failed ({exc}); using fallback reader.")

    sample_format = int(header.get("sample_format_code", 0) or 0)
    if sample_format not in (0, 5):
        print(
            "Warning: fallback reader assumes big-endian IEEE float32. "
            f"SEGY sample format code is {sample_format}."
        )

    data = np.empty((trace_count, n_samples), dtype=np.float32)
    bytes_per_trace = 240 + n_samples * 4
    data_start = int(header.get("data_start_byte", 3600) or 3600)
    with open_binary(path) as fid:
        fid.seek(data_start + trace_start * bytes_per_trace)
        for itrace in range(trace_count):
            trace_header = fid.read(240)
            if len(trace_header) < 240:
                raise EOFError(f"Trace header ended early at trace {trace_start + itrace}.")
            buf = fid.read(n_samples * 4)
            if len(buf) < n_samples * 4:
                raise EOFError(f"Trace data ended early at trace {trace_start + itrace}.")
            data[itrace, :] = np.frombuffer(buf, dtype=">f4", count=n_samples)
    return data


def fill_missing_sample_count(
    header: dict[str, object],
    n_samples_override: int | None,
    record_seconds: float | None,
) -> None:
    n_samples = int(header.get("n_samples", 0) or 0)
    fs = float(header.get("fs_hz", 0.0) or 0.0)

    if n_samples_override is not None:
        header["n_samples"] = int(n_samples_override)
        header["n_samples_source"] = "command-line --n-samples"
    elif n_samples == 0 and record_seconds is not None and fs > 0:
        header["n_samples"] = int(round(record_seconds * fs))
        header["n_samples_source"] = f"command-line/default --record-seconds={record_seconds:g}"

    n_samples = int(header.get("n_samples", 0) or 0)
    size = file_size_uncompressed(Path(str(header["source_file"])))
    data_start = int(header.get("data_start_byte", 3600) or 3600)
    if size is not None and n_samples > 0:
        bytes_per_trace = 240 + 4 * n_samples
        n_traces = (size - data_start) // bytes_per_trace
        leftover = (size - data_start) % bytes_per_trace
        header["n_traces"] = int(n_traces)
        header["trace_size_bytes"] = int(bytes_per_trace)
        header["file_size_leftover_bytes"] = int(leftover)


def preprocess_for_plot(data: np.ndarray, median_subtract: bool, trace_normalize: bool) -> np.ndarray:
    data = np.asarray(data, dtype=np.float32)
    data = data - np.nanmean(data, axis=1, keepdims=True)
    if median_subtract:
        data = data - np.nanmedian(data, axis=0, keepdims=True)
    if trace_normalize:
        scale = np.nanpercentile(np.abs(data), 99.0, axis=1, keepdims=True)
        scale[~np.isfinite(scale) | (scale == 0)] = 1.0
        data = data / scale
    return np.nan_to_num(data)


def plot_time_distance(
    data: np.ndarray,
    fs: float,
    trace_start: int,
    title: str,
    out_png: Path,
    clip_percentile: float,
    time_max: float | None,
) -> None:
    n_traces, n_samples = data.shape
    if time_max is not None and fs > 0:
        keep = max(1, min(n_samples, int(time_max * fs)))
        data = data[:, :keep]
        n_samples = keep

    clip = np.nanpercentile(np.abs(data), clip_percentile)
    if not np.isfinite(clip) or clip == 0:
        clip = np.nanmax(np.abs(data)) or 1.0

    t1 = n_samples / fs if fs > 0 else n_samples
    extent = [0, t1, trace_start + n_traces - 1, trace_start]

    fig, ax = plt.subplots(figsize=(12, 8), constrained_layout=True)
    im = ax.imshow(
        data,
        aspect="auto",
        cmap="seismic",
        vmin=-clip,
        vmax=clip,
        extent=extent,
        interpolation="nearest",
    )
    ax.set_title(title)
    ax.set_xlabel("Time since file start (s)" if fs > 0 else "Sample")
    ax.set_ylabel("Trace index")
    fig.colorbar(im, ax=ax, label="Strain-rate amplitude (clipped)")
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_rms(data: np.ndarray, trace_start: int, title: str, out_png: Path) -> None:
    traces = trace_start + np.arange(data.shape[0])
    rms = np.sqrt(np.nanmean(data**2, axis=1))

    fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
    ax.plot(traces, rms, lw=1.2)
    ax.set_title(title)
    ax.set_xlabel("Trace index")
    ax.set_ylabel("RMS amplitude")
    ax.grid(True, alpha=0.25)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def write_header(header: dict[str, object], path: Path, source_file: Path) -> None:
    with path.open("w", encoding="utf-8") as fid:
        fid.write(f"source_file: {source_file}\n")
        for key in sorted(header):
            fid.write(f"{key}: {header[key]}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read a small preview from a SAFOD DAS SEGY file and save plots."
    )
    parser.add_argument("input", help="SEGY file or directory containing SEGY files")
    parser.add_argument("--pattern", default="*.sgy*", help="Glob pattern if input is a directory")
    parser.add_argument("--recursive", action="store_true", help="Search below input directory")
    parser.add_argument("--dasutils-path", default=None, help="Directory containing DASutils.py")
    parser.add_argument("--trace-start", type=int, default=0, help="First trace/channel to read")
    parser.add_argument("--trace-count", type=int, default=400, help="Number of traces to read")
    parser.add_argument("--n-samples", type=int, default=None, help="Override samples per trace")
    parser.add_argument(
        "--record-seconds",
        type=float,
        default=60.0,
        help="Record length used when SEGY sample-count headers are zero",
    )
    parser.add_argument("--time-max", type=float, default=10.0, help="Seconds to plot from file start")
    parser.add_argument("--clip-percentile", type=float, default=99.0, help="Symmetric plot clip percentile")
    parser.add_argument("--median-subtract", action="store_true", help="Subtract median trace at each time sample")
    parser.add_argument(
        "--trace-normalize",
        action="store_true",
        help="Scale each trace by its own 99th percentile absolute amplitude for visual QC",
    )
    parser.add_argument("--out-dir", default="segy_preview_out", help="Output directory for PNGs/header")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).expanduser()
    segy_file = find_segy_file(input_path, args.pattern, args.recursive)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dasutils = import_dasutils(args.dasutils_path)
    header = header_with_dasutils(segy_file, dasutils)
    header["source_file"] = str(segy_file)
    fill_missing_sample_count(header, args.n_samples, args.record_seconds)

    if int(header.get("n_samples", 0) or 0) <= 0:
        raise ValueError(
            "Could not determine samples per trace. Rerun with --n-samples or --record-seconds."
        )

    n_traces_total = header.get("n_traces")
    trace_count = args.trace_count
    if n_traces_total is not None:
        trace_count = min(trace_count, max(0, int(n_traces_total) - args.trace_start))
    if trace_count <= 0:
        raise ValueError("Requested trace range is empty.")

    print(f"Selected file: {segy_file}")
    print(f"Header: {header}")
    print(f"Reading traces {args.trace_start}:{args.trace_start + trace_count}")

    data = read_preview_data(segy_file, header, dasutils, args.trace_start, trace_count)
    data_plot = preprocess_for_plot(data, args.median_subtract, args.trace_normalize)

    stem = segy_file.name.replace(".gz", "").replace(".segy", "").replace(".sgy", "")
    suffix = f"tr{args.trace_start}-{args.trace_start + trace_count}"
    header_txt = out_dir / f"{stem}_{suffix}_header.txt"
    td_png = out_dir / f"{stem}_{suffix}_timedist.png"
    rms_png = out_dir / f"{stem}_{suffix}_rms.png"

    write_header(header, header_txt, segy_file)
    plot_time_distance(
        data_plot,
        float(header.get("fs_hz", 0.0) or 0.0),
        args.trace_start,
        f"{segy_file.name} traces {args.trace_start}-{args.trace_start + trace_count}",
        td_png,
        args.clip_percentile,
        args.time_max,
    )
    plot_rms(data_plot, args.trace_start, f"{segy_file.name} trace RMS", rms_png)

    print(f"Saved {header_txt}")
    print(f"Saved {td_png}")
    print(f"Saved {rms_png}")


if __name__ == "__main__":
    main()
