#!/usr/bin/env python3
"""Synthetic audit of deep_cc_steps.read_records -- HDF5 fancy indexing + memory.

Writes small HDF5 files into a temporary directory (no $OAK), so it runs in
seconds anywhere.

  R1  ascending channels: rows come back in the requested order
  R2  DESCENDING channels (the Deep return limb): rows come back in the
      requested order, not the sorted order. This is the bug class that
      silently broke every return-limb read before the argsort was added.
  R3  arbitrary/shuffled order
  R4  duplicate channels are rejected rather than silently deduplicated
  R5  records of UNEQUAL length concatenate correctly (the preallocation reads
      each dataset's shape rather than assuming the first record's length)
  R6  values match a direct h5py read, exactly
  R7  peak allocation: the returned array is the only full-size copy

Run:  python audit_test_readrecords.py
"""
from __future__ import annotations

import sys
import tempfile
import tracemalloc
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import deep_cc_steps as steps

RESULTS = []


def record(name, verdict, detail):
    RESULTS.append((name, verdict, detail))
    print("%-6s %-48s %s" % (verdict, name, detail), flush=True)


def write_file(path, data, fs=500.0, dx=2.0419):
    """data is (n_time, n_channel), matching the OptaSense layout."""
    with h5py.File(path, "w") as h:
        acq = h.create_group("Acquisition")
        acq.attrs["SpatialSamplingInterval"] = dx
        raw = acq.create_group("Raw[0]")
        raw.attrs["OutputDataRate"] = fs
        raw.create_dataset("RawData", data=data)


def main():
    print("=" * 84)
    print("audit_test_readrecords.py -- deep_cc_steps.read_records")
    print("=" * 84)
    rng = np.random.default_rng(20260820)
    n_ch = 40
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        lens = [500, 500, 300]                       # deliberately unequal
        truth = []
        paths = []
        for i, n in enumerate(lens):
            d = rng.integers(-30000, 30000, size=(n, n_ch)).astype(np.int32)
            truth.append(d)
            p = td / ("rec%d.h5" % i)
            write_file(p, d)
            paths.append(str(p))
        full = np.concatenate(truth, axis=0).T.astype(np.float64)   # (ch, time)

        for label, chans in (("R1 ascending", [3, 7, 11, 19]),
                             ("R2 DESCENDING (return limb)", [19, 11, 7, 3]),
                             ("R3 shuffled", [11, 3, 19, 7])):
            got, fs, dx = steps.read_records(paths, chans)
            want = full[np.asarray(chans)]
            ok = got.shape == want.shape and np.array_equal(got, want)
            record(label, "PASS" if ok else "FAIL",
                   "shape %s, exact match %s, fs %.0f, dx %.4f"
                   % (got.shape, ok, fs, dx))

        # R4 duplicates
        try:
            steps.read_records(paths, [5, 9, 5])
            record("R4 duplicate channels rejected", "FAIL", "no exception raised")
        except ValueError as e:
            record("R4 duplicate channels rejected", "PASS", "ValueError: %s" % e)

        # R5 unequal record lengths
        got, _, _ = steps.read_records(paths, [2, 1, 0])
        record("R5 unequal record lengths concatenate",
               "PASS" if got.shape[1] == sum(lens) else "FAIL",
               "got %d samples, expected %d %s" % (got.shape[1], sum(lens), lens))

        # R6 exact against a direct h5py read
        direct = np.concatenate(
            [h5py.File(p, "r")["Acquisition/Raw[0]/RawData"][:, [8]].T for p in paths],
            axis=1).astype(np.float64)
        got, _, _ = steps.read_records(paths, [8])
        record("R6 values match a direct h5py read",
               "PASS" if np.array_equal(got, direct) else "FAIL",
               "max abs diff %.3e" % float(np.max(np.abs(got - direct))))

        # R7 peak allocation
        big_n, big_ch, n_files = 40000, 24, 4
        bp = []
        for i in range(n_files):
            d = rng.integers(-30000, 30000, size=(big_n, big_ch)).astype(np.int32)
            p = td / ("big%d.h5" % i)
            write_file(p, d)
            bp.append(str(p))
        want_bytes = big_ch * big_n * n_files * 8
        tracemalloc.start()
        out, _, _ = steps.read_records(bp, list(range(big_ch)))
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        ratio = peak / want_bytes
        del out
        record("R7 peak allocation ~ 1x the returned array",
               "PASS" if ratio < 1.6 else "FAIL",
               "returned %.1f MB, traced peak %.1f MB (%.2fx)"
               % (want_bytes / 2**20, peak / 2**20, ratio))

    n_fail = sum(1 for _, v, _ in RESULTS if v == "FAIL")
    print()
    print("%d checks, %d FAIL" % (len(RESULTS), n_fail))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
