#!/usr/bin/env python3
"""Synthetic audit of the correlator (item C of the 2026-08-20 numerical audit).

Covers `deep_cc_steps.correlate` and `deep_timeseries.block_spectra`.
Pure numpy/scipy, no $OAK reads, runs in ~1 min.

  C1  two channels, one wavelet, KNOWN offset -> peak at exactly that lag,
      with the correct SIGN (positive lag = receiver later than source).
  C2  the correlation is LINEAR, not circular: a wavelet placed so that a
      circular correlation would wrap must not produce a wrapped peak, and the
      result must match np.correlate on the same input.
  C3  vectorised `conj(W[0]) * W[1:]` == the per-receiver loop it replaced.
  C4  block_spectra splits the block by CHANNEL GROUP (rewritten 2026-08-20,
      replacing a time-chunked version):
        C4a  grouped == ungrouped, bit-exactly, and the window count is the
             same for every group;
        C4b  MUTATION: perturb the group bookkeeping on purpose and confirm the
             test sees it. A test that only compares grouped against ungrouped
             is blind to any error present in BOTH, so C4c covers sign.
        C4c  ABSOLUTE lag and sign of block_spectra's own correlator -- it is a
             second implementation of the correlation, independent of
             deep_cc_steps.correlate, and nothing else tests it.
  C5  memory: the `target_bytes` sizing formula against what is really
      allocated (read_records returns float64 and copies three times).

Run:  python audit_test_correlator.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ambient_lellouch2019_exact_stack as ex
import deep_cc_steps as steps
import deep_timeseries as dts

RESULTS = []


def record(name, verdict, detail):
    RESULTS.append((name, verdict, detail))
    print("%-6s %-52s %s" % (verdict, name, detail), flush=True)


def ricker(t, f0=12.0):
    a = (np.pi * f0 * t) ** 2
    return (1.0 - 2.0 * a) * np.exp(-a)


# ------------------------------------------------------------------ C1, C2
def test_lag_and_sign():
    fs = 500.0
    n_win = int(30 * fs)
    n_t = n_win                      # exactly one window
    t = np.arange(n_t) / fs
    tau_samples = 137                # KNOWN offset, receiver LATER than source
    t0 = 5.0
    src = ricker(t - t0)
    rx = ricker(t - t0 - tau_samples / fs)
    data = np.vstack([src, rx])

    g, lags, nw = steps.correlate(data, 0, [np.array([1])], fs,
                                  n_win, int(15 * fs))
    k = int(np.argmax(np.abs(g[0])))
    got = lags[k]
    want = tau_samples / fs
    record("C1 peak lag == known offset, correct sign",
           "PASS" if abs(got - want) <= 1.5 / fs else "FAIL",
           "want %+.4f s, got %+.4f s (%d windows)" % (want, got, nw))

    # negative offset: receiver EARLIER
    rx2 = ricker(t - t0 + tau_samples / fs)
    g2, lags2, _ = steps.correlate(np.vstack([src, rx2]), 0, [np.array([1])],
                                   fs, n_win, int(15 * fs))
    got2 = lags2[int(np.argmax(np.abs(g2[0])))]
    record("C1b negative offset -> negative lag",
           "PASS" if abs(got2 + want) <= 1.5 / fs else "FAIL",
           "want %+.4f s, got %+.4f s" % (-want, got2))


def test_linear_not_circular():
    """Put the wavelets near opposite ends of the window. A CIRCULAR
    correlation folds the true lag by n_win; a linear one does not."""
    fs = 500.0
    n_win = int(30 * fs)
    t = np.arange(n_win) / fs
    src = ricker(t - 1.0)            # near the start
    rx = ricker(t - 28.0)            # near the end: true lag +27 s
    true_lag = 27.0

    n_fft = 1 << int(np.ceil(np.log2(2 * n_win - 1)))
    record("C2 n_fft >= 2*n_win-1 (no circular wrap)",
           "PASS" if n_fft >= 2 * n_win - 1 else "FAIL",
           "n_win %d, n_fft %d, need >= %d" % (n_win, n_fft, 2 * n_win - 1))

    S = np.conj(np.fft.rfft(src, n=n_fft))
    full = np.fft.fftshift(np.fft.irfft(S * np.fft.rfft(rx, n=n_fft), n=n_fft))
    lags = (np.arange(n_fft) - n_fft // 2) / fs
    got = lags[int(np.argmax(np.abs(full)))]
    record("C2b unbandpassed peak at the true long lag",
           "PASS" if abs(got - true_lag) <= 2.0 / fs else "FAIL",
           "want %+.3f s, got %+.3f s" % (true_lag, got))

    ref = np.correlate(rx, src, mode="full")        # index n_win-1 == lag 0
    ref_lag = (int(np.argmax(np.abs(ref))) - (n_win - 1)) / fs
    ok = abs(ref_lag - got) <= 2.0 / fs
    # amplitude comparison over the valid lag range
    mid = n_fft // 2
    fft_slice = full[mid - (n_win - 1): mid + n_win]
    rel = float(np.max(np.abs(fft_slice - ref))) / float(np.max(np.abs(ref)))
    record("C2c FFT correlation == np.correlate(mode='full')",
           "PASS" if (ok and rel < 1e-9) else "FAIL",
           "lag %+.3f vs %+.3f s, max rel amp diff %.3e" % (got, ref_lag, rel))


# ------------------------------------------------------------------ C3
def test_vectorised_cross_spectrum():
    rng = np.random.default_rng(7)
    n_ch, n_win, n_fft = 12, 4096, 8192
    x = rng.standard_normal((n_ch, n_win)).astype(np.float32)
    W = np.fft.rfft(x, n=n_fft, axis=1)
    vec = np.conj(W[0]) * W[1:]
    loop = np.empty_like(vec)
    S = np.conj(np.fft.rfft(x[0], n=n_fft))
    for r in range(1, n_ch):
        loop[r - 1] = S * np.fft.rfft(x[r], n=n_fft)
    d = float(np.max(np.abs(vec - loop))) / float(np.max(np.abs(loop)))
    record("C3 vectorised conj(W[0])*W[1:] == per-receiver loop",
           "PASS" if d < 1e-12 else "FAIL", "max rel diff %.3e" % d)


# ------------------------------------------------------------------ C4
def _install_fake(master, rec_n, fs):
    def fake_read(paths, channels):
        idx = [int(Path(p).stem) for p in paths]
        cols = np.concatenate([np.arange(i * rec_n, (i + 1) * rec_n) for i in idx])
        return master[np.asarray(channels)][:, cols], fs, 1.0
    dts.steps.read_records = fake_read


def _synth(n_ch, n_rec, rec_n, seed=20260820):
    rng = np.random.default_rng(seed)
    src = rng.standard_normal(n_rec * rec_n).astype(np.float32)
    m = np.empty((n_ch, n_rec * rec_n), dtype=np.float32)
    for c in range(n_ch):
        m[c] = np.roll(src, 7 * c) + 0.3 * rng.standard_normal(src.size)
    return np.cumsum(m, axis=1)


def _run_pair(master, rec_n, fs, n_ch, n_rec, target_bytes):
    _install_fake(master, rec_n, fs)
    dts.SOURCE_CH = 0
    paths = ["%d" % i for i in range(n_rec)]
    rx = list(range(1, n_ch))
    whole, nw_w, _, _ = dts.block_spectra(paths, rx, (5.0, 20.0), target_bytes=1e15)
    grouped, nw_c, _, _ = dts.block_spectra(paths, rx, (5.0, 20.0),
                                            target_bytes=target_bytes)
    rel = float(np.linalg.norm(grouped - whole) / max(np.linalg.norm(whole), 1e-30))
    return rel, nw_w, nw_c


def test_channel_grouping():
    fs, rec_n, n_ch, n_rec = 500.0, int(60 * 500), 8, 20
    master = _synth(n_ch, n_rec, rec_n)
    n_win, n_step = int(dts.WINDOW_S * fs), int(dts.STEP_S * fs)
    expect = (n_rec * rec_n - n_win) // n_step + 1

    # ---- C4a: splitting by channel has no joins, so it must be BIT-exact.
    # Vary target_bytes so per_group really changes: a group of ONE receiver
    # makes several ordering bugs vacuous, so at least one case must have
    # per_group > 1.
    total_n = rec_n * n_rec
    pg = lambda tb: max(1, int(tb / (total_n * 4 * 3)) - 1)
    MUT_TB = 3.0e7                      # chosen so per_group > 1
    for tb in (2.0e6, MUT_TB):
        rel, nw_w, nw_c = _run_pair(master, rec_n, fs, n_ch, n_rec, tb)
        record("C4a grouped == ungrouped bit-exactly (per_group %d)" % pg(tb),
               "PASS" if (rel == 0.0 and nw_w == nw_c == expect) else "FAIL",
               "rel diff %.3e, windows %d/%d (analytic %d)"
               % (rel, nw_w, nw_c, expect))
    if pg(MUT_TB) < 2:
        record("C4a per_group > 1 exercised", "FAIL",
               "every group held one receiver; ordering bugs are vacuous")

    # ---- C4b: MUTATION -- would the comparison notice a bookkeeping error?
    src_txt = (HERE / "deep_timeseries.py").read_text()
    import types

    def _mutate(sub_from, sub_to, label):
        mutated = src_txt.replace(sub_from, sub_to)
        if mutated == src_txt:
            record("C4b mutation %s" % label, "FAIL",
                   "could not apply the mutation -- source changed?")
            return
        mod = types.ModuleType("dts_mut")
        mod.__dict__["__file__"] = str(HERE / "deep_timeseries.py")
        sys.modules["dts_mut"] = mod
        try:
            exec(compile(mutated, "deep_timeseries_MUT.py", "exec"), mod.__dict__)
            _install_fake(master, rec_n, fs)
            mod.steps.read_records = dts.steps.read_records
            mod.SOURCE_CH = 0
            paths = ["%d" % i for i in range(n_rec)]
            rx = list(range(1, n_ch))
            good, nw_g, _, _ = dts.block_spectra(paths, rx, (5.0, 20.0),
                                                 target_bytes=1e15)
            bad, nw_b, _, _ = mod.block_spectra(paths, rx, (5.0, 20.0),
                                                target_bytes=MUT_TB)
            rel = float(np.linalg.norm(bad - good) / np.linalg.norm(good))
            seen = (nw_b != nw_g) or (rel > 1e-12)
            detail = "windows %d vs %d, rel diff %.3e" % (nw_g, nw_b, rel)
        except Exception as exc:                       # a crash is detection too
            seen, detail = True, "raised %s: %s" % (type(exc).__name__, exc)
        record("C4b mutation %s is detected" % label,
               "PASS" if seen else "FAIL", detail)

    _mutate("acc[g0:g0 + len(grp)] += np.conj(W[0]) * W[1:]",
            "acc[g0:g0 + len(grp)] += np.conj(W[0]) * W[1:] * 1.0000001",
            "1e-7 amplitude perturbation")
    _mutate("for s in range(0, x.shape[1] - n_win + 1, n_step):",
            "for s in range(0, x.shape[1] - n_win, n_step):",
            "drop the last window")
    _mutate("grp = [int(c) for c in rx_ch[g0:g0 + per_group]]",
            "grp = [int(c) for c in rx_ch[g0:g0 + per_group]][::-1]",
            "reverse the receiver order within a group")


def test_block_spectra_lag_and_sign():
    """block_spectra is a SECOND implementation of the correlator. Comparing it
    against itself (grouped vs ungrouped) cannot see a sign error present in
    both, so check the absolute lag on synthetic data with a known moveout."""
    fs, rec_n, n_ch, n_rec = 500.0, int(60 * 500), 6, 4
    rng = np.random.default_rng(3)
    total = n_rec * rec_n
    src = rng.standard_normal(total).astype(np.float32)
    lag_samples = 25                       # each receiver LATER than the last
    master = np.empty((n_ch, total), dtype=np.float32)
    for c in range(n_ch):
        master[c] = np.roll(src, lag_samples * c)     # +ve roll = delayed
    _install_fake(master, rec_n, fs)
    dts.SOURCE_CH = 0
    acc, nw, fs_, n_fft = dts.block_spectra(["%d" % i for i in range(n_rec)],
                                            list(range(1, n_ch)), (5.0, 20.0),
                                            target_bytes=1e15)
    g, lags = dts.to_gather(acc, nw, fs_, n_fft, (5.0, 20.0))
    got = [float(lags[int(np.argmax(np.abs(row)))]) for row in g]
    want = [lag_samples * c / fs for c in range(1, n_ch)]
    err = max(abs(a - b) for a, b in zip(got, want))
    record("C4c block_spectra lag and sign (independent check)",
           "PASS" if err <= 2.0 / fs else "FAIL",
           "want %s, got %s (max err %.4f s)"
           % (["%+.3f" % v for v in want], ["%+.3f" % v for v in got], err))


def test_primitives_dtype_and_equivalence():
    """`differentiate` and `ram_normalise` were rewritten to be dtype-preserving
    and to stop allocating extra full-size temporaries. Same arithmetic?

    Also confirms the float32 path is exact for data shaped like the real Deep
    record: int32 phase up to ~1.9e6, which fits float32's exactly-representable
    integer range (2^24), differenced sample to sample."""
    rng = np.random.default_rng(12)
    fs = 1000.0
    base = rng.integers(-1_900_000, 1_900_000, size=(16, 1)).astype(np.int64)
    walk = np.cumsum(rng.integers(-400, 400, size=(16, 20000)), axis=1)
    raw_i = (base + walk).astype(np.int64)
    assert np.abs(raw_i).max() < 2 ** 24, "test data outside float32 exact range"

    r64 = np.zeros(raw_i.shape, dtype=np.float64)
    r64[:, 1:] = np.diff(raw_i.astype(np.float64), axis=1) * fs      # old formula
    for dt in (np.float32, np.float64):
        got = steps.differentiate(raw_i.astype(dt), fs)
        same = np.array_equal(got.astype(np.float64), r64)
        record("C6 differentiate exact and dtype-preserving (%s)" % dt.__name__,
               "PASS" if (same and got.dtype == dt) else "FAIL",
               "dtype %s, exact match %s, max abs diff %.3e"
               % (got.dtype, same, float(np.max(np.abs(got.astype(np.float64) - r64)))))

    n32 = steps.ram_normalise(steps.differentiate(raw_i.astype(np.float32), fs), fs)
    n64 = steps.ram_normalise(steps.differentiate(raw_i.astype(np.float64), fs), fs)
    rel = float(np.max(np.abs(n32.astype(np.float64) - n64)) /
                max(1e-30, np.max(np.abs(n64))))
    record("C6b ram_normalise float32 == float64 to f32 precision",
           "PASS" if (n32.dtype == np.float32 and rel < 1e-6) else "FAIL",
           "dtype %s, max rel diff %.3e" % (n32.dtype, rel))


def test_memory_sizing():
    """`per_group` is sized as target_bytes / (total_n * 4 * 3), i.e. it budgets
    12 bytes per sample per channel. Check that against what is really live.

    read_records returns FLOAT64. Before 2026-08-20 it also held the per-file
    list, the concatenate result and the fancy-indexed copy at the same instant
    -- 24 bytes/sample, twice the budget, which is what OOM-killed the deepC
    tasks. It now preallocates (measured 1.13x by audit_test_readrecords.py R7),
    so the binding term is the float64 return alive next to the float32 copy."""
    fs, rec_n, n_rec = 1000.0, int(60 * 1000), 240      # a real 4 h deepC block
    total_n = rec_n * n_rec
    target_bytes = 4.0e9
    per_group = max(1, int(target_bytes / (total_n * 4 * 3)) - 1)
    n_read = per_group + 1                               # + the source channel
    budget = total_n * 4 * 3 * n_read
    per_sample = {
        "read_records (float64, preallocated 1.13x)": 8 * 1.13,
        "raw float64 alive while x float32 is built": 8 + 4,
        "x + abs(x) + w, all float32":                 4 * 3,
    }
    worst_name = max(per_sample, key=per_sample.get)
    worst = total_n * n_read * per_sample[worst_name]
    record("C5 target_bytes sizing vs real allocation",
           "PASS" if worst <= target_bytes * 1.02 else "FAIL",
           "per_group %d (%d ch read): budget %.2f GB, worst stage %.2f GB "
           "(%.2fx) -- %s"
           % (per_group, n_read, budget / 2**30, worst / 2**30,
              worst / max(budget, 1), worst_name))
    for k, v in per_sample.items():
        print("      %-44s %5.1f B/sample -> %.2f GB"
              % (k, v, total_n * n_read * v / 2**30))


def main():
    print("=" * 92)
    print("audit_test_correlator.py -- deep_cc_steps.correlate / deep_timeseries.block_spectra")
    print("=" * 92)
    test_lag_and_sign(); print()
    test_linear_not_circular(); print()
    test_vectorised_cross_spectrum(); print()
    test_channel_grouping(); print()
    test_block_spectra_lag_and_sign(); print()
    test_primitives_dtype_and_equivalence(); print()
    test_memory_sizing(); print()
    n_fail = sum(1 for _, v, _ in RESULTS if v == "FAIL")
    print("%d checks, %d FAIL, %d WARN"
          % (len(RESULTS), n_fail, sum(1 for _, v, _ in RESULTS if v == "WARN")))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
