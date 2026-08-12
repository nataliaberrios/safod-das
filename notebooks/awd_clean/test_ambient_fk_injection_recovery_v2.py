#!/usr/bin/env python3
"""Deterministic unit tests for ambient F-K injection recovery."""
from __future__ import annotations

import numpy as np

from ambient_fk_injection_recovery_v2 import (
    broadband_plane_wave,
    parse_amplitudes,
    preprocess_injected,
    prefilter_enrichment,
)
from ambient_transfer_test import preprocess


def test_amplitude_parser_requires_zero() -> None:
    assert parse_amplitudes("0.1,0,0.01") == (0.0, 0.01, 0.1)
    try:
        parse_amplitudes("0.01,0.1")
    except ValueError:
        pass
    else:
        raise AssertionError("missing zero baseline was accepted")


def test_broadband_wave_unit_rms_and_reproducibility() -> None:
    args = (128, 4096, 250.0, 2.0, 3200.0, 1, 42)
    first = broadband_plane_wave(*args)
    second = broadband_plane_wave(*args)
    assert np.array_equal(first, second)
    channel_rms = np.sqrt(np.mean(first ** 2, axis=1))
    assert np.allclose(channel_rms, 1.0, rtol=2e-6, atol=2e-6)


def test_prefilter_recovers_both_directions_and_rejects_off_target() -> None:
    fs, dx, nz, nt = 250.0, 2.0, 256, 4096
    increasing = broadband_plane_wave(nz, nt, fs, dx, 3200.0, 1, 1)
    decreasing = broadband_plane_wave(nz, nt, fs, dx, 3200.0, -1, 2)
    off_target = broadband_plane_wave(nz, nt, fs, dx, 1800.0, 1, 3)
    inc = prefilter_enrichment(preprocess(increasing, fs), fs, dx)
    dec = prefilter_enrichment(preprocess(decreasing, fs), fs, dx)
    off = prefilter_enrichment(preprocess(off_target, fs), fs, dx)
    assert inc["negative"] > inc["positive"] + 0.5
    assert dec["positive"] > dec["negative"] + 0.5
    assert inc["negative"] > off["negative"] + 1.0


def test_subcount_injection_survives_large_float32_offset() -> None:
    """A small wave must not disappear against a large raw integer offset."""
    fs, dx, nz, nt = 100.0, 2.0, 16, 2048
    coordinate = np.arange(nz, dtype=np.float32)[:, None]
    time = np.arange(nt, dtype=np.float32)[None, :]
    raw = np.asarray(2.0e8 + 32.0 * coordinate + 64.0 * time / nt, dtype=np.float32)
    wave = broadband_plane_wave(nz, nt, fs, dx, 3200.0, 1, 44)
    baseline = preprocess_injected(raw, wave, 0.0, 50.0, fs)
    injected = preprocess_injected(raw, wave, 0.003, 50.0, fs)
    assert baseline.dtype == np.float64 and injected.dtype == np.float64
    assert np.sqrt(np.mean((injected - baseline) ** 2)) > 1e-5


def test_zero_injection_matches_production_preprocess() -> None:
    rng = np.random.default_rng(45)
    raw = rng.normal(size=(32, 4096)).astype(np.float32)
    wave = broadband_plane_wave(32, 4096, 250.0, 2.0, 3200.0, 1, 46)
    expected = preprocess(raw, 250.0)
    actual = preprocess_injected(raw, wave, 0.0, 1.0, 250.0)
    assert np.allclose(actual, expected, rtol=1e-6, atol=1e-7)


def main() -> None:
    tests = (
        test_amplitude_parser_requires_zero,
        test_broadband_wave_unit_rms_and_reproducibility,
        test_prefilter_recovers_both_directions_and_rejects_off_target,
        test_subcount_injection_survives_large_float32_offset,
        test_zero_injection_matches_production_preprocess,
    )
    for test in tests:
        test()
        print(f"PASS {test.__name__}")


if __name__ == "__main__":
    main()
