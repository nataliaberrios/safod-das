#!/usr/bin/env python3
"""Deterministic unit tests for ambient F-K injection recovery."""
from __future__ import annotations

import numpy as np

from ambient_fk_injection_recovery_v1 import (
    broadband_plane_wave,
    parse_amplitudes,
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


def main() -> None:
    tests = (
        test_amplitude_parser_requires_zero,
        test_broadband_wave_unit_rms_and_reproducibility,
        test_prefilter_recovers_both_directions_and_rejects_off_target,
    )
    for test in tests:
        test()
        print(f"PASS {test.__name__}")


if __name__ == "__main__":
    main()
