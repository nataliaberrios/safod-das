#!/usr/bin/env python3
"""Deterministic regression tests for ambient correlation and signed F-K lags."""
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from ambient_fk_transfer_test import fk_filter
from ambient_transfer_test import normalized_corr_pairs


def test_impulse_lag_reversal():
    fs = 100.0
    data = np.zeros((2, 400), dtype=float)
    data[0, 100] = 1.0
    data[1, 120] = 1.0
    lags, forward = normalized_corr_pairs(data, [(0, 1)], fs, max_lag=0.3)
    _, reverse = normalized_corr_pairs(data, [(1, 0)], fs, max_lag=0.3)
    forward_peak = float(lags[np.argmax(forward[0])])
    reverse_peak = float(lags[np.argmax(reverse[0])])
    assert np.isclose(forward_peak, 0.2), forward_peak
    assert np.isclose(reverse_peak, -0.2), reverse_peak


def broadband_plane_wave(direction, fs, dx, nz, nt, velocity):
    rng = np.random.default_rng(20260805)
    frequency = np.fft.rfftfreq(nt, 1.0 / fs)
    spectrum = np.fft.rfft(rng.standard_normal(nt))
    spectrum[(frequency < 7.0) | (frequency > 17.0)] = 0.0
    coordinate = np.arange(nz) * dx
    phase = np.exp(
        -1j * 2.0 * np.pi * frequency[None, :] * direction
        * coordinate[:, None] / velocity
    )
    return np.fft.irfft(spectrum[None, :] * phase, n=nt, axis=1)


def branch_check(direction, retained_mode, rejected_mode, expected_sign):
    fs, dx, velocity = 200.0, 1.0, 3200.0
    wave = broadband_plane_wave(direction, fs, dx, nz=401, nt=4096, velocity=velocity)
    retained, filtered_fs, filtered_dx = fk_filter(wave, fs, dx, retained_mode)
    rejected, _, _ = fk_filter(wave, fs, dx, rejected_mode)
    retained_rms = float(np.sqrt(np.mean(retained ** 2)))
    rejected_rms = float(np.sqrt(np.mean(rejected ** 2)))
    assert retained_rms > 8.0 * rejected_rms, (retained_rms, rejected_rms)

    targets = [20, 40, 80]
    lags, correlations = normalized_corr_pairs(
        retained, [(0, target) for target in targets], filtered_fs, max_lag=0.15
    )
    observed = lags[np.argmax(correlations, axis=1)]
    expected = expected_sign * np.asarray(targets) * filtered_dx / velocity
    assert np.all(np.abs(observed - expected) <= 1.0 / filtered_fs + 1e-12), (
        observed,
        expected,
    )


def test_increasing_coordinate_negative_mask_positive_lag():
    branch_check(direction=1.0, retained_mode="negative", rejected_mode="positive", expected_sign=1.0)


def test_decreasing_coordinate_positive_mask_negative_lag():
    branch_check(direction=-1.0, retained_mode="positive", rejected_mode="negative", expected_sign=-1.0)


def main():
    tests = [
        test_impulse_lag_reversal,
        test_increasing_coordinate_negative_mask_positive_lag,
        test_decreasing_coordinate_positive_mask_negative_lag,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")


if __name__ == "__main__":
    main()
