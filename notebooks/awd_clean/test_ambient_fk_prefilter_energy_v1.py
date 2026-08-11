#!/usr/bin/env python3
"""Deterministic tests for the independent pre-filter F-K energy statistic."""
from __future__ import annotations

import numpy as np

from ambient_fk_prefilter_energy_v1 import (
    enrichment,
    segment_statistics,
    spatial_power,
    supports,
    time_spectrum,
)


def plane_wave(
    velocity_m_s: float,
    *,
    frequency_hz: float = 10.0,
    nchannels: int = 256,
    nsamples: int = 4096,
    dx_m: float = 2.0,
    fs_hz: float = 250.0,
    direction: float = 1.0,
) -> tuple[np.ndarray, float, float]:
    """Return a unit plane wave using the production coordinate convention."""
    distance = np.arange(nchannels, dtype=float)[:, None] * dx_m
    time = np.arange(nsamples, dtype=float)[None, :] / fs_hz
    phase = 2.0 * np.pi * frequency_hz * (
        time - direction * distance / velocity_m_s
    )
    return np.cos(phase), fs_hz, dx_m


def direct_enrichment(x: np.ndarray, fs: float, dx: float, branch: str) -> float:
    """Evaluate the pre-filter statistic without any surrogate construction."""
    frequency, channel_spectrum = time_spectrum(x, fs)
    wavenumber, power = spatial_power(channel_spectrum, dx)
    return enrichment(power, *supports(frequency, wavenumber, branch))


def test_target_direction_and_off_target_rejection() -> None:
    target, fs, dx = plane_wave(3200.0)
    off_target, _, _ = plane_wave(1800.0)
    target_negative = direct_enrichment(target, fs, dx, "negative")
    target_positive = direct_enrichment(target, fs, dx, "positive")
    off_target_negative = direct_enrichment(off_target, fs, dx, "negative")
    assert target_negative > target_positive + 0.5
    assert target_negative > off_target_negative + 2.0


def test_target_survives_prefilter_null_test() -> None:
    target, fs, dx = plane_wave(3200.0)
    observed_values = []
    null_values = {
        "channel_permutation": [],
        "circular_time_shift": [],
    }
    for file_index in range(8):
        observed, null, _ = segment_statistics(
            target,
            fs,
            dx,
            f"synthetic-target-{file_index}",
            seed=20260810,
            null_start=0,
            nulls=20,
        )
        observed_values.append(observed["negative"])
        for method in null_values:
            null_values[method].append(null[method]["negative"])
    observed_mean = float(np.mean(observed_values))
    for method, values in null_values.items():
        null_mean = np.mean(values, axis=0)
        assert observed_mean > np.max(null_mean), method


def test_noise_has_no_forced_target_enrichment() -> None:
    rng = np.random.default_rng(23)
    noise = rng.standard_normal((256, 4096))
    value = direct_enrichment(noise, 250.0, 2.0, "negative")
    assert abs(value) < 0.5


def main() -> None:
    tests = (
        test_target_direction_and_off_target_rejection,
        test_target_survives_prefilter_null_test,
        test_noise_has_no_forced_target_enrichment,
    )
    for test in tests:
        test()
        print(f"PASS {test.__name__}")


if __name__ == "__main__":
    main()
