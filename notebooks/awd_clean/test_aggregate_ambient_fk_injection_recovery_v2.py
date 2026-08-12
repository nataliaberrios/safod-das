#!/usr/bin/env python3
"""Unit tests for ambient F-K injection-recovery aggregation helpers."""
from __future__ import annotations

import numpy as np

from aggregate_ambient_fk_injection_recovery_v2 import (
    paired_bootstrap_mean,
    paired_bootstrap_stack_score,
    physical_trajectory_values,
)


def test_physical_trajectory_sampling() -> None:
    lags = np.linspace(-0.2, 0.2, 81)
    distance = np.asarray([100.0, 200.0, 300.0])
    correlations = np.zeros((4, 3, lags.size))
    for receiver, offset in enumerate(distance):
        correlations[:, receiver, np.argmin(abs(lags - offset / 2000.0))] = receiver + 1
    values = physical_trajectory_values(correlations, lags, distance, 2000.0, 1.0)
    assert np.array_equal(values, np.tile([1.0, 2.0, 3.0], (4, 1)))


def test_paired_bootstraps_detect_known_uplift() -> None:
    rng = np.random.default_rng(7)
    baseline = rng.normal(0.0, 0.01, 100)
    injected = baseline + 0.2
    estimate, low, high = paired_bootstrap_mean(
        injected, baseline, np.random.default_rng(8), 500
    )
    assert low > 0.0 and low <= estimate <= high

    baseline_stack = rng.normal(0.0, 0.01, (100, 14))
    injected_stack = baseline_stack + 0.1
    estimate, low, high = paired_bootstrap_stack_score(
        injected_stack, baseline_stack, np.random.default_rng(9), 500
    )
    assert low > 0.0 and low <= estimate <= high


def test_stack_uplift_uses_difference_of_scores() -> None:
    """Guard against replacing score(injected)-score(baseline) by median(delta)."""
    baseline = np.tile([0.0, 100.0, 100.0], (20, 1))
    injected = np.tile([1.0, 0.0, 101.0], (20, 1))
    estimate, low, high = paired_bootstrap_stack_score(
        injected, baseline, np.random.default_rng(10), 100
    )
    expected = np.median(injected.mean(axis=0)) - np.median(baseline.mean(axis=0))
    old_wrong_value = np.median((injected - baseline).mean(axis=0))
    assert expected == -99.0 and old_wrong_value == 1.0
    assert estimate == low == high == expected


def main() -> None:
    tests = (
        test_physical_trajectory_sampling,
        test_paired_bootstraps_detect_known_uplift,
        test_stack_uplift_uses_difference_of_scores,
    )
    for test in tests:
        test()
        print(f"PASS {test.__name__}")


if __name__ == "__main__":
    main()
