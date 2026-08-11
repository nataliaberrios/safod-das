#!/usr/bin/env python3
"""Small deterministic tests for the pre-F-K surrogate null machinery."""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from ambient_fk_full_pipeline_null_v2 import (
    channel_permutation,
    circular_time_shift,
    stable_rng,
)


def test_stable_rng() -> None:
    first = stable_rng(7, "channel_permutation", 3, "file-a").integers(0, 10_000, 20)
    second = stable_rng(7, "channel_permutation", 3, "file-a").integers(0, 10_000, 20)
    different = stable_rng(7, "channel_permutation", 4, "file-a").integers(0, 10_000, 20)
    assert np.array_equal(first, second)
    assert not np.array_equal(first, different)


def test_channel_permutation_preserves_source_and_trace_inventory() -> None:
    x = np.arange(8 * 12, dtype=float).reshape(8, 12)
    surrogate = channel_permutation(x, np.random.default_rng(11))
    assert np.array_equal(surrogate[0], x[0])
    assert sorted(map(tuple, surrogate[1:])) == sorted(map(tuple, x[1:]))
    assert not np.array_equal(surrogate[1:], x[1:])


def test_circular_shift_preserves_source_and_channel_spectra() -> None:
    rng = np.random.default_rng(12)
    x = rng.standard_normal((10, 128))
    surrogate = circular_time_shift(x, np.random.default_rng(13))
    assert np.array_equal(surrogate[0], x[0])
    assert np.allclose(np.abs(np.fft.rfft(surrogate, axis=1)), np.abs(np.fft.rfft(x, axis=1)))
    assert not np.array_equal(surrogate[1:], x[1:])


def main() -> None:
    tests = [
        test_stable_rng,
        test_channel_permutation_preserves_source_and_trace_inventory,
        test_circular_shift_preserves_source_and_channel_spectra,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")


if __name__ == "__main__":
    main()
