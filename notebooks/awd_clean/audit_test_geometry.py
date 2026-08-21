#!/usr/bin/env python3
"""Synthetic/registration audit of safod_geometry (item D of the audit).

`tvd()` and `describe()` use np.searchsorted, which returns an INSERTION POINT,
not a nearest index. That is only safe if the channel axis is contiguous and the
requested channel is present. These checks establish whether it is, and pin the
numbers the module docstring asserts.

  D1  the channel axis is contiguous 0..N-1 with unit spacing
  D2  searchsorted returns the exact index for every integer channel
  D3  tvd(211) == 2 m and tvd(1700) == 2550 m, as documented
  D4  out-of-range and non-integer channels: what does tvd() actually return?
  D5  in_hole = md > 0.01 puts the first downhole channel at 211
  D6  near_vertical is 211-949, TVD 2-1513 m, and is CONTIGUOUS -- a disjoint
      near-vertical set would let a receiver mask silently pick up a second,
      far-away segment
  D7  return_limb_partner: TVD really matches, and it round-trips
  D8  return_limb_partner on a SURFACE channel -- outbound[j] is False for
      surface fibre too, so the branch is taken for the wrong reason

Run:  python audit_test_geometry.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import safod_geometry as geo

RESULTS = []


def record(name, verdict, detail):
    RESULTS.append((name, verdict, detail))
    print("%-6s %-46s %s" % (verdict, name, detail), flush=True)


def main():
    print("=" * 88)
    print("audit_test_geometry.py -- safod_geometry channel -> depth")
    print("=" * 88)
    d = geo.load()
    ch = d["channel"]

    # D1
    diffs = np.unique(np.diff(ch))
    ok = ch[0] == 0 and diffs.size == 1 and diffs[0] == 1
    record("D1 channel axis contiguous 0..N-1", "PASS" if ok else "FAIL",
           "n=%d, first=%g, last=%g, unique diffs=%s"
           % (ch.size, ch[0], ch[-1], diffs))

    # D2 -- searchsorted exactness for every integer channel
    idx = np.searchsorted(ch, np.arange(ch.size))
    bad = int(np.sum(idx != np.arange(ch.size)))
    record("D2 searchsorted exact for every channel",
           "PASS" if bad == 0 else "FAIL",
           "%d of %d channels map to the wrong index" % (bad, ch.size))

    # D3 -- the documented values
    t211, t1700 = float(geo.tvd(211)), float(geo.tvd(1700))
    record("D3 tvd(211) == 2 m", "PASS" if abs(t211 - 2.0) < 0.5 else "FAIL",
           "tvd(211) = %.3f m" % t211)
    record("D3b tvd(1700) == 2550 m",
           "PASS" if abs(t1700 - 2550.0) < 1.0 else "FAIL",
           "tvd(1700) = %.3f m" % t1700)

    # D4 -- out of range and non-integer
    hi = float(geo.tvd(ch.size + 500))
    lo = float(geo.tvd(-5))
    half = float(geo.tvd(400.5))
    record("D4 out-of-range channel is CLIPPED, not rejected", "WARN",
           "tvd(%d) = %.1f m (== tvd(last) = %.1f); tvd(-5) = %.1f m"
           % (ch.size + 500, hi, float(d["tvd_m"][-1]), lo))
    record("D4b non-integer channel rounds UP (insertion point)", "WARN",
           "tvd(400.5) = %.3f m, tvd(400) = %.3f, tvd(401) = %.3f"
           % (half, float(geo.tvd(400)), float(geo.tvd(401))))

    # D5 -- in_hole threshold
    first_in = int(ch[d["in_hole"]][0])
    md_prev = float(d["md_m"][first_in - 1])
    record("D5 first in-hole channel is 211",
           "PASS" if first_in == 211 else "FAIL",
           "first in_hole = %d (md %.4f m); md[%d] = %.4f m, threshold 0.01"
           % (first_in, float(d["md_m"][first_in]), first_in - 1, md_prev))

    # D6 -- near-vertical span and contiguity
    nv = geo.near_vertical_channels()
    contiguous = bool(np.all(np.diff(nv) == 1))
    tv = geo.tvd(nv)
    record("D6 near_vertical span matches the docstring",
           "PASS" if (nv.min() == 211 and abs(nv.max() - 949) <= 2) else "FAIL",
           "ch %d-%d, TVD %.0f-%.0f m (docstring: 211-949, 2-1513 m)"
           % (nv.min(), nv.max(), tv.min(), tv.max()))
    record("D6b near_vertical is CONTIGUOUS",
           "PASS" if contiguous else "FAIL",
           "%d channels, %d gaps" % (nv.size, int(np.sum(np.diff(nv) != 1))))
    # would a >5 deg channel reappear below 5 deg further down the outbound limb?
    out_ch = ch[d["outbound"]].astype(int)
    inc_out = d["inclination_deg"][d["outbound"]]
    below = out_ch[inc_out < geo.NEAR_VERTICAL_MAX_INC]
    record("D6c no second near-vertical patch deeper down",
           "PASS" if np.all(np.diff(below) == 1) else "FAIL",
           "outbound channels under %.0f deg: %d-%d, %d gaps"
           % (geo.NEAR_VERTICAL_MAX_INC, below.min(), below.max(),
              int(np.sum(np.diff(below) != 1))))

    # D7 -- return_limb_partner
    errs = []
    for c in (300, 400, 600, 800, 949):
        p = geo.return_limb_partner(c)
        errs.append((c, p, float(geo.tvd(c)), float(geo.tvd(p))))
    worst = max(abs(a - b) for _, _, a, b in errs)
    record("D7 return_limb_partner matches TVD",
           "PASS" if worst < 3.0 else "FAIL",
           "worst |TVD_out - TVD_ret| = %.2f m over %d probes" % (worst, len(errs)))
    for c, p, a, b in errs:
        print("      ch %4d (TVD %7.1f m) -> ch %4d (TVD %7.1f m)  d=%.2f m"
              % (c, a, p, b, abs(a - b)))
    back = geo.return_limb_partner(errs[1][1])
    record("D7b partner round-trips back to the outbound limb",
           "PASS" if abs(geo.tvd(back) - errs[1][2]) < 3.0 else "FAIL",
           "ch %d -> %d -> %d (TVD %.1f vs %.1f)"
           % (errs[1][0], errs[1][1], back, float(geo.tvd(back)), errs[1][2]))

    # D8 -- surface channel falls through the outbound branch
    j = int(np.searchsorted(ch, 98))
    p98 = geo.return_limb_partner(98)
    record("D8 return_limb_partner(98) on SURFACE fibre", "WARN",
           "in_hole=%s outbound=%s -> returns ch %s (TVD %.1f m); the function "
           "has no surface guard and treats lead-in as return limb"
           % (bool(d["in_hole"][j]), bool(d["outbound"][j]), p98,
              float(geo.tvd(p98)) if p98 is not None else float("nan")))

    n_fail = sum(1 for _, v, _ in RESULTS if v == "FAIL")
    print()
    print("%d checks, %d FAIL, %d WARN"
          % (len(RESULTS), n_fail, sum(1 for _, v, _ in RESULTS if v == "WARN")))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
