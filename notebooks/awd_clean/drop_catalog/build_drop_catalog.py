"""Build the clean AWD weight-drop catalogue: node timing -> DAS detection.

What this assembles is the *DAS-side* work. The drop times themselves were
delivered, not derived here (see README.md "Provenance"). This joins them to
what the DAS actually recorded.

Inputs (all read-only):
  ../Check shots/p26.cc9.txt   node 453009664 CC picks, 989 drops   [DELIVERED]
  ../Check shots/p26.cc4.txt   node 453001432 CC picks, 989 drops   [DELIVERED]
  ../awd_manifest.csv          drop times x DAS file coverage       (build_manifest.py)
  ../nano_drop_repeatability.csv  per-drop Nano detection metrics   (nano_drop_repeatability.py)

Outputs (this directory):
  awd_drop_catalog.csv         one row per drop: timing, coverage, Nano detection
  awd_burst_summary.csv        one row per burst
  timing_uncertainty.txt       two-node agreement = empirical timing error

Times are UTC; local is PDT (UTC-7) for the June 2026 survey. Survey extent is
reported as duration, not file counts.
"""
import csv
import datetime as dt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
AWD = HERE.parent
CC9 = AWD / "Check shots" / "p26.cc9.txt"
CC4 = AWD / "Check shots" / "p26.cc4.txt"
MANIFEST = AWD / "awd_manifest.csv"
NANO = AWD / "nano_drop_repeatability.csv"
LOCAL_OFFSET_H = -7.0  # PDT during the 2026-06-16/17 survey

# Detection thresholds, fixed here so the flags are reproducible.
NCC_DETECT = 0.90       # leave-one-out signal NCC against the burst stack
SNR_DETECT_DB = 10.0    # beam SNR


def load_picks(path):
    times, ccs, srcs = [], [], []
    with path.open() as stream:
        next(stream)
        for line in stream:
            f = line.split()
            if len(f) >= 3:
                srcs.append(f[0])
                times.append(dt.datetime.fromisoformat(f[1]))
                ccs.append(float(f[2]))
    return np.array(times), np.array(ccs), srcs


PICK_QUANTUM_MS = 0.5  # node picks are quantised to 0.5 ms


def match_nearest(t_a, t_b, tol_s=1.0):
    """Signed offset t_a - nearest(t_b), in ms, snapped to the pick quantum.

    The picks are quantised to PICK_QUANTUM_MS, so raw float offsets carry
    ~1e-4 ms of representation noise. Left unsnapped, a threshold that lands on
    a grid value (2, 10, 50, 100 ms all do) splits the drops sitting exactly on
    it arbitrarily -- that is the whole difference between the 68 and 162 counts
    two implementations of this comparison have produced. Snap first, then use
    >= so the count is well defined.
    """
    sec_b = np.array([x.timestamp() for x in t_b])
    out = []
    for a in t_a:
        d = sec_b - a.timestamp()
        j = int(np.argmin(np.abs(d)))
        if abs(d[j]) > tol_s:
            out.append(None)
            continue
        ms = -d[j] * 1e3
        out.append(round(ms / PICK_QUANTUM_MS) * PICK_QUANTUM_MS)
    return out


def main():
    t9, c9, src9 = load_picks(CC9)
    t4, c4, _ = load_picks(CC4)
    offsets_ms = match_nearest(t9, t4)

    manifest = {}
    with MANIFEST.open() as s:
        for row in csv.DictReader(s):
            manifest[dt.datetime.fromisoformat(row["utc_time"]).replace(tzinfo=None)] = row

    nano = {}
    with NANO.open() as s:
        for row in csv.DictReader(s):
            k = dt.datetime.fromisoformat(row["utc_date_from_p26"]).replace(tzinfo=None)
            nano[k] = row

    catalog = []
    for i, (t, cc, src) in enumerate(zip(t9, c9, src9)):
        m = manifest.get(t)
        n = nano.get(t)
        off = offsets_ms[i]
        ncc = float(n["loo_signal_ncc"]) if n else None
        snr = float(n["beam_snr_db"]) if n else None
        catalog.append({
            "drop_index": i,
            "utc_time": t.isoformat(),
            "local_time_pdt": (t + dt.timedelta(hours=LOCAL_OFFSET_H)).isoformat(),
            "burst_id": m["burst_id"] if m else "",
            "drop_id_in_burst": m["drop_id"] if m else "",
            # --- delivered node timing, and its quality ---
            "node9_max_cc": f"{cc:.4f}",
            "node4_offset_ms": "" if off is None else f"{off:+.2f}",
            # --- DAS file coverage (coverage is NOT detection) ---
            "nano_covered": m["nano_available"] if m else "0",
            "deep_covered": m["deep_available"] if m else "0",
            "paired_covered": m["paired_available"] if m else "0",
            # --- DAS detection, Nano ---
            "nano_loo_signal_ncc": "" if ncc is None else f"{ncc:.4f}",
            "nano_loo_noise_ncc": "" if n is None else f"{float(n['loo_noise_ncc']):.4f}",
            "nano_beam_snr_db": "" if snr is None else f"{snr:.2f}",
            "nano_relative_amplitude": "" if n is None else f"{float(n['relative_amplitude']):.4f}",
            "nano_detected": "" if ncc is None else int(ncc > NCC_DETECT and snr > SNR_DETECT_DB),
            "source_miniseed": src,
        })

    with (HERE / "awd_drop_catalog.csv").open("w", newline="") as s:
        w = csv.DictWriter(s, fieldnames=list(catalog[0]))
        w.writeheader()
        w.writerows(catalog)

    bursts = {}
    for row in catalog:
        if row["burst_id"] != "":
            bursts.setdefault(int(row["burst_id"]), []).append(row)

    with (HERE / "awd_burst_summary.csv").open("w", newline="") as s:
        w = csv.writer(s)
        w.writerow(["burst_id", "n_drops", "start_utc", "start_local_pdt", "duration_s",
                    "median_node9_cc", "n_nano_covered", "n_deep_covered",
                    "n_paired_covered", "n_nano_detected", "median_nano_beam_snr_db"])
        for b in sorted(bursts):
            rs = bursts[b]
            ts = [dt.datetime.fromisoformat(r["utc_time"]) for r in rs]
            cc = sorted(float(r["node9_max_cc"]) for r in rs)
            snr = sorted(float(r["nano_beam_snr_db"]) for r in rs if r["nano_beam_snr_db"])
            w.writerow([
                b, len(rs), min(ts).isoformat(),
                (min(ts) + dt.timedelta(hours=LOCAL_OFFSET_H)).isoformat(),
                f"{(max(ts) - min(ts)).total_seconds():.1f}",
                f"{cc[len(cc) // 2]:.4f}",
                sum(int(r["nano_covered"] or 0) for r in rs),
                sum(int(r["deep_covered"] or 0) for r in rs),
                sum(int(r["paired_covered"] or 0) for r in rs),
                sum(int(r["nano_detected"] or 0) for r in rs),
                f"{snr[len(snr) // 2]:.2f}" if snr else "",
            ])

    d = np.array([o for o in offsets_ms if o is not None])
    med = float(np.median(d))
    span_h = (max(t9) - min(t9)).total_seconds() / 3600.0
    det = [int(r["nano_detected"]) for r in catalog if r["nano_detected"] != ""]
    lines = [
        "AWD weight-drop timing uncertainty — two independent nodes",
        "=" * 58, "",
        f"node 453009664 (p26.cc9) : {t9.size} picks, median Max_CC {np.median(c9):.4f}",
        f"node 453001432 (p26.cc4) : {t4.size} picks, median Max_CC {np.median(c4):.4f}",
        f"matched within 1 s       : {d.size}", "",
        f"median offset            : {med:+.2f} ms",
        f"MAD                      : {float(np.median(np.abs(d - med))):.2f} ms",
        f"5th / 95th percentile    : {np.percentile(d, 5):+.2f} / {np.percentile(d, 95):+.2f} ms",
    ]
    lines.append(f"(picks quantised to {PICK_QUANTUM_MS} ms; counts use >= so that a")
    lines.append(" threshold landing on a grid value is unambiguous)")
    for th in (2, 10, 50, 100):
        k = int((np.abs(d) >= th).sum())
        lines.append(f"|offset| >= {th:3d} ms       : {k:4d} drops ({100 * k / d.size:.1f}%)")
    lines += [
        "", f"drops with node9 CC < 0.9: {int((c9 < 0.9).sum())}/{c9.size} "
            f"({100 * (c9 < 0.9).mean():.1f}%)",
        "", f"survey span (UTC)        : {min(t9):%Y-%m-%d %H:%M} -> {max(t9):%Y-%m-%d %H:%M}"
            f"  ({span_h:.1f} h)",
        f"survey span (local PDT)  : "
        f"{min(t9) + dt.timedelta(hours=LOCAL_OFFSET_H):%Y-%m-%d %H:%M} -> "
        f"{max(t9) + dt.timedelta(hours=LOCAL_OFFSET_H):%Y-%m-%d %H:%M}",
        "", f"Nano drops detected      : {sum(det)}/{len(det)} "
            f"(LOO NCC > {NCC_DETECT} AND beam SNR > {SNR_DETECT_DB:.0f} dB)",
        "",
        "Both nodes sit a few tens of metres from the source, so both picks carry the",
        "SAME unknown source-to-node travel time. This measures RELATIVE pick quality,",
        "not the absolute offset. The absolute offset is the ~90 ms measured against",
        "the check shot; see faultzone/digitize_checkshot.py.",
    ]
    (HERE / "timing_uncertainty.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwrote awd_drop_catalog.csv ({len(catalog)} rows), "
          f"awd_burst_summary.csv ({len(bursts)} bursts), timing_uncertainty.txt")


if __name__ == "__main__":
    main()
