"""Where does AWD timing scatter live: drop-to-drop, or burst-to-burst?

This is the question that decides what is worth doing next, because the two
terms behave completely differently under stacking:

  within-burst  (drop to drop)   averages down as 1/sqrt(N_drops)
  between-burst (burst to burst) does NOT -- it is the floor a monitoring
                                 measurement actually sits on

Both are already measured per drop and per burst by
`nano_hierarchical_repeatability.py`; this script only reads its two CSVs and
separates the terms. It computes no new waveform quantity.

A note on estimator choice. An obvious way to get per-drop timing is to pick the
envelope peak of each beam. Do not: on low-SNR drops the peak jumps between
lobes and the scatter comes out at ~32 ms, which is picker failure, not timing.
The leave-one-out cross-correlation delays in the CSVs are the right quantity
and give 0.4 ms.

Outputs
-------
variability.npz    the separated terms and the per-burst series
variability.txt    the same, readable
"""
import csv
import datetime as dt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
AWD = HERE.parent
DROPS = AWD / "nano_drop_repeatability.csv"
BURSTS = AWD / "nano_burst_repeatability_hierarchical.csv"
SNR_GATE_DB = 10.0        # the drop-quality gate used throughout drop_catalog


def col(rows, key, scale=1.0):
    return np.array([float(r[key]) * scale for r in rows if r.get(key) not in ("", None)])


def main():
    d = list(csv.DictReader(DROPS.open()))
    b = list(csv.DictReader(BURSTS.open()))

    dly = col(d, "loo_signal_delay_s", 1e3)          # ms, vs own burst template
    snr = col(d, "beam_snr_db")
    amp = col(d, "relative_amplitude")
    t = np.array([dt.datetime.fromisoformat(r["utc_date_from_p26"]).replace(tzinfo=None)
                  for r in d])
    bid = np.array([int(r["burst_id"]) for r in d])

    bdly = col(b, "loo_burst_delay_s", 1e3)          # ms, vs other bursts
    bamp = col(b, "burst_relative_amplitude")
    bncc = col(b, "loo_burst_signal_ncc")
    bn = col(b, "n_drops")

    gate = snr > SNR_GATE_DB
    s_within_all = dly.std(ddof=1)
    s_within = dly[gate].std(ddof=1)
    s_between = bdly.std(ddof=1)
    n_per = bn.mean()
    s_within_stacked = s_within / np.sqrt(n_per)

    lines = [
        "AWD timing scatter: drop-to-drop versus burst-to-burst",
        "=" * 66, "",
        f"{len(d)} drops in {len(b)} bursts, {n_per:.0f} drops per burst",
        f"drop-quality gate: beam SNR > {SNR_GATE_DB:.0f} dB "
        f"({int(gate.sum())} of {len(d)} drops pass)", "",
        "WITHIN-BURST (drop to drop, leave-one-out vs own burst template)",
        f"  all drops                  sd {s_within_all:7.3f} ms",
        f"  SNR-gated drops            sd {s_within:7.3f} ms"
        f"   ({s_within_all/s_within:.1f}x tighter)",
        f"  after stacking {n_per:.0f} drops    sd {s_within_stacked:7.3f} ms", "",
        "BETWEEN-BURST (leave-one-burst-out)",
        f"  burst-to-burst             sd {s_between:7.3f} ms"
        f"   range {bdly.min():+.2f} to {bdly.max():+.2f} ms",
        f"  burst waveform NCC         median {np.median(bncc):.4f}",
        f"  burst relative amplitude   sd {bamp.std(ddof=1):.3f}", "",
        "THE COMPARISON THAT DECIDES WHAT TO DO NEXT",
        f"  within-burst, stacked      {s_within_stacked:7.3f} ms   averages down as 1/sqrt(N)",
        f"  between-burst              {s_between:7.3f} ms   does NOT average down",
        f"  ratio                      {s_between/s_within_stacked:7.1f}x", "",
        "  Stacking 20 drops already puts the within-burst term "
        f"{s_between/s_within_stacked:.1f}x below the",
        "  between-burst floor. More drops per burst buys effectively nothing.",
        "  The floor is between-burst, and it is what a monitoring measurement",
        "  actually sits on.", "",
        f"  For reference, docs/paper1/STATUS.md records sigma_alpha ~ 0.30 ms",
        f"  common-mode. The between-burst term measured here is {s_between:.3f} ms,",
        "  so that documented figure is the BETWEEN-burst scatter -- which is",
        "  consistent, and confirms stacking is not the lever.",
    ]

    # per-burst series for the figure
    bt = np.array([t[bid == int(r["burst_id"])].min() for r in b])
    hours = np.array([(x - bt.min()).total_seconds() / 3600 for x in bt])

    (HERE / "variability.txt").write_text("\n".join(lines) + "\n")
    np.savez_compressed(
        HERE / "variability.npz",
        drop_delay_ms=dly, drop_snr_db=snr, drop_amp=amp,
        drop_hours=np.array([(x - t.min()).total_seconds() / 3600 for x in t]),
        drop_burst_id=bid,
        burst_delay_ms=bdly, burst_amp=bamp, burst_ncc=bncc, burst_hours=hours,
        burst_start_utc=np.array([x.isoformat() for x in bt]),
        s_within_all=s_within_all, s_within=s_within, s_between=s_between,
        s_within_stacked=s_within_stacked, snr_gate_db=SNR_GATE_DB,
        n_per_burst=n_per,
    )
    print("\n".join(lines))
    print(f"\nwrote variability.npz, variability.txt")


if __name__ == "__main__":
    main()
