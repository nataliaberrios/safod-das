"""Stack-size monitoring sensitivity summary for the validated Nano observable."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
TARGET_EARTH_TIDE = 5e-4  # illustrative benchmark; not a measured tide signal


def main():
    d = np.load(HERE / "nano_dvv_injection_recovery.npz", allow_pickle=True)
    injected = d["injected_dvv"].astype(float)
    estimated = d["estimated_dvv"].astype(float)
    null = estimated[np.isclose(injected, 0.)]
    sigma = float(np.std(null, ddof=1))
    stack_sizes = np.array([1, 2, 4, 8, 16, 32, 46])
    # Independence scaling is explicitly a model, not a new empirical claim.
    threshold95 = 1.96 * sigma / np.sqrt(stack_sizes)
    power95 = 3.24 * sigma / np.sqrt(stack_sizes)
    conv = pd.read_csv(HERE / "nano_stack_convergence.csv")
    summary = (conv.groupby("n_drops_per_substack")["independent_substack_ncc"]
               .agg(["median", lambda x: np.percentile(x, 16),
                     lambda x: np.percentile(x, 84)]).reset_index())
    summary.columns = ["n_drops", "median_ncc", "ncc_p16", "ncc_p84"]
    out = pd.DataFrame({"stack_size": stack_sizes,
                        "95pct_null_threshold": threshold95,
                        "approx_95pct_power_threshold": power95})
    out.to_csv(HERE / "monitoring_sensitivity.csv", index=False)
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.4), constrained_layout=True)
    ax[0].loglog(stack_sizes, threshold95, "o-", label="95% null threshold")
    ax[0].loglog(stack_sizes, power95, "s--", label="approx. 95% power")
    ax[0].axhline(TARGET_EARTH_TIDE, color="#d95f02", ls=":",
                  label="illustrative $5\\times10^{-4}$ benchmark")
    ax[0].set(xlabel="number of independent stacks", ylabel="fractional apparent-moveout change",
              title="A  Modelled sensitivity versus stack size")
    ax[0].grid(True, which="both", alpha=.2); ax[0].legend(frameon=False, fontsize=8)
    ax[1].errorbar(summary.n_drops, summary.median_ncc,
                   yerr=[summary.median_ncc-summary.ncc_p16,
                         summary.ncc_p84-summary.median_ncc], fmt="o-",
                   color="#0b5d8b", capsize=3)
    ax[1].set(xlabel="drops per independent within-burst substack",
              ylabel="NCC to independent reference", title="B  Empirical repeatability convergence")
    ax[1].grid(True, alpha=.2)
    fig.suptitle("Nano AWD monitoring sensitivity — observable-level result", fontsize=13)
    fig.savefig(HERE / "monitoring_sensitivity.png", dpi=220)
    (HERE / "monitoring_sensitivity.txt").write_text(
        "Nano AWD monitoring sensitivity\n"
        "===============================\n"
        f"Null trials: {len(null)}; null standard deviation: {sigma:.6g}.\n"
        "Threshold curves use sigma/sqrt(N) independence scaling and are not a replacement for new synthetic injections.\n"
        f"Illustrative Earth-tide benchmark: {TARGET_EARTH_TIDE:.1e}; this is a comparison target, not a detected signal.\n"
        "The threshold applies to the accepted Nano apparent-moveout observable, not directly to formation dV_P/V_P.\n"
    )
    print(f"Saved {HERE/'monitoring_sensitivity.png'}; null sigma={sigma:.6g}")


if __name__ == "__main__":
    main()
