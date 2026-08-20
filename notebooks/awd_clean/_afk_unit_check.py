"""Unit checks for adaptive_fk before it is used on real data."""
import sys, numpy as np
sys.path.insert(0, ".")
from ambient_lellouch2019_exact_stack import adaptive_fk

rng = np.random.default_rng(7)
x = rng.standard_normal((64, 600)).astype(np.float32)

# 1. alpha = 0 must be a bit-identical no-op (the option must be inert by default)
y = adaptive_fk(x, 250.0, 1.0, 0.0, 32, 200, False)
assert y is x or np.array_equal(x, y), "alpha=0 is not a no-op"
print("1. alpha=0 no-op                      : bit-identical  OK")

# 2. Bartlett 50%-overlap weights must sum to unity everywhere they cover, so a
#    unity mask reconstructs the input. Test with alpha tiny -> mask ~ const.
z = adaptive_fk(x, 250.0, 1.0, 1e-12, 32, 200, False)
rel = np.abs(z - x).max() / np.abs(x).max()
print("2. overlap-add reconstruction         : max rel err %.3e %s"
      % (rel, "OK" if rel < 1e-6 else "FAIL"))

# 3. It must actually enhance a coherent plane wave against noise.
nch, nt, fs, dx, v, f0 = 64, 600, 250.0, 1.0, 3200.0, 12.0
t = np.arange(nt) / fs
ch = np.arange(nch) * dx
wave = np.sin(2 * np.pi * f0 * (t[None, :] - ch[:, None] / v)).astype(np.float32)
noise = rng.standard_normal((nch, nt)).astype(np.float32) * 3.0
mixed = wave + noise
def snr(d):
    num = float(np.sum(d * wave)) ** 2
    den = float(np.sum(d * d)) * float(np.sum(wave * wave))
    return num / den if den > 0 else 0.0
for a in (0.5, 1.0, 2.0):
    for norm in (False, True):
        out = adaptive_fk(mixed, fs, dx, a, 32, 200, norm)
        print("3. alpha=%.1f %-5s: coherence with truth %.4f -> %.4f"
              % (a, "NAFK" if norm else "AFK", snr(mixed), snr(out)))

# 4. It must NOT create moveout from noise alone (the recurring failure here).
pure = rng.standard_normal((nch, nt)).astype(np.float32)
for a in (1.0, 2.0):
    out = adaptive_fk(pure, fs, dx, a, 32, 200, False)
    print("4. alpha=%.1f on pure noise: coherence with a 3200 m/s wave %.4f "
          "(input %.4f)" % (a, snr(out), snr(pure)))
