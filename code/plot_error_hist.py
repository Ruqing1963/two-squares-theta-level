"""
Figure: distribution of the normalised cuspidal error (T - M)/sqrt(N) over admissible N in (X/2, X],
standardised by its root mean square, for several pairs (m, m'), against the standard normal density.

Reuses the tables and the FFT convolution of sup_error_scan.py (same admissibility: (N, mm') = 1, v_2(N) <= 3).
Floating point, not certified; this is the statistical evidence of Section 5 only.

Usage:  python plot_error_hist.py [X] [m,m' ...]        (default X = 2^21, pairs (1,101) (13,29) (1,1009) (61,67))
For mm' ~ 10^6 (e.g. (1009,1013)) the main term M is of size 30 while T is a multiple of 8, so the distribution
is dominated by this discreteness (kurtosis ~ 24); such pairs are deliberately not in the default list.
"""
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sup_error_scan import r_table, sigma_over_n_odd, delta2_table, T_all, M_all, admissible

args = sys.argv[1:]
X = int(args[0]) if args and ',' not in args[0] else 2 ** 21
pairs = [tuple(int(t) for t in a.split(',')) for a in args if ',' in a] or [(1, 101), (13, 29), (1, 1009), (61, 67)]

t0 = time.time()
r = r_table(X)
f_sigma, v2 = sigma_over_n_odd(X)
print(f"X = 2^{X.bit_length() - 1}; tables in {time.time() - t0:.1f}s")

colors = ["#2a78d6", "#eb6834", "#1baf7a", "#8e5bd6", "#c9a400"]
edges = np.linspace(-6, 6, 97)
mid = 0.5 * (edges[1:] + edges[:-1])

fig, ax = plt.subplots(figsize=(6.4, 4.0))
for (m, mp), col in zip(pairs, colors):
    dens2 = delta2_table(m, mp)
    T = T_all(r, m, mp, X)
    M = M_all(m, mp, X, f_sigma, v2, dens2)
    ok = admissible(X, m, mp, v2)
    n = np.arange(X + 1)[ok]
    e = (T[ok] - M[ok]) / np.sqrt(n)
    rms = np.sqrt(np.mean(e ** 2))
    z = e / rms
    hist, _ = np.histogram(z, bins=edges, density=True)
    kurt = np.mean(z ** 4)
    frac3 = np.mean(np.abs(z) > 3)
    print(f"({m},{mp}): n={len(e)}, rms={rms:.3f}, sup/rms={np.max(np.abs(z)):.2f}, "
          f"kurtosis={kurt:.2f}, P(|z|>3)={frac3:.4f} (normal: 0.0027)")
    ax.step(mid, hist, where="mid", color=col, lw=1.2, label=rf"$(m,m')=({m},{mp})$,  rms $={rms:.3f}$")
ax.plot(mid, np.exp(-mid ** 2 / 2) / np.sqrt(2 * np.pi), "k--", lw=1, label="standard normal")
ax.set_yscale("log")
ax.set_ylim(1e-5, 1)
ax.set_xlim(-6, 6)
ax.set_xlabel(r"$(T-M)/\sqrt{N}$ divided by its root mean square,  $N\in(2^{20},2^{21}]$")
ax.set_ylabel("density")
ax.grid(True, which="both", alpha=0.25)
ax.legend(fontsize=8, loc="upper right")
fig.tight_layout()
fig.savefig("../figures/error_hist.pdf")
fig.savefig("../figures/error_hist.png", dpi=150)
print("wrote ../figures/error_hist.pdf/.png")
