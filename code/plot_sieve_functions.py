"""
Figure: the linear-sieve functions F_1, f_1 as enclosed by certify_modular.table(60) (closed forms in mpmath
interval arithmetic on (0,3] and [2,4], delay-differential recursion with outward rounding beyond), the
arguments s_A, s_B at which L is evaluated, the range of arguments entering R, and the width of the enclosures.
Usage:  python plot_sieve_functions.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from certify_modular import lin, H_DEN

n = lin["n"]
s = np.arange(n + 1) / H_DEN
Flo, Fhi, flo, fhi = lin["Flo"], lin["Fhi"], lin["flo"], lin["fhi"]
# the parameters of the first line of Table 2 of the paper
theta, alpha, beta1, tA = 0.4, 9 / 250, 0.3, 0.1998
sA, sB = tA / alpha, (theta - tA) / alpha
# arguments t'/alpha and (theta-u-t')/alpha used in R: u in [alpha, beta1], t' in [0.3, 0.7] (theta-u)
s_R_min, s_R_max = 0.3 * (theta - beta1) / alpha, 0.7 * (theta - alpha) / alpha

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 5.6), height_ratios=[3, 2])
mask = (s >= 0.8) & (s <= 10)
ax1.plot(s[mask], Fhi[mask], color="#2a78d6", lw=1.4, label=r"$F_1(s)$ (upper enclosure)")
mf = mask & (s >= 2)
ax1.plot(s[mf], flo[mf], color="#eb6834", lw=1.4, label=r"$f_1(s)$ (lower enclosure)")
ax1.axhline(1, color="k", lw=0.6, ls=":")
ax1.axvspan(s_R_min, s_R_max, color="#1baf7a", alpha=0.12, label=r"arguments entering $R$")
for x, lab in ((sA, r"$s_A$"), (sB, r"$s_B$")):
    ax1.axvline(x, color="#8e5bd6", lw=0.9, ls="--")
ax1.annotate(r"$s_A\approx5.55,\ s_B\approx5.56$", (sB, 1.25), xytext=(8, 0), textcoords="offset points",
             fontsize=8, color="#8e5bd6")
for x in (3, 4):
    ax1.axvline(x, color="#999", lw=0.6)
ax1.text(3.02, 2.3, "closed form / recursion", fontsize=7, color="#666")
ax1.set_xlim(0.8, 10); ax1.set_ylim(0.6, 2.6)
ax1.set_ylabel(r"$F_1(s),\ f_1(s)$")
ax1.grid(True, alpha=0.25); ax1.legend(fontsize=8, loc="upper right")

wF = np.where(Fhi[1:] < np.inf, Fhi[1:] - Flo[1:], np.nan)
wf = fhi[2 * H_DEN:] - flo[2 * H_DEN:]
ax2.semilogy(s[1:], np.maximum(wF, 1e-19), color="#2a78d6", lw=1, label=r"$F_1^+-F_1^-$")
ax2.semilogy(s[2 * H_DEN:], np.maximum(wf, 1e-19), color="#eb6834", lw=1, label=r"$f_1^+-f_1^-$")
ax2.set_xlim(0, 60); ax2.set_ylim(1e-17, 1e-2)
ax2.text(4.5, 3e-16, "closed forms (80-bit intervals)", fontsize=7, color="#666")
ax2.text(20, 3e-6, "delay-differential recursion, outward rounding", fontsize=7, color="#666")
ax2.set_xlabel(r"$s$"); ax2.set_ylabel("width of enclosure")
ax2.grid(True, which="both", alpha=0.25); ax2.legend(fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig("../figures/sieve_functions.pdf")
fig.savefig("../figures/sieve_functions.png", dpi=150)
print(f"F_1(s_A) <= {Fhi[int(sA*H_DEN)]:.6f}, f_1(s_A) >= {flo[int(sA*H_DEN)]:.6f}; "
      f"max width F: {np.nanmax(wF):.2e}, f: {np.nanmax(wf):.2e}; wrote ../figures/sieve_functions.pdf/.png")
