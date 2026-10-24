"""
Figure: the sieve inequality L - lam R > 0 in the (lam, beta1) plane, for the level theta = 1/2 (alpha = 6/125;
Theorem 3.2, giving Theorem 1.1) and theta = 2/5 (alpha = 9/250; what Waibel's norm bound alone gives,
Remark 4.5), together with the curves 1/lam + 2/beta1 = k + 1 (Omega(n(N-n)) <= k) and
1/lam + 1/beta1 = k' + 1 (max Omega <= k').

Floating-point sieve functions of sieve_core.py (same model as modular_sieve_scan.py); the four points are the
parameter sets of Table 2 of the paper, certified by certify_modular.py.
Usage:  python plot_sieve_landscape.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from modular_sieve_scan import best_L, R_value

panels = [
    dict(theta=0.4, alpha=9 / 250, title=r"$\theta=2/5$ (Waibel's bound alone, Remark 4.5), $\alpha=9/250$",
         points=[(3 / 10, 3 / 10, "Table 2, line 2", (-62, -12)), (7 / 30, 7 / 20, "line 3", (6, 4)),
                 (1 / 3, 3 / 10, "line 4", (6, -12))]),
    dict(theta=0.5, alpha=6 / 125, title=r"$\theta=1/2$ (Theorem 3.2), $\alpha=6/125$",
         points=[(1 / 3, 17 / 50, "Table 2, line 1", (6, -12))]),
]
lam = np.linspace(0.15, 0.5, 141)

fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), sharey=False)
for ax, P in zip(axes, panels):
    theta, alpha = P["theta"], P["alpha"]
    L, tA = best_L(alpha, theta)
    beta1 = np.linspace(0.20, theta - 0.005, 90)
    R = np.array([R_value(alpha, theta, b, nu=200) for b in beta1])
    net = L - lam[None, :] * R[:, None]
    cf = ax.contourf(lam, beta1, net, levels=np.linspace(-1.2, 0.6, 19), cmap="RdBu", extend="both", zorder=-5)
    ax.set_rasterization_zorder(-1)  # rasterise the filled contours: no white seams between levels in PDF viewers
    ax.contour(lam, beta1, net, levels=[0.0], colors="k", linewidths=1.6)
    # constraint curves: beta1 = 2/(c - 1/lam) and beta1 = 1/(c - 1/lam); the bound holds ABOVE the curve
    LAM, B1 = np.meshgrid(lam, beta1)
    for c, k, ls in ((10, 9, "-"), (9, 8, "--")):
        b = np.where(c - 1 / lam > 0, 2 / np.maximum(c - 1 / lam, 1e-9), np.nan)
        ax.plot(lam, b, color="#444", ls=ls, lw=1.1, label=rf"$1/\lambda+2/\beta_1={c}$  ($\Omega(n(N-n))\leq {k}$)")
    for c, k, ls in ((7, 6, "-"), (6, 5, "--")):
        b = np.where(c - 1 / lam > 0, 1 / np.maximum(c - 1 / lam, 1e-9), np.nan)
        ax.plot(lam, b, color="#8e5bd6", ls=ls, lw=1.1, label=rf"$1/\lambda+1/\beta_1={c}$  ($\max\Omega\leq {k}$)")
    # hatched: parameters giving (9, 6) resp. (8, 5) with L - lam R > 0
    feas96 = (net > 0) & (1 / LAM + 2 / B1 <= 10) & (1 / LAM + 1 / B1 <= 7)
    feas85 = (net > 0) & (1 / LAM + 2 / B1 <= 9) & (1 / LAM + 1 / B1 <= 6)
    for mask, hatch, lab in ((feas96, "///", "(9, 6) admissible"), (feas85, "xxx", "(8, 5) admissible")):
        if mask.any():
            ax.contourf(lam, beta1, mask.astype(float), levels=[0.5, 1.5], colors="none", hatches=[hatch])
            ax.plot([], [], color="k", lw=0, label=lab + f"  [{hatch}]")
    for (l0, b0, name, off) in P["points"]:
        ax.plot(l0, b0, "o", color="#eb6834", ms=6, mec="k")
        ax.annotate(name, (l0, b0), xytext=off, textcoords="offset points", fontsize=7,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    ax.set_xlim(lam[0], lam[-1]); ax.set_ylim(beta1[0], beta1[-1])
    ax.set_xlabel(r"$\lambda$"); ax.set_ylabel(r"$\beta_1$")
    ax.set_title(P["title"] + rf",  $L={L:.4f}$", fontsize=9)
    ax.grid(True, alpha=0.2)
    print(f"theta={theta}: L={L:.5f} at thetaA={tA:.4f}; R(beta1) from {R[0]:.3f} to {R[-1]:.3f}")
cb = fig.colorbar(cf, ax=axes, fraction=0.03, pad=0.02)
cb.set_label(r"$L-\lambda R$  (black: $L-\lambda R=0$)")
h0, l0 = axes[0].get_legend_handles_labels()
h1, l1 = axes[1].get_legend_handles_labels()
for h, l in zip(h1, l1):
    if l not in l0:
        h0.append(h); l0.append(l)
fig.legend(h0, l0, fontsize=7, ncol=3, loc="lower center", bbox_to_anchor=(0.45, -0.1), frameon=False)
fig.savefig("../figures/sieve_landscape.pdf", bbox_inches="tight", dpi=300)
fig.savefig("../figures/sieve_landscape.png", dpi=150, bbox_inches="tight")
print("wrote ../figures/sieve_landscape.pdf/.png")
