"""
Floating-point parameter scan for the sieve of Paper III (holomorphic route).

Sequence: n odd, n <= N, N = 2 mod 4, weights w(n) = r(n) r(N-n); level of distribution theta = 2/5 (Theorem 2).
Sift all odd primes < z = N^alpha on both sides: vector sieve of two linear components (side A: p | n, side B:
p | N-n; each has dimension 1 since sum_{split p<z} h(p) ~ log log z, h(p) = (2p-1)/(p(p+1))).

    L   = f1(sA) F1(sB) + F1(sA) f1(sB) - F1(sA) F1(sB),   sX = thetaX/alpha,  thetaA + thetaB = theta
    R   = 2 int_alpha^beta1 (1 - u/beta1) U(theta - u) du/u,   U(t) = min_{tA+tB=t} F1(tA/alpha) F1(tB/alpha)
    positivity L - lam R > 0  =>  Omega(n(N-n)) < 1/lam + 2/beta1,  max(Omega(n), Omega(N-n)) < 1/lam + 1/beta1.

Heuristic main-term model, NOT certified (certify_modular.py does that for the chosen parameters).
Uses the linear sieve functions of ../bilateral-two-squares-sieve/code/sieve_core.py.
"""
import sys
import os
import numpy as np

from sieve_core import F_lin, f_lin

THETA = 0.4


def L_lower(alpha, theta, tA):
    sA, sB = tA / alpha, (theta - tA) / alpha
    return float(f_lin(sA) * F_lin(sB) + F_lin(sA) * f_lin(sB) - F_lin(sA) * F_lin(sB))


def U_upper(alpha, t):
    ts = np.linspace(0.02 * t, 0.98 * t, 97)
    return float(np.min(F_lin(ts / alpha) * F_lin((t - ts) / alpha)))


def R_value(alpha, theta, beta1, nu=400):
    u = np.linspace(alpha, beta1, nu + 1)
    U = np.array([U_upper(alpha, theta - x) for x in u])
    return 2 * float(np.trapz((1 - u / beta1) * U / u, u))


def best_L(alpha, theta):
    best = (-np.inf, theta / 2)
    for tA in np.linspace(0.3 * theta, 0.7 * theta, 81):
        v = L_lower(alpha, theta, tA)
        if v > best[0]:
            best = (v, tA)
    return best


if __name__ == "__main__":
    theta = float(sys.argv[1]) if len(sys.argv) > 1 else THETA
    print(f"level theta = {theta}")
    print("no weights: Omega(n(N-n)) < 2/alpha requires L > 0")
    for alpha in (0.04, 0.05, 0.06, 0.07, 0.075, 0.08, 0.09, 0.10):
        L, tA = best_L(alpha, theta)
        print(f"  alpha={alpha:.3f}: L={L:+.4f} (thetaA={tA:.4f}, s={tA/alpha:.3f}), 2/alpha={2/alpha:.1f}")
    print()
    print("Richert weights: targets (joint k, side k') from 1/lam + 2/beta1 <= k+1 and 1/lam + 1/beta1 <= k'+1")
    rows = []
    for alpha in np.arange(0.010, 0.0801, 0.002):
        L, tA = best_L(alpha, theta)
        if L <= 0:
            continue
        for beta1 in np.arange(0.22, 0.401, 0.01):
            R = R_value(alpha, theta, beta1, nu=200)
            lam_max = L / R                       # largest admissible lambda (net -> 0)
            joint = 1 / lam_max + 2 / beta1       # smallest achievable 1/lam + 2/beta1 (not integer-adjusted)
            side = 1 / lam_max + 1 / beta1
            rows.append((joint, side, alpha, beta1, L, R, lam_max, tA))
    rows.sort()
    print("best (joint) configurations:")
    for joint, side, alpha, beta1, L, R, lam, tA in rows[:12]:
        print(f"  1/lam+2/beta1={joint:6.3f} 1/lam+1/beta1={side:6.3f}  alpha={alpha:.3f} beta1={beta1:.2f} "
              f"L={L:.4f} R={R:.4f} lam_max={lam:.4f} thetaA={tA:.4f}")
    rows.sort(key=lambda r: r[1])
    print("best (max side) configurations:")
    for joint, side, alpha, beta1, L, R, lam, tA in rows[:12]:
        print(f"  1/lam+2/beta1={joint:6.3f} 1/lam+1/beta1={side:6.3f}  alpha={alpha:.3f} beta1={beta1:.2f} "
              f"L={L:.4f} R={R:.4f} lam_max={lam:.4f} thetaA={tA:.4f}")
    # concrete candidates with rational lam, beta1 and 5% margin
    print()
    print("candidates with exact (lam, beta1) and margin L - lam R:")
    for (lam, beta1, label) in ((0.2, 0.4, "joint 9 (1/lam+2/beta1 = 10)"),
                                (1 / 6, 1 / 3, "joint 9 / side 6 (1/lam+2/beta1 = 12?)"),
                                (0.25, 1 / 3, "1/lam+2/beta1 = 10, side 7"),
                                (0.2, 0.5, "joint 8 / side 6 (9, 7)"),
                                (0.25, 0.5, "joint 7 / side 5 (8, 6)")):
        best = None
        for alpha in np.arange(0.010, 0.0801, 0.002):
            L, tA = best_L(alpha, theta)
            if L <= 0:
                continue
            R = R_value(alpha, theta, beta1, nu=200)
            net = L - lam * R
            if best is None or net > best[0]:
                best = (net, alpha, L, R, tA)
        net, alpha, L, R, tA = best
        print(f"  lam={lam:.4f} beta1={beta1:.4f} [{label}]: best alpha={alpha:.3f} thetaA={tA:.4f} "
              f"L={L:.4f} R={R:.4f} net={net:+.4f}  -> Omega(n(N-n)) < {1/lam + 2/beta1:.3f}, side < {1/lam + 1/beta1:.3f}")

