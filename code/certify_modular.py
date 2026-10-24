"""
Rigorous (outward-rounded) check of the sieve inequality of Paper III (holomorphic route):

    N = n1 + n2,  n1, n2 odd sums of two squares,  N = 2 mod 4,  weights r(n1) r(n2),
    level of distribution theta (= 1/2 by Theorem 3.2; 2/5 from Waibel's bound alone),  z = N^alpha,
    Richert weights up to y1 = N^beta1,

    net = L - lam * R  > 0,   in units of X V  (V = prod_{p < z} (1 - g(p))),
      L = f1(sA) F1(sB) + F1(sA) f1(sB) - F1(sA) F1(sB),   sA = thetaA/alpha, sB = (theta - thetaA)/alpha,
      R = 2 int_alpha^beta1 (1 - u/beta1) F1(t/alpha) F1((theta - u - t)/alpha) du/u   (any split t = t(u)).
    Then  Omega(n1 n2) < 1/lam + 2/beta1  and  max(Omega(n1), Omega(n2)) < 1/lam + 1/beta1.

This certifies the numerical inequality only (not the sieve lemmas, not the o(1) terms).
Linear sieve functions (kappa = 1, beta = 2): F = 2e^g/s on (0,3], f = 2e^g log(s-1)/s on [2,4],
(sF)' = f(s-1), (sf)' = F(s-1); F decreasing, f increasing. Closed forms in mpmath interval arithmetic
(80 bits); the recurrence in IEEE doubles with one outward nextafter step after each operation
(same scheme as bilateral-two-squares-sieve/code/certify.py, from which table() is taken).
"""
import math
import sys
from fractions import Fraction
import numpy as np
from mpmath import iv
iv.prec = 80

dn = lambda x: math.nextafter(x, -math.inf)
up = lambda x: math.nextafter(x, math.inf)
EG_iv = iv.exp(iv.euler)
H_DEN = 1000


def table(s_top):
    """Rigorous enclosures of F1, f1 on the grid i/H_DEN, 0 <= i <= s_top*H_DEN (kappa = 1)."""
    n = s_top * H_DEN
    A = 2 * EG_iv
    Fc = lambda s: A / s
    fc = lambda s: A * iv.log(s - 1) / s
    Flo = np.zeros(n + 1); Fhi = np.full(n + 1, math.inf); flo = np.zeros(n + 1); fhi = np.zeros(n + 1)
    pl = [0.0] * (n + 1); ph = [0.0] * (n + 1)
    for i in range(n + 1):
        v = iv.mpf([i, i]) / H_DEN; pl[i], ph[i] = dn(float(v.a)), up(float(v.b))
    iF, if0, if1 = 3 * H_DEN, 2 * H_DEN, 4 * H_DEN
    for i in range(1, if1 + 1):
        s = iv.mpf([i, i]) / H_DEN
        if i <= iF:
            v = Fc(s); Flo[i], Fhi[i] = dn(float(v.a)), up(float(v.b))
        if if0 < i <= if1:
            v = fc(s); flo[i], fhi[i] = max(0.0, dn(float(v.a))), up(float(v.b))
    gF_lo, gF_hi = dn(pl[iF] * Flo[iF]), up(ph[iF] * Fhi[iF])
    gf_lo, gf_hi = dn(pl[if1] * flo[if1]), up(ph[if1] * fhi[if1])
    for i in range(iF, n):
        j = i - H_DEN
        w_lo, w_hi = dn(pl[i + 1] - ph[i]), up(ph[i + 1] - pl[i])
        gF_lo = dn(gF_lo + dn(flo[j] * w_lo)); gF_hi = up(gF_hi + up(fhi[j + 1] * w_hi))
        Flo[i + 1] = dn(gF_lo / ph[i + 1]); Fhi[i + 1] = up(gF_hi / pl[i + 1])
        if i >= if1:
            gf_lo = dn(gf_lo + dn(Flo[j + 1] * w_lo)); gf_hi = up(gf_hi + up(Fhi[j] * w_hi))
            flo[i + 1] = dn(gf_lo / ph[i + 1]); fhi[i + 1] = up(gf_hi / pl[i + 1])
    return dict(Flo=Flo, Fhi=Fhi, flo=flo, fhi=fhi, n=n, A_hi=up(float(A.b)))


def F_up(T, s):
    if s <= 3 - 1e-9:
        return up(T["A_hi"] / dn(s) * (1 + 1e-15))
    i = int(math.floor(s * H_DEN)) - 1
    return T["Fhi"][min(i, T["n"])]


def F_dn(T, s):
    i = int(math.ceil(s * H_DEN)) + 1
    return T["Flo"][i] if i <= T["n"] else 0.0


def f_dn(T, s):
    i = int(math.floor(s * H_DEN)) - 1
    return T["flo"][min(i, T["n"])] if i >= 0 else 0.0


lin = table(60)
print(f"linear sieve: F(3) in [{lin['Flo'][3000]:.6f},{lin['Fhi'][3000]:.6f}], f(4) in [{lin['flo'][4000]:.6f},"
      f"{lin['fhi'][4000]:.6f}], F(60) <= {lin['Fhi'][-1]:.6f}, f(60) >= {lin['flo'][-1]:.6f}")


def certify(theta_s, alpha_s, lam_s, beta1_s, label, NU=4000):
    theta, alpha, lam, beta1 = (float(Fraction(x)) for x in (theta_s, alpha_s, lam_s, beta1_s))
    bound = 1 / Fraction(lam_s) + 2 / Fraction(beta1_s)
    k = int(-(-bound // 1)) - 1
    side = 1 / Fraction(lam_s) + 1 / Fraction(beta1_s)
    ks = min(int(-(-side // 1)) - 1, k - 1)
    print(f"\n=== {label}: theta={theta_s}, alpha={alpha_s}, lam={lam_s}, beta1={beta1_s};"
          f" 1/lam + 2/beta1 = {bound} => Omega(n1 n2) <= {k};"
          f" 1/lam + 1/beta1 = {side} ~ {float(side):.4f} => max Omega(n_i) <= {ks} ===")
    assert beta1 < theta, "Richert weights must stop below the level"

    # ---- L (lower bound); the split thetaA is a free parameter, take the best grid point ----
    best = None
    for tA in np.linspace(0.3 * theta, 0.7 * theta, 801):
        sA, sB = tA / alpha, (theta - tA) / alpha
        v = f_dn(lin, sA) * F_dn(lin, sB) + F_up(lin, sA) * (f_dn(lin, sB) - F_up(lin, sB))
        if best is None or v > best[0]:
            best = (v, tA)
    tA = best[1]; sA, sB = tA / alpha, (theta - tA) / alpha
    L_lo = dn(dn(f_dn(lin, sA) * F_dn(lin, sB)) + dn(F_up(lin, sA) * dn(f_dn(lin, sB) - F_up(lin, sB))))

    # ---- R (upper bound): on [ua, ub] the factor (1 - u/beta1)/u is largest at ua and the upper sieve
    #      bound at level theta - u is largest at ub; any split t gives an upper bound for the min over t ----
    R_hi = 0.0
    for kk in range(NU):
        ua, ub = dn(alpha + (beta1 - alpha) * kk / NU), up(alpha + (beta1 - alpha) * (kk + 1) / NU)
        lev = theta - ub
        bestU = min(F_up(lin, t / alpha) * F_up(lin, dn((lev - t) / alpha))
                    for t in np.linspace(lev * 0.3, lev * 0.7, 81))
        weight = up(up(1 - ua / beta1) / ua)
        R_hi = up(R_hi + up(up(weight * bestU) * up(ub - ua)))
    R_hi = up(2 * R_hi)

    net_lo = dn(L_lo - up(up(lam * (1 + 2.0 ** -50)) * R_hi))
    print(f"thetaA = {tA:.4f} (sA = {sA:.4f}, sB = {sB:.4f})")
    print(f"L   >= {L_lo:.6f}")
    print(f"R   <= {R_hi:.6f}")
    print(f"certified:  net = L - lam R >= {net_lo:+.6f}   ->  {'POSITIVE' if net_lo > 0 else 'NOT certified'}")
    return net_lo


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "1/2":
        # level 1/2 (Theorem 3.2, using the norm bound of Part II): Table 2, first line; Theorem 1.1
        certify("1/2", "6/125", "1/3", "17/50", "level 1/2: joint 8 and max side 5")
    else:
        # level 2/5 (what Waibel's norm bound alone gives, Remark 4.5): Table 2, lines 2-4
        certify("2/5", "9/250", "3/10", "3/10", "level 2/5: joint 9 and max side 6")
        certify("2/5", "4/125", "7/30", "7/20", "level 2/5: joint 9 (side 7), more margin")
        certify("2/5", "9/250", "1/3", "3/10", "level 2/5: joint 9 and max side 6, lambda = 1/3")
