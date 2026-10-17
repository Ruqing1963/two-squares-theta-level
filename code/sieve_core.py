"""
Floating-point beta-sieve functions and the switching density c(delta, alpha).
(Adapted from almost-prime-sum-of-two-squares/code/exploration; not interval-certified.)

  kappa = 1/2, beta = 1 (semi-linear):  F(s) = A/sqrt(s) on (0,2],  f(s) = A arccosh(sqrt s)/sqrt s on [1,3]
  kappa = 1,   beta = 2 (linear):       F(s) = 2e^g/s  on (0,3],    f(s) = 2e^g log(s-1)/s    on [2,4]
  (s^k F)' = k s^(k-1) f(s-1),  (s^k f)' = k s^(k-1) F(s-1)
"""
import numpy as np

EG = np.exp(np.euler_gamma)
H = 1e-4
S_MAX = 30.0
_grid = np.arange(0, S_MAX + H / 2, H)
_n1 = int(round(1 / H))


def _solve(kappa, beta, A):
    F = np.zeros_like(_grid); f = np.zeros_like(_grid)
    s = _grid
    i0F = int(round((beta + 1) / H)); i0f = int(round(beta / H))
    F[1:i0F + 1] = A * s[1:i0F + 1] ** -kappa
    gF = s ** kappa * F; gf = np.zeros_like(s)
    if kappa == 0.5:   # exact on [1,3]; avoids the integrable singularity of F(s-1) at s = 1
        j = int(round(3 / H))
        gf[i0f:j + 1] = A * np.arccosh(np.sqrt(s[i0f:j + 1])); f[i0f:j + 1] = gf[i0f:j + 1] / s[i0f:j + 1] ** kappa
        i0f = j
    for i in range(min(i0f, i0F) + 1, len(s)):
        s0, s1 = s[i - 1], s[i]
        if i > i0F:
            gF[i] = gF[i - 1] + H * kappa * (s0 ** (kappa - 1) * f[i - 1 - _n1] + s1 ** (kappa - 1) * f[i - _n1]) / 2
            F[i] = gF[i] / s1 ** kappa
        if i > i0f:
            gf[i] = gf[i - 1] + H * kappa * (s0 ** (kappa - 1) * F[i - 1 - _n1] + s1 ** (kappa - 1) * F[i - _n1]) / 2
            f[i] = gf[i] / s1 ** kappa
    return F, f


_A_HALF = 2 * np.sqrt(EG / np.pi)
_A_ONE = 2 * EG
_Fh, _fh = _solve(0.5, 1, _A_HALF)
_F1, _f1 = _solve(1.0, 2, _A_ONE)


def _at(arr, s):
    return np.interp(np.minimum(s, S_MAX), _grid, arr)


def F_half(s):
    s = np.asarray(s, dtype=float)
    return np.where(s <= 2, _A_HALF / np.sqrt(np.maximum(s, 1e-12)), _at(_Fh, s))


def f_half(s):
    return _at(_fh, np.asarray(s, dtype=float))


def F_lin(s):
    s = np.asarray(s, dtype=float)
    return np.where(s <= 3, _A_ONE / np.maximum(s, 1e-12), _at(_F1, s))


def f_lin(s):
    return _at(_f1, np.asarray(s, dtype=float))


def c_switch(delta, alpha, h=5e-4):
    """|T|/(N/log N): triples q1 < q2 inert, q1 >= N^(1/2-delta), q1 q2 m <= N, m made of split primes >= N^alpha."""
    mu = np.arange(0, 2 * delta + h / 2, h)
    nu1 = np.where(mu >= alpha, 0.5 / np.maximum(mu, 1e-12), 0.0) * h
    total = np.zeros_like(mu); total[0] = 1.0
    term = total.copy()
    for j in range(1, int(2 * delta / alpha) + 1):
        term = np.convolve(term, nu1)[: len(mu)] / j
        total += term
    lo, hi = 0.5 - delta, (1 - mu) / 2
    I = np.where(hi > lo, (np.log(hi / lo) - np.log((1 - mu - hi) / (1 - mu - lo))) / (1 - mu), 0.0)
    return float(np.sum(total * 0.25 * I))


if __name__ == "__main__":
    print("checks: F_half(30)=%.5f f_half(30)=%.5f F_lin(30)=%.5f f_lin(30)=%.5f"
          % (F_half(30.0), f_half(30.0), F_lin(30.0), f_lin(30.0)))
    print("f_lin(3) vs 2e^g ln2/3:", float(f_lin(3.0)), 2 * EG * np.log(2) / 3)
    print("f_half(2) vs closed form:", float(f_half(2.0)), _A_HALF * np.arccosh(np.sqrt(2)) / np.sqrt(2))
