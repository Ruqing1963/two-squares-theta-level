"""
First-principles measurement of the level dependence of the error in

    T(m, m') = #{(x, y, z, w) in Z^4 : m(x^2 + y^2) + m'(z^2 + w^2) = N} = sum_{m a + m' b = N} r(a) r(b),

for coprime odd m, m' and all N <= X at once (FFT convolution), against Siegel's main term

    M(N) = pi^2 N / (m m') * prod_p delta_p(N),
    delta_p(N) = 1 - chi_{-4}(p)/p                     (p | m m', p odd, p nmid N;  Hensel from z^2 + w^2 = N/m' mod p)
    delta_p(N) = (1 - p^{-2}) sigma_{-1}(p^{v_p(N)})    (p nmid 2 m m')
    delta_2(N) = exact count mod 2^nu, nu = v_2(N) + 6 (checked against nu + 1).

T - M is the N-th coefficient of the cusp part of theta^2(mz) theta^2(m'z), a holomorphic form of weight 2 on
Gamma_0(4 m m').  Deligne gives |T - M| <= C(m, m') d(N) N^{1/2}; the question is how C(m, m') grows with the
level.  We measure, over all admissible N in (X/2, X]:

    sup |T - M| / sqrt(N),  rms,  99.9 % quantile,  and the N where the sup is attained.

Usage:  python sup_error_scan.py [X] [m,m' ...]      (X a power of two, default 2^21)
"""
import sys
import time
import numpy as np
from math import gcd, pi
from scipy.signal import fftconvolve


# ---------------------------------------------------------------- arithmetic tables
def r_table(X):
    """r(n) for 0 <= n <= X, r(0) = 1."""
    chi = np.zeros(X + 1)
    idx = np.arange(1, X + 1)
    chi[1:] = np.where(idx % 2 == 0, 0.0, np.where(idx % 4 == 1, 1.0, -1.0))
    s = np.zeros(X + 1)
    for d in range(1, X + 1, 2):
        s[d::d] += chi[d]
    r = 4 * s
    r[0] = 1.0
    return r


def sigma_over_n_odd(X):
    """f(N) = sigma(N_odd)/N_odd for 0 <= N <= X (N_odd = odd part of N); f(0) := 0."""
    sig = np.zeros(X + 1)
    for d in range(1, X + 1):
        sig[d::d] += d
    v2 = np.zeros(X + 1, dtype=np.int64)
    n = np.arange(X + 1)
    t = n.copy()
    while True:
        even = (t % 2 == 0) & (t > 0)
        if not even.any():
            break
        t[even] //= 2
        v2[even] += 1
    f = np.zeros(X + 1)
    f[1:] = sig[t[1:]] / t[1:]
    return f, v2


def chi4(p):
    return 0 if p % 2 == 0 else (1 if p % 4 == 1 else -1)


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def counts_mod(m, mp, q):
    """#{(x,y,z,w) mod q : m(x^2+y^2) + m'(z^2+w^2) = t mod q} for all t mod q (exact for q^4 < 2^53)."""
    sq = np.bincount((np.arange(q, dtype=np.int64) ** 2) % q, minlength=q).astype(float)
    two = np.round(np.fft.irfft(np.fft.rfft(sq) ** 2, q))
    A = np.bincount((np.arange(q) * m) % q, weights=two, minlength=q)
    B = np.bincount((np.arange(q) * mp) % q, weights=two, minlength=q)
    return np.round(np.fft.irfft(np.fft.rfft(A) * np.fft.rfft(B), q))


def delta2_table(m, mp):
    """delta_2(N) as a function of N mod 2^9 for v_2(N) <= 3, with a stability check nu vs nu + 1."""
    tabs = {nu: counts_mod(m, mp, 2 ** nu) for nu in range(6, 11)}
    dens = np.zeros(2 ** 9)
    for N in range(2 ** 9):
        v = 0
        t = N
        while t and t % 2 == 0:
            t //= 2
            v += 1
        if N == 0 or v > 3:
            continue
        nu = v + 6
        d1 = tabs[nu][N % 2 ** nu] / 8.0 ** nu
        d2 = tabs[nu + 1][N % 2 ** (nu + 1)] / 8.0 ** (nu + 1)
        assert abs(d1 - d2) < 1e-12 * max(1.0, d1), (m, mp, N, d1, d2)
        dens[N] = d1
    return dens


def check_odd_bad_prime_density(m, mp):
    """Compare delta_p = 1 - chi4(p)/p (p | mm') with exact counts mod p^2 for the small primes dividing mm'."""
    for p in factor(m * mp):
        if p > 13:
            continue
        q = p * p
        cnt = counts_mod(m, mp, q)
        for N in range(1, q):
            if N % p:
                assert abs(cnt[N] / q ** 3 - (1 - chi4(p) / p)) < 1e-9, (p, N, cnt[N] / q ** 3)


# ---------------------------------------------------------------- T and M for all N <= X
def T_all(r, m, mp, X):
    A = np.zeros(X + 1)
    A[::m] = r[: X // m + 1]
    B = np.zeros(X + 1)
    B[::mp] = r[: X // mp + 1]
    T = fftconvolve(A, B)[: X + 1]
    return np.round(T)


def M_all(m, mp, X, f_sigma, v2, dens2):
    N = np.arange(X + 1, dtype=np.float64)
    bad = factor(m * mp)
    good = (6 / pi ** 2) / np.prod([1 - p ** -2.0 for p in list(bad) + [2]])
    local_bad = np.prod([1 - chi4(p) / p for p in bad])
    M = pi ** 2 * N / (m * mp) * good * local_bad * f_sigma * dens2[np.arange(X + 1) % 512]
    return M


def admissible(X, m, mp, v2):
    n = np.arange(X + 1)
    ok = (n > X // 2) & (v2 <= 3)
    for p in factor(m * mp):
        ok &= n % p != 0
    return ok


def report(label, N, e, T, M):
    if len(e) == 0:
        print(f"  {label}: no N")
        return
    i = np.argmax(np.abs(e))
    Nmax = int(N[i])
    fa = factor(Nmax)
    dN = int(np.prod([v + 1 for v in fa.values()]))
    q999 = np.quantile(np.abs(e), 0.999)
    print(f"  {label:<12} n={len(e):>7}  rms={np.sqrt(np.mean(e ** 2)):7.3f}  q99.9={q999:7.3f}  "
          f"sup={abs(e[i]):7.3f}  at N={Nmax} = {'*'.join(f'{p}^{v}' if v > 1 else str(p) for p, v in sorted(fa.items()))}"
          f"  d(N)={dN}  T={int(T[i])}  M={M[i]:.1f}")


if __name__ == "__main__":
    args = sys.argv[1:]
    X = int(args[0]) if args and ',' not in args[0] else 2 ** 21
    pairs = [tuple(int(t) for t in a.split(',')) for a in args if ',' in a]
    if not pairs:
        pairs = [(1, 1), (1, 3), (1, 5), (1, 7), (1, 13), (3, 5), (3, 7), (5, 7), (9, 1), (25, 1), (9, 25),
                 (1, 101), (13, 29), (15, 77), (61, 67), (1, 1009), (105, 143), (101, 103), (1, 15015),
                 (1009, 1013), (3, 100003)]
    t0 = time.time()
    r = r_table(X)
    f_sigma, v2 = sigma_over_n_odd(X)
    print(f"X = {X} = 2^{X.bit_length() - 1};  tables in {time.time() - t0:.1f}s.  "
          f"Statistics of (T - M)/sqrt(N) over admissible N in (X/2, X] (gcd(N, mm') = 1, v_2(N) <= 3).")
    for (m, mp) in pairs:
        assert gcd(m, mp) == 1 and m % 2 == 1 and mp % 2 == 1
        t1 = time.time()
        check_odd_bad_prime_density(m, mp)
        dens2 = delta2_table(m, mp)
        T = T_all(r, m, mp, X)
        M = M_all(m, mp, X, f_sigma, v2, dens2)
        ok = admissible(X, m, mp, v2)
        n = np.arange(X + 1)
        print(f"(m, m') = ({m}, {mp}), mm' = {m * mp}, level 4mm' = {4 * m * mp}  [{time.time() - t1:.1f}s]")
        for label, mask in (("all", ok), ("N odd", ok & (n % 2 == 1)), ("N = 2 mod 4", ok & (n % 4 == 2))):
            NN = n[mask]
            e = (T[mask] - M[mask]) / np.sqrt(NN)
            report(label, NN, e, T[mask], M[mask])
        sys.stdout.flush()
