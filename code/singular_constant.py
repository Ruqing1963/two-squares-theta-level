"""
The Mertens-type constant of the sieve of Paper III:

    c_h := lim_{z -> oo} (log z) prod_{2 < p < z} (1 - h(p)),
    h(p) = (2p-1)/(p(p+1))  for p = 1 mod 4,   h(q) = 1/(q(q-1))  for q = 3 mod 4.

Method.  prod_{2<p<z}(1-h(p)) = P_1(z)^2 * prod_{2<p<z} (1-h(p)) (1-1/p)^{-2[p = 1 (4)]},
where P_1(z) = prod_{p<z, p=1(4)} (1 - 1/p) ~ (2 e^{-gamma}/kappa)^{1/2} (log z)^{-1/2},
kappa = (pi/4) prod_{q = 3 (4)} (1 - q^{-2})    (from L(1, chi_{-4}) = pi/4 and Mertens' theorem).
Hence c_h = (2 e^{-gamma}/kappa) * prod_{2<p} (1-h(p)) (1-1/p)^{-2[p=1(4)]}.
Both products converge absolutely (terms 1 + O(p^{-2})); we truncate at P and bound the tails by
exp(+-sum_{p > P} c/p^2) with sum_{p>P} 1/p^2 < 1/(P log P) ... crude: < 2/P (P >= 10^6), and print
rigorous-style enclosures (floating point with generous tails).
"""
import numpy as np
from math import exp, log, pi
from mpmath import mp, euler

mp.dps = 20
P = 10 ** 7


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


pr = primes_upto(P)
odd = pr[pr > 2]
split = odd[odd % 4 == 1].astype(float)
inert = odd[odd % 4 == 3].astype(float)

# kappa
log_kappa_prod = np.sum(np.log1p(-1.0 / inert ** 2))
tail = 2.0 / P                      # |sum_{q > P} q^{-2}| < 2/P, generous
kappa_lo, kappa_hi = (pi / 4) * exp(log_kappa_prod - tail), (pi / 4) * exp(log_kappa_prod)

# correction product: split factors (1-h)/(1-1/p)^2 = p (p^2 - p + 1)/((p+1)(p-1)^2) = 1 + (2p-1)/((p+1)(p-1)^2)
# inert factors 1 - 1/(q(q-1))
log_corr = np.sum(np.log1p((2 * split - 1) / ((split + 1) * (split - 1) ** 2))) + np.sum(np.log1p(-1.0 / (inert * (inert - 1))))
tail_corr = 2 * 2.0 / P             # split tail < sum 2/p^2 < 2/P ; inert tail < sum 1/q^2 < 2/P
corr_lo, corr_hi = exp(log_corr - tail_corr), exp(log_corr + tail_corr)

two_eg = 2 * exp(-float(euler))
c_lo, c_hi = two_eg / kappa_hi * corr_lo, two_eg / kappa_lo * corr_hi
print(f"primes up to {P}: kappa in [{kappa_lo:.8f}, {kappa_hi:.8f}], correction product in [{corr_lo:.8f}, {corr_hi:.8f}]")
print(f"c_h in [{c_lo:.6f}, {c_hi:.6f}]")

# direct check: (log z) prod_{2<p<z}(1-h(p)) at z = P (converges slowly, O(1/log z) error)
h = np.where(odd % 4 == 1, (2 * odd - 1) / (odd * (odd + 1.0)), 1.0 / (odd * (odd - 1.0)))
print(f"direct: (log P) * prod_(2<p<P) (1-h(p)) = {log(P) * exp(np.sum(np.log1p(-h))):.6f}  (slow convergence, for orientation only)")

# constant c_1 of Theorem 1.1:  c_1 = net * c_h^2 / (1024 alpha^2)  with net >= 0.153967, alpha = 9/250
net, alpha = 0.153967, 9 / 250
print(f"c_1 = net c_h^2/(1024 alpha^2) >= {net * c_lo ** 2 / (1024 * alpha ** 2):.4f}")
