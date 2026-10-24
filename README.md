# Sums of two squares with few prime factors via theta series of quaternary forms

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23618710.svg)](https://doi.org/10.5281/zenodo.23618710)

Code, data and manuscript for

> R. Chen, *Sums of two squares with few prime factors via theta series of quaternary forms*, preprint, 2026.

## What the paper does

Let N ≡ 2 (mod 4) be large. The paper shows that there are ≫ N/(log N)² integers n ≤ N such that n and N − n
are both sums of two squares, Ω(n(N − n)) ≤ 9 and max(Ω(n), Ω(N − n)) ≤ 6. The previous bounds (obtained by
sifting the primes ≡ 3 mod 4 and a switching argument, [bilateral-two-squares-sieve](https://github.com/Ruqing1963/bilateral-two-squares-sieve))
were 11 and 8.

The new input is a level of distribution 2/5 for the weighted sequence r(n) r(N − n), where r counts
representations as a sum of two squares. For coprime m, m′ the sum T(m, m′; N) = Σ_{ma + m′b = N} r(a) r(b) is the
N-th coefficient of the weight-2 theta series θ²(mz)θ²(m′z) of the quaternary form m(x² + y²) + m′(z² + w²).
Siegel's theorem gives the main term as a product of local densities; Waibel's uniform bound for the Petersson
norm of the cuspidal part of a theta series (Math. Proc. Camb. Phil. Soc. 173 (2022)) gives ‖G‖ ≪ (mm′)^{1/4+ε}
for the cusp part; and Deligne's bound, applied through the Hecke operator T_{n₁} rather than through the
Petersson formula directly, turns this into

    T(m, m′; N) = M(m, m′; N) + O(d(N) N^{1/2} (mm′)^{1/4+ε})      uniformly in m, m′.

The sieve is then a vector sieve of two linear components with Richert weights; its numerical inequality is
certified with interval arithmetic. The lower bound is c₁ 𝔖(N) N/(log N)² with c₁ = 0.19 and
𝔖(N) = ∏_{p | N, p > 2} (p − 1)/(p − 2), the Goldbach singular series (the events p | n and p | N − n are mutually
exclusive, so the sifted density is ∏(1 − 2h(p)); the small primes are handled by the fundamental lemma).
Numerically the error does not grow with mm′ at all, which suggests level 1/2 and the bounds 8 and 5
(Conjecture 5.1 of the paper).

## Repository layout

```
paper/      two_squares_theta.tex, two_squares_theta.pdf (manuscript) and the four figure PDFs
code/       sup_error_scan.py         T(m,m';N) for all N <= X by FFT, exact Siegel main term, sup/rms statistics
            certify_modular.py        interval-arithmetic certification of the sieve inequality (Table 2, Remark 4.5)
            modular_sieve_scan.py     floating-point parameter scan for the sieve (uses sieve_core.py)
            sieve_core.py             linear/semi-linear sieve functions (from bilateral-two-squares-sieve)
            singular_constant.py      the Mertens-type constants c_h, c_2h = kappa_0 c_h^2 and c_1 of Theorem 1.1
            plot_sup_error.py         Figure (sup/rms against mm') from the scan output
            plot_error_hist.py        Figure (distribution of the normalised error, Section 5.3)
            plot_sieve_landscape.py   Figure (L - lam R in the (lam, beta1) plane, Section 4)
            plot_sieve_functions.py   Figure (enclosures of F_1, f_1, Appendix A.1)
figures/    sup_error, error_hist, sieve_landscape, sieve_functions  (.pdf and .png)
results/    sup_error_scan_x21.txt, sup_error_scan_x23.txt      (Table 3)
            certify_modular_output.txt, certify_modular_output_half.txt   (Table 2, Remark 4.5)
            modular_sieve_scan_output.txt, singular_constant_output.txt
```

## Reproducing the numbers in the paper

Requirements: Python ≥ 3.10 with `numpy`, `scipy`, `mpmath`, `matplotlib`. The output files in `results/` and the
figures were produced with Python 3.12.4, numpy 1.26.4, scipy 1.13.1, mpmath 1.3.0, matplotlib 3.8.4 (Windows 11);
these versions are pinned in `requirements.txt` (`pip install -r requirements.txt`). Each script runs in a few
minutes on a laptop.

```
cd code
python certify_modular.py              # Table 2: level 2/5, three parameter sets (certified)
python certify_modular.py 1/2          # Remark 4.5: level 1/2, conditional
python singular_constant.py            # c_h in [1.501360, 1.501362], c_2h in [1.66485, 1.66487], c_1 >= 0.193
python sup_error_scan.py               # Table 3, X = 2^21, all pairs (~1 min)
python sup_error_scan.py 8388608 1,5 9,1 3,5 1,101 13,29 61,67 1,1009 101,103 1009,1013   # X = 2^23 (~3 min)
python plot_sup_error.py ../results/sup_error_scan_x21.txt    # figure: sup and rms against mm'
python plot_error_hist.py              # figure: distribution of (T-M)/sqrt(N) (~1 min)
python plot_sieve_landscape.py         # figure: L - lam R in the (lam, beta1) plane (floating point)
python plot_sieve_functions.py         # figure: enclosures of F_1, f_1
python modular_sieve_scan.py 0.4       # exploration of the sieve parameters (not certified)
python modular_sieve_scan.py 0.333333  # same at level 1/3 (what the Petersson formula alone would give): 10 and 6
```

### What is rigorous and what is not

| Script | Arithmetic | Enters the proof? |
|---|---|---|
| `certify_modular.py` | closed forms in `mpmath` interval arithmetic (80 bits); recurrences and sums in IEEE doubles with one outward `nextafter` step per operation; 1/λ + 2/β₁ and 1/λ + 1/β₁ exact rationals | yes: Table 2 and Remark 4.5 are read off its output lines `L >=`, `R <=`, `certified: net = L - lam R >=` |
| `singular_constant.py` | floating point with explicit tail and rounding margins | only the value of the constant c₁ in Theorem 1.1 |
| `sup_error_scan.py`, `plot_*.py`, `modular_sieve_scan.py` | floating point, no error control | no: evidence (Section 5) and illustration (figures) |

`certify_modular.py` evaluates the linear-sieve functions F₁, f₁ from their closed forms in `mpmath` interval
arithmetic and continues them by the delay-differential equations with outward rounding; L is bounded below and
R above (Riemann sum with monotonicity), and 1/λ + 2/β₁, 1/λ + 1/β₁ are computed in exact rational arithmetic.

`sup_error_scan.py` computes T(m, m′; N) for all N ≤ X at once by an FFT convolution and the Siegel main term from
the exact local densities (p ∤ 2mm′: (1 − p⁻²)σ₋₁(p^{v_p(N)}); p | mm′: 1 − χ₋₄(p)/p; p = 2: exact count modulo 2⁹).
The identities T = M for the levels 4 and 12 (no cusp forms) are built-in consistency checks.

## Citation

See `CITATION.cff`. The manuscript is archived at Zenodo: DOI [10.5281/zenodo.23618710](https://doi.org/10.5281/zenodo.23618710).

## License

Code: MIT (see `LICENSE`). Manuscript (`paper/`): CC BY 4.0.
