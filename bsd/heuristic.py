"""The original 1960s Birch–Swinnerton-Dyer heuristic.

Birch and Swinnerton-Dyer computed, on EDSAC-2, the partial products
    f_E(X) = ∏_{p ≤ X, p ∤ N} N_p / p,      N_p = #E(F_p) = p + 1 - a_p,
and observed f_E(X) ≈ C (log X)^r with r the rank.  This heuristic was the
origin of the conjecture.  It is *not* a theorem and is numerically rather
poor (the convergence is slow and oscillatory); Goldfeld (1982) showed that
if f_E(X) ~ C (log X)^r then r equals the analytic rank, and the constant C
involves L^{(r)}(E,1)/r! times a factor √2 explained by Conrad (2005).
We reproduce the experiment for illustration only.
"""
from __future__ import annotations

import math

from .pari_session import get_pari


def partial_products(ainvs, X: int, n_samples: int = 40) -> list[tuple[int, float]]:
    """Return [(X_i, log f_E(X_i))] for X_i spaced geometrically up to X."""
    pari = get_pari()
    E = pari.ellinit(list(ainvs))
    N = int(pari.ellglobalred(E)[0])
    an = pari.ellan(E, X)                      # a_n for n ≤ X (a_p = p + 1 - N_p at good p)
    primes = pari.primes(int(pari.primepi(X)))
    checkpoints = sorted({int(round(math.exp(math.log(X) * (k + 1) / n_samples))) for k in range(n_samples)})
    out = []
    logf = 0.0
    ci = 0
    for p in primes:
        p = int(p)
        while ci < len(checkpoints) and checkpoints[ci] < p:
            out.append((checkpoints[ci], logf))
            ci += 1
        if N % p == 0:
            continue
        Np = p + 1 - int(an[p - 1])
        logf += math.log(Np / p)
    while ci < len(checkpoints):
        out.append((checkpoints[ci], logf))
        ci += 1
    return [(x, v) for (x, v) in out if x >= 30]


def fitted_exponent(points: list[tuple[int, float]]) -> float:
    """Least-squares slope of log f(X) against log log X (should be ≈ rank)."""
    xs = [math.log(math.log(x)) for x, _ in points]
    ys = [v for _, v in points]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / sxx
