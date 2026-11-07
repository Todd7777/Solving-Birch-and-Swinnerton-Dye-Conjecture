"""Rigorous ball enclosures of the Taylor coefficients of L(E, s) at s = 1.

For an elliptic curve E/Q of conductor N this module computes Arb balls
(python-flint) that provably contain

    c_k = L^{(k)}(E, 1) / k!            and      lambda_k = Lambda^{(k)}(1) / k!,

for k = 0..kmax, and uses them to certify orders of vanishing.  It does not
assume that lower derivatives vanish.  Nothing in this module proves anything
about BSD itself; it gives certified statements about the analytic side only.

Mathematics (re-derived in research/rigorous-l/report.md, section 1)
--------------------------------------------------------------------
By the modularity theorem, L(E, s) = L(f, s) for a newform f = sum a_n q^n of
weight 2 on Gamma_0(N).  Put

    Lambda(s) = N^{s/2} (2 pi)^{-s} Gamma(s) L(E, s),
    g(u)      = f(i u / sqrt N) = sum_n a_n exp(-x_n u),   x_n = 2 pi n / sqrt N.

The Fricke involution gives g(1/u) = w u^2 g(u), where w = +-1 is the global
root number (w equals minus the Atkin-Lehner eigenvalue of f).  Splitting the
Mellin integral at u = 1 then gives the entire function

    Lambda(1 + t) = int_1^oo g(u) (u^t + w u^{-t}) du,

so lambda_k = (1 + w (-1)^k) * sum_n a_n K_k(x_n), where

    K_k(x) = J_k(x) / k!,     J_k(x) = int_1^oo e^{-x u} (log u)^k du.

The K_k(x) are the Taylor coefficients in t of int_1^oo e^{-xu} u^t du =
x^{-1-t} Gamma(1+t, x).  From Gamma(1+t, x) = Gamma(1+t) - gamma(1+t, x) and
Kummer's transformation 1F1(1+t; 2+t; -x) = e^{-x} 1F1(1; 2+t; x):

    x^{-1-t} Gamma(1+t, x) = x^{-1-t} Gamma(1+t) - e^{-x} sum_{m>=0} x^m / (1+t)_{m+1}.

Write q_{k,m} = |[t^k] 1/(1+t)_{m+1}| = h_k(1, 1/2, ..., 1/(m+1)) / (m+1)!
(h_k is the complete homogeneous symmetric polynomial).  The coefficient of
t^k in 1/(1+t)_{m+1} is (-1)^k q_{k,m}, so

    K_k(x) = (1/x) sum_{i=0}^k gamma_i (-log x)^{k-i} / (k-i)!  -  (-1)^k e^{-x} Q_k(x),
    Q_k(x) = sum_{m>=0} q_{k,m} x^m,                  Gamma(1+t) = sum_i gamma_i t^i.

Q_k has positive coefficients, so evaluating it at x > 0 does not cancel.  The
form with e^{-x} pulled out keeps absolute accuracy about 2^-prec even for
x ~ 100.  The alternating series for gamma(1+t, x) would lose about 1.44 x bits.

Finally L(E, 1+t) = C(t) Lambda(1+t), where

    C(t) = (2 pi / sqrt N)^{1+t} / Gamma(1+t)
         = (2 pi / sqrt N) exp( t (log(2 pi / sqrt N) + gamma) - sum_{m>=2} (-1)^m zeta(m) t^m / m ),

and c_k = sum_{j<=k} lambda_j C_{k-j}.  The exponent is truncated at t^kmax.
That truncation is exact for the coefficients up to t^kmax.

Proven error bounds (proofs in the report, section 2)
-----------------------------------------------------
(B1) |a_n| <= d(n) sqrt(n) <= 2n         (Hasse bound, multiplicativity, d(n) <= 2 sqrt n).
(B2) 0 <= K_k(x) <= e^{-x} / x^{k+1}     (log u <= u - 1 on [1, oo)).
(B3) n-tail:  |sum_{n>M} a_n K_k(x_n)| <= 2 e^{-a(M+1)} / (a^{k+1} (M+1)^k (1 - e^{-a})),  a = 2 pi / sqrt N.
(B4) m-tail of Q_k:  q_{k,m} <= H_{m+1}^k / (m+1)!, hence for rho = (1 + 1/(D+3))^k x / (D+3) < 1,
     sum_{m>D} q_{k,m} x^m <= H_{D+2}^k x^{D+1} / ((D+2)! (1 - rho)).
Every bound is evaluated in ball arithmetic and added to the result as a radius.
All constants (pi, Euler's gamma, zeta(m), logs) are Arb balls.

Inputs and the assumptions behind them
--------------------------------------
* a_n: PARI ellan, exact integers.  N: PARI ellglobalred (Tate's algorithm), exact.
  w: PARI ellrootno, an exact algorithm (product of local root numbers).
* Modularity of E/Q (Wiles, Taylor-Wiles, Breuil-Conrad-Diamond-Taylor) gives
  the newform f, the analytic continuation and the functional equation with sign w.
* Correctness of the Arb/FLINT ball arithmetic (python-flint) and of the PARI
  routines above.  This is a software-correctness assumption, as in any
  computer-assisted proof.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Iterable, Sequence

from flint import arb, arb_poly, arb_series, ctx

DEFAULT_PREC = 128          # target precision in bits (radii ~ 2^-128 ~ 3e-39)
DEGREE_STEP = 8             # granularity of the truncation degrees D used for Q_k

ASSUMPTIONS = [
    "Modularity theorem (Wiles; Taylor-Wiles; Breuil-Conrad-Diamond-Taylor): L(E,s) = L(f,s) for a weight-2 "
    "newform f of level N, so Lambda(s) is entire and Lambda(s) = w Lambda(2-s) with w the global root number; "
    "equivalently g(1/u) = w u^2 g(u) for g(u) = f(iu/sqrt N).",
    "a_n from PARI ellan (exact integers); conductor N from PARI ellglobalred (Tate's algorithm, exact); "
    "root number w from PARI ellrootno (exact algorithm: product of local root numbers).",
    "Hasse's bound |a_p| <= 2 sqrt p (theorem), giving |a_n| <= d(n) sqrt n <= 2n for all n >= 1.",
    "Correctness of the Arb ball arithmetic (FLINT 3 via python-flint) and of the PARI routines used.",
]


# ----------------------------------------------------------------------------
# JSON helpers
# ----------------------------------------------------------------------------

def _dec_upper(x: arb, digits: int = 6) -> str:
    """A decimal string that is a rigorous upper bound for |x|."""
    if x.is_zero():
        return "0"
    m, r, e = x.mid_rad_10exp(digits)
    u = abs(int(m)) + int(r)
    s = str(u)
    exp10 = int(e) + len(s) - 1
    mant = s[0] + ("." + s[1:digits + 2] if len(s) > 1 else "")
    # truncating digits of u lowers the value, so bump the last kept digit up
    if len(s) > digits + 2:
        kept = int(s[: digits + 2]) + 1
        ks = str(kept)
        if len(ks) > digits + 2:          # carry, e.g. 999 -> 1000
            exp10 += 1
            ks = ks[:-1]
        mant = ks[0] + "." + ks[1:]
    return f"{mant}e{exp10}"


def ball_to_json(x: arb, digits: int = 40) -> dict:
    """JSON-friendly form of a ball.

    'ball' is Arb's decimal string, which is guaranteed to enclose x (arb_get_str).
    'mid' is a decimal midpoint for display.  'rad' is a rigorous decimal upper bound
    for the radius of the ball as stored.  'abs_upper' is a rigorous upper bound for |x|.
    'exact_zero' is True only for an exact zero, such as a coefficient that vanishes
    by parity.
    """
    exact_zero = bool(x.is_zero())
    out = {
        "ball": x.str(digits, radius=True),
        "mid": x.mid().str(digits, radius=False) if not exact_zero else "0",
        "rad": _dec_upper(arb(x.rad())) if not exact_zero else "0",
        "rad_float": float(x.rad()) if not exact_zero else 0.0,
        "abs_upper": _dec_upper(x),
        "contains_zero": bool(x.contains(0)),
        "exact_zero": exact_zero,
    }
    return out


# ----------------------------------------------------------------------------
# The special functions K_k(x) = J_k(x)/k!
# ----------------------------------------------------------------------------

class _QTable:
    """Coefficients q_{k,m} = |[t^k] 1/(1+t)_{m+1}| as balls, plus data for the tail bound (B4).

    Recurrence: q_{k,0} = 1 and q_{k,m+1} = sum_{i=0}^k q_{k-i,m} / (m+2)^{i+1}, from
    1/(1+t)_{m+2} = (1/(1+t)_{m+1}) * (1/(m+2+t)).
    """

    def __init__(self, K: int, prec: int):
        self.K = K
        self.prec = prec
        self.q: list[list[arb]] = [[] for _ in range(K)]
        self.D = -1
        self.polys: dict[tuple[int, int], arb_poly] = {}
        self.fac: dict[int, arb] = {}        # D -> (D+2)!
        self.harm: dict[int, arb] = {}       # D -> H_{D+2}
        self._cur: list[arb] | None = None

    def extend(self, D: int) -> None:
        if D <= self.D:
            return
        with ctx.workprec(self.prec + 32):
            if self._cur is None:
                self._cur = [arb(1)] * self.K
                for k in range(self.K):
                    self.q[k].append(arb(1))
                self.D = 0
            cur = self._cur
            for m in range(self.D, D):
                c = arb(m + 2)
                inv = [c ** (-(i + 1)) for i in range(self.K)]
                cur = [sum((cur[k - i] * inv[i] for i in range(k + 1)), arb(0)) for k in range(self.K)]
                for k in range(self.K):
                    self.q[k].append(cur[k])
            self._cur = cur
            self.D = D

    def poly(self, k: int, D: int) -> arb_poly:
        key = (k, D)
        p = self.polys.get(key)
        if p is None:
            self.extend(D)
            with ctx.workprec(self.prec):
                p = arb_poly(self.q[k][: D + 1])
            self.polys[key] = p
        return p

    def tail_consts(self, D: int) -> tuple[arb, arb]:
        """((D+2)!, H_{D+2}) as balls."""
        if D not in self.fac:
            with ctx.workprec(self.prec):
                self.fac[D] = arb.fac_ui(D + 2)
                h = arb(0)
                for i in range(1, D + 3):
                    h += arb(1) / i
                self.harm[D] = h
        return self.fac[D], self.harm[D]


_QTABLES: dict[tuple[int, int], _QTable] = {}


def _qtable(K: int, prec: int) -> _QTable:
    key = (K, prec)
    t = _QTABLES.get(key)
    if t is None:
        t = _QTable(K, prec)
        _QTABLES[key] = t
    return t


def _choose_degree(x: float, k: int, prec_w: int) -> int:
    """Pick D (a multiple of DEGREE_STEP) whose tail bound (B4), times e^{-x}, is below about 2^-prec_w.

    Only a heuristic: the rigorous tail bound for the chosen D is added to the result in any case.
    """
    D = DEGREE_STEP
    target = -prec_w * math.log(2) + x
    lx = math.log(x)
    while True:
        rho = (1 + 1 / (D + 3)) ** k * x / (D + 3)
        if rho < 0.75:
            H = math.log(D + 2) + 0.5772156649 + 1 / (2 * (D + 2))
            est = k * math.log(H) + (D + 1) * lx - math.lgamma(D + 3) - math.log(1 - rho)
            if est < target:
                return D
        D += DEGREE_STEP


def _tail_bound_Q(x: arb, k: int, D: int, table: _QTable) -> arb:
    """Rigorous upper bound (B4) for sum_{m>D} q_{k,m} x^m, as a ball (take its upper end)."""
    fac, H = table.tail_consts(D)
    rho = (1 + arb(1) / (D + 3)) ** k * x / (D + 3)
    if not (rho < 1):
        raise ValueError("tail bound (B4) needs rho < 1; increase D")
    return H ** k * x ** (D + 1) / (fac * (1 - rho))


def K_values(x: arb, ks: Sequence[int], prec_w: int, gamma_coeffs: Sequence[arb] | None = None) -> dict[int, arb]:
    """Balls for K_k(x) = J_k(x)/k! = (1/k!) int_1^oo e^{-xu} (log u)^k du, for each k in ks (x > 0).

    Rigorous: Arb ball arithmetic plus the truncation bound (B4) added as a radius.
    """
    K = max(ks) + 1
    table = _qtable(K, prec_w)
    with ctx.workprec(prec_w):
        x = arb(x)
        if not (x > 0):
            raise ValueError("K_values needs x > 0")
        if gamma_coeffs is None:
            gamma_coeffs = arb_series([1, 1], prec=K).gamma().coeffs()
        lx = x.log()
        invx = 1 / x
        ex = (-x).exp()
        xf = float(x.mid())
        mlx_pows = [arb(1)]
        for i in range(1, K):
            mlx_pows.append(mlx_pows[-1] * (-lx) / i)          # (-log x)^i / i!
        out = {}
        for k in ks:
            A = sum((gamma_coeffs[i] * mlx_pows[k - i] for i in range(k + 1)), arb(0)) * invx
            D = _choose_degree(xf, k, prec_w)
            Q = table.poly(k, D)(x)
            Q += arb(0, _tail_bound_Q(x, k, D, table))
            val = A - ex * Q if k % 2 == 0 else A + ex * Q
            out[k] = val
        return out


# ----------------------------------------------------------------------------
# The L-function
# ----------------------------------------------------------------------------

def _ntail_bound(a: arb, M: int, k: int) -> arb:
    """(B3): rigorous bound for |sum_{n>M} a_n K_k(x_n)|, x_n = a n."""
    M1 = M + 1
    return 2 * (-a * M1).exp() / (a ** (k + 1) * arb(M1) ** k * (1 - (-a).exp()))


def choose_M(N: int, ks: Iterable[int], prec: int) -> int:
    """Smallest M with 2 * (B3) <= 2^-(prec+4) for every k in ks.

    This uses a float estimate, which is fine: the rigorous bound (B3) for the chosen M is
    added to the result in any case, so M only affects the size of the radius.
    """
    a = 2 * math.pi / math.sqrt(N)
    goal = -(prec + 4) * math.log(2)

    def logbound(M, k):
        M1 = M + 1
        return math.log(4) - a * M1 - (k + 1) * math.log(a) - k * math.log(M1) - math.log(-math.expm1(-a))

    ks = list(ks) or [0]
    M = max(1, int(((prec + 4) * math.log(2)) / a))
    while any(logbound(M, k) > goal for k in ks):
        M = int(M * 1.02) + 1
    while M > 1 and all(logbound(M - 1, k) <= goal for k in ks):
        M -= 1
    return M


def C_series(N: int, K: int, prec_w: int) -> list[arb]:
    """Taylor coefficients C_0..C_{K-1} of C(t) = (2 pi/sqrt N)^{1+t} / Gamma(1+t), via the zeta-exponent form."""
    with ctx.workprec(prec_w):
        L0 = (2 * arb.pi() / arb(N).sqrt()).log()
        coeffs = [arb(0)] * K
        if K > 1:
            coeffs[1] = L0 + arb.const_euler()
        for m in range(2, K):
            coeffs[m] = -((-1) ** m) * arb(m).zeta() / m
        expo = arb_series(coeffs, prec=K).exp()
        pref = 2 * arb.pi() / arb(N).sqrt()
        return [pref * c for c in expo.coeffs()] + [arb(0)] * (K - len(expo.coeffs()))


def C_series_rgamma(N: int, K: int, prec_w: int) -> list[arb]:
    """The same coefficients via Arb's rgamma power series (independent cross-check)."""
    with ctx.workprec(prec_w):
        L0 = (2 * arb.pi() / arb(N).sqrt()).log()
        s = arb_series([L0, L0], prec=K).exp() * arb_series([1, 1], prec=K).rgamma()
        cs = s.coeffs()
        return cs + [arb(0)] * (K - len(cs))


@dataclass
class TaylorEnclosure:
    ainvs: list[int]
    conductor: int
    root_number: int
    kmax: int
    prec: int
    prec_work: int
    M: int
    lam: list[arb]                 # Lambda^{(k)}(1)/k!
    c: list[arb]                   # L^{(k)}(E,1)/k!
    C: list[arb]                   # Taylor coefficients of C(t)
    ntail: dict[int, arb] = field(default_factory=dict)     # bound (B3) used for each k, already doubled for lambda
    timings: dict = field(default_factory=dict)

    def derivative(self, k: int) -> arb:
        """Ball for L^{(k)}(E,1) itself."""
        return self.c[k] * math.factorial(k)

    def to_json(self, digits: int = 40) -> dict:
        return {
            "ainvs": list(self.ainvs),
            "conductor": self.conductor,
            "root_number": self.root_number,
            "kmax": self.kmax,
            "prec_bits": self.prec,
            "prec_work_bits": self.prec_work,
            "M_terms": self.M,
            "c": [ball_to_json(x, digits) for x in self.c],
            "lambda": [ball_to_json(x, digits) for x in self.lam],
            "ntail_bound": {str(k): _dec_upper(v) for k, v in self.ntail.items()},
            "timings": self.timings,
        }


def _curve_data(ainvs, E=None):
    from .pari_session import get_pari
    pari = get_pari()
    if E is None:
        E = pari.ellinit([int(a) for a in ainvs])
    N = int(pari.ellglobalred(E)[0])
    w = int(pari.ellrootno(E))
    return pari, E, N, w


def taylor_enclosures(ainvs: Sequence[int], kmax: int = 4, prec: int = DEFAULT_PREC, *,
                      M: int | None = None, E=None, an: Sequence[int] | None = None) -> TaylorEnclosure:
    """Rigorous balls for c_k = L^{(k)}(E,1)/k! and lambda_k = Lambda^{(k)}(1)/k!, k = 0..kmax.

    The radii come out at about 2^-prec.  Larger prec gives smaller radii.
    """
    t0 = time.time()
    if kmax < 0:
        raise ValueError("kmax must be >= 0")
    K = kmax + 1
    pari, E, N, w = _curve_data(ainvs, E)
    ks = [k for k in range(K) if (-1) ** k == w]          # coefficients not killed by parity
    if M is None:
        M = choose_M(N, ks, prec) if ks else 0
    prec_w = prec + 24 + int(math.log2(M + 1)) + 4 * K
    t1 = time.time()
    if not ks:
        an = []
    elif an is None:
        an = [int(v) for v in pari.ellan(E, M)]
    else:
        an = [int(v) for v in an[:M]]
        if len(an) < M:
            raise ValueError("not enough a_n supplied")
    t2 = time.time()
    with ctx.workprec(prec_w):
        a = 2 * arb.pi() / arb(N).sqrt()
        gam = arb_series([1, 1], prec=K).gamma().coeffs()
        S = {k: arb(0) for k in ks}
        for n in range(1, M + 1):
            c = an[n - 1]
            if c == 0 or not ks:
                continue
            vals = K_values(a * n, ks, prec_w, gam)
            for k in ks:
                S[k] += c * vals[k]
        lam = []
        ntail = {}
        for k in range(K):
            if k in S:
                tb = 2 * _ntail_bound(a, M, k)              # lambda_k = 2 * S_k when (-1)^k = w
                ntail[k] = tb
                lam.append(2 * S[k] + arb(0, tb))
            else:
                lam.append(arb(0))                           # exact zero: parity
        C = C_series(N, K, prec_w)
        c_coeffs = []
        for k in range(K):
            c_coeffs.append(sum((lam[j] * C[k - j] for j in range(k + 1)), arb(0)))
    t3 = time.time()
    return TaylorEnclosure(ainvs=[int(v) for v in ainvs], conductor=N, root_number=w, kmax=kmax, prec=prec,
                           prec_work=prec_w, M=M, lam=lam, c=c_coeffs, C=C, ntail=ntail,
                           timings={"setup_s": round(t1 - t0, 4), "ellan_s": round(t2 - t1, 4),
                                    "series_s": round(t3 - t2, 4), "total_s": round(t3 - t0, 4)})


def lambda_by_integration(ainvs: Sequence[int], kmax: int = 4, prec: int = DEFAULT_PREC, *, E=None,
                          M: int | None = None) -> list[arb]:
    """Independent rigorous cross-check: lambda_k = Lambda^{(k)}(1)/k! via Arb's certified integrator.

    lambda_k = ((1 + w(-1)^k)/k!) int_1^oo g(u) (log u)^k du.  g is truncated at n <= M, so that
    g_M(u) = P(e^{-a u}) with P = sum_{n<=M} a_n X^n an integer polynomial.  The integral over
    [1, U] is computed by acb.integral (Petras' Gauss-Legendre scheme with rigorous error bounds).
    Two proven tails are added:
      u-tail: |int_U^oo g_M(u) (log u)^k du| <= sum_{n<=M} |a_n| Gamma(k+1, x_n U) / x_n^{k+1},
              using |g_M(u)| <= sum |a_n| e^{-x_n u} and log u <= u;
      n-tail: |sum_{n>M} a_n J_k(x_n)| <= k! * (B3).
    It shares none of the K_k machinery of taylor_enclosures, so overlapping balls check both.
    It is slower and slightly less tight, and is used for validation only.
    """
    from flint import acb, acb_poly
    pari, E, N, w = _curve_data(ainvs, E)
    K = kmax + 1
    ks = [k for k in range(K) if (-1) ** k == w]
    if M is None:
        M = choose_M(N, ks, prec) if ks else 0
    an = [int(v) for v in pari.ellan(E, M)] if ks else []
    prec_w = prec + 24 + int(math.log2(M + 1)) + 4 * K
    out = []
    with ctx.workprec(prec_w):
        a = 2 * arb.pi() / arb(N).sqrt()
        P = acb_poly([0] + an)
        aa = acb(a)
        U = int(math.ceil(((prec + 30) * math.log(2) + 10) / float(a.mid()))) + 1
        for k in range(K):
            if k not in ks:
                out.append(arb(0))
                continue

            def f(u, analytic, k=k):
                v = P((-aa * u).exp())
                return v * u.log(analytic=analytic) ** k if k else v

            I = acb.integral(f, 1, U).real
            utail = arb(0)
            for n, c in enumerate(an, start=1):
                if c:
                    xn = a * n
                    utail += abs(c) * (xn * U).gamma_upper(k + 1) / xn ** (k + 1)
            S = I + arb(0, utail) + arb(0, _ntail_bound(a, M, k) * math.factorial(k))
            out.append(2 * S / math.factorial(k))
    return out


def conditioned_coefficients(enc: TaylorEnclosure, zero_lambdas: Iterable[int]) -> list[arb]:
    """Recompute c_k after setting lambda_j = 0 exactly for j in zero_lambdas.

    Use this only for lambda_j proven to be zero, for instance by modular symbols.
    Then c_r for the first non-zero coefficient becomes exactly lambda_r C_0.
    """
    zl = set(zero_lambdas)
    lam = [arb(0) if j in zl else enc.lam[j] for j in range(enc.kmax + 1)]
    with ctx.workprec(enc.prec_work):
        return [sum((lam[j] * enc.C[k - j] for j in range(k + 1)), arb(0)) for k in range(enc.kmax + 1)]


# ----------------------------------------------------------------------------
# Certification of the order of vanishing
# ----------------------------------------------------------------------------

def modsym_L1(E, timeout: float | None = None) -> dict:
    """Exact status of L(E,1) from modular symbols (bsd.analytic.exact_L1_ratio).

    Returns {"computed": True, "ratio_abs": "p/q", "zero": bool, ...} or
    {"computed": False, "error": ...}.  The value |L(E,1)/Omega^+| comes from exact linear
    algebra over Q, so "zero" is an exact statement.  It does not depend on numerical precision.
    """
    from .analytic import exact_L1_ratio
    from .pari_session import time_limit
    t = time.time()
    try:
        with time_limit(timeout):
            r = exact_L1_ratio(E)
        return {"computed": True, "ratio_abs": str(r), "zero": r == 0, "seconds": round(time.time() - t, 3),
                "source": "modular symbols (PARI msfromell/mseval, exact over Q)"}
    except Exception as exc:  # pragma: no cover - depends on PARI resources
        return {"computed": False, "error": f"{type(exc).__name__}: {str(exc)[:160]}",
                "seconds": round(time.time() - t, 3)}


def certify_order_rigorous(ainvs: Sequence[int], kmax: int = 4, prec: int = DEFAULT_PREC, *,
                           use_modular_symbols: bool = True, modsym: dict | None = None,
                           modsym_timeout: float | None = None, E=None, label: str | None = None,
                           digits: int = 40, external_zeros: dict[int, str] | None = None,
                           external_lower_bound: tuple[int, str] | None = None) -> dict:
    """Certify ord_{s=1} L(E,s) as far as rigorous enclosures plus exact side information allow.

    Exact zeros come from parity (lambda_k = 0 when (-1)^k = -w) and from modular symbols
    (L(E,1) = 0 exactly; computed only when w = +1 and the ball for lambda_0 contains 0).
    Non-zero coefficients come from balls that exclude 0.
    Since C(0) = 2 pi/sqrt N != 0, c_0 = ... = c_{r-1} = 0 holds iff lambda_0 = ... = lambda_{r-1} = 0.
    So ord = r is certified when lambda_0..lambda_{r-1} are proven zero and c_r (or lambda_r)
    is proven non-zero.

    A ball around 0 never proves vanishing.  'undetermined' lists each k below the first
    certified non-zero coefficient for which lambda_k is neither proven zero nor proven non-zero.
    For each it gives rigorous bounds |c_k| <= eps, |lambda_k| <= eps' and |L^{(k)}(E,1)| <= k! eps.
    For root number -1 and analytic rank 3 this is k = 1: the output field 'L_prime_abs_upper'
    is the bound on |L'(E,1)| that a Gross-Zagier argument must turn into exact vanishing.

    `modsym` may carry a precomputed result of modsym_L1 (same format), e.g. one per isogeny class.

    Exact facts proven elsewhere can be supplied.  They are recorded as external and named
    as such in the conclusion:
      * external_zeros = {k: source}: Lambda^{(k)}(1) = 0 exactly, e.g. {1: "Gross-Zagier: Heegner point torsion"};
      * external_lower_bound = (r0, source): ord >= r0 proven, e.g. (3, "GZK contrapositive: rank >= 2 proven,
        root number -1"); this means lambda_j = 0 for every j < r0.
    A RuntimeError is raised if an enclosure contradicts such an input (its ball excludes 0).
    """
    t0 = time.time()
    pari, E, N, w = _curve_data(ainvs, E)
    enc = taylor_enclosures(ainvs, kmax=kmax, prec=prec, E=E)
    K = kmax + 1
    lam_zero: dict[int, str] = {}
    for k in range(K):
        if (-1) ** k == -w:
            lam_zero[k] = f"parity: Lambda(1+t) = w Lambda(1-t) with w = {w:+d}, so Lambda^({k})(1) = 0"
    l1_nonzero_exact = False
    if w == 1 and enc.lam[0].contains(0) and use_modular_symbols and modsym is None:
        modsym = modsym_L1(E, timeout=modsym_timeout)
    if w == 1 and modsym is not None:
        if modsym.get("computed"):
            if modsym["zero"]:
                lam_zero[0] = "modular symbols: L(E,1)/Omega^+ = 0 exactly (PARI msfromell/mseval)"
            else:
                l1_nonzero_exact = True
    external = {}
    if external_lower_bound is not None:
        r0, src = external_lower_bound
        for j in range(min(int(r0), K)):
            if j not in lam_zero:
                lam_zero[j] = f"external: ord >= {int(r0)} ({src})"
        external["lower_bound"] = {"ord_at_least": int(r0), "source": src}
    if external_zeros:
        for j, src in external_zeros.items():
            j = int(j)
            if j < K and j not in lam_zero:
                lam_zero[j] = f"external: {src}"
        external["zeros"] = {str(j): s for j, s in external_zeros.items()}
    if l1_nonzero_exact and 0 in lam_zero:
        raise RuntimeError("external input says L(E,1) = 0 but modular symbols say L(E,1) != 0")
    c_nonzero = [k for k in range(K) if not enc.c[k].contains(0)]
    lam_nonzero = [k for k in range(K) if not enc.lam[k].contains(0)]
    nonzero = sorted(set(c_nonzero) | set(lam_nonzero) | ({0} if l1_nonzero_exact else set()))
    # exact zeros of c_k: c_k = 0 iff lambda_0..lambda_k are all 0
    def _short(reason: str) -> str:
        head = reason.split(":")[0]
        return reason if head == "external" else head

    c_zero = {}
    for k in range(K):
        if all(j in lam_zero for j in range(k + 1)):
            c_zero[k] = "; ".join(sorted({_short(lam_zero[j]) for j in range(k + 1)}))
        else:
            break
    r_lower = len(c_zero)                                    # ord >= r_lower (proven)
    r_upper = nonzero[0] if nonzero else None                # ord <= r_upper (proven)
    parity = 0 if w == 1 else 1
    # consistency: an exactly-vanishing lambda_k must have a ball containing 0
    inconsistent = [k for k in lam_zero if not enc.lam[k].contains(0)]
    if inconsistent or (r_upper is not None and r_upper < r_lower):
        raise RuntimeError(f"inconsistent certificate data (lambda balls exclude 0 at {inconsistent}); "
                           "this indicates a bug or wrong input data")
    if r_upper is not None:
        candidates = [r for r in range(r_lower, r_upper + 1) if r % 2 == parity and r not in lam_zero]
    else:
        candidates = None
    certified = candidates is not None and len(candidates) == 1
    undetermined = {}
    if not certified:
        top = max(candidates) if candidates else K
        for k in range(r_lower, top):
            if k in lam_zero or k % 2 != parity:
                continue
            undetermined[str(k)] = {
                "c_abs_upper": _dec_upper(enc.c[k]),
                "lambda_abs_upper": _dec_upper(enc.lam[k]),
                "derivative_abs_upper": _dec_upper(enc.c[k] * math.factorial(k)),
                "meaning": f"rigorous: |L^({k})(E,1)/{k}!| <= c_abs_upper; a ball cannot prove L^({k})(E,1) = 0",
            }
    if certified:
        r = candidates[0]
        if r_upper in c_nonzero:
            how_nz = "ball for c_%d excludes 0" % r_upper
        elif r_upper in lam_nonzero:
            how_nz = "ball for lambda_%d excludes 0" % r_upper
        else:
            how_nz = "modular symbols: L(E,1)/Omega^+ is a non-zero rational"
        reasons = "; ".join(f"c_{k} = 0 ({c_zero[k]})" for k in range(r)) if r else ""
        if r == r_upper:
            conclusion = (f"ord_(s=1) L(E,s) = {r} rigorously: " + (reasons + "; " if reasons else "")
                          + f"c_{r} != 0 ({how_nz}).")
        else:
            conclusion = (f"ord_(s=1) L(E,s) = {r} rigorously: " + (reasons + "; " if reasons else "")
                          + f"ord <= {r_upper} ({how_nz}); ord = {parity} mod 2 (root number {w:+d}) "
                          + f"leaves only {r}.")
    elif r_upper is not None:
        conclusion = (f"rigorous: {r_lower} <= ord <= {r_upper} and ord = {parity} mod 2, so ord in {candidates}. "
                      + " ".join(f"|c_{k}| <= {v['c_abs_upper']} (not proven zero)." for k, v in undetermined.items()))
        if w == -1 and candidates == [1, 3]:
            conclusion += (" ord = 3 iff L'(E,1) = 0; exact vanishing needs another argument "
                           "(e.g. Gross-Zagier: the Heegner point is torsion).")
    else:
        conclusion = (f"rigorous: ord >= {r_lower}; no coefficient up to k = {kmax} is certified non-zero "
                      "at this precision.")
    cond = conditioned_coefficients(enc, list(lam_zero))
    out = {
        "label": label,
        "ainvs": [int(v) for v in ainvs],
        "conductor": N,
        "root_number": w,
        "kmax": kmax,
        "prec_bits": prec,
        "M_terms": enc.M,
        "c": [ball_to_json(x, digits) for x in enc.c],
        "lambda": [ball_to_json(x, digits) for x in enc.lam],
        "c_conditioned": [ball_to_json(x, digits) for x in cond],
        "c_conditioned_note": "c_k recomputed with lambda_j := 0 exactly for every j in lambda_zero_exact "
                              "(rigorous given those exact zeros; tighter than 'c')",
        "lambda_zero_exact": {str(k): v for k, v in lam_zero.items()},
        "c_zero_exact": {str(k): v for k, v in c_zero.items()},
        "c_nonzero_certified": c_nonzero,
        "lambda_nonzero_certified": lam_nonzero,
        "parity_zero_lambda": [k for k in range(K) if (-1) ** k == -w],
        "parity_zero_c": [k for k in range(K) if all((-1) ** j == -w for j in range(k + 1))],
        "modular_symbols": modsym,
        "external_inputs": external or None,
        "order_lower": r_lower,
        "order_upper": r_upper,
        "order": candidates[0] if certified else None,
        "order_certified": bool(certified),
        "order_candidates": candidates,
        "undetermined": undetermined,
        "conclusion": conclusion,
        "assumptions": ASSUMPTIONS,
        "timings": dict(enc.timings, certify_total_s=round(time.time() - t0, 4)),
    }
    if w == -1 and K > 1 and "1" in undetermined:
        out["L_prime_abs_upper"] = undetermined["1"]["derivative_abs_upper"]
    return out


def _cli(argv=None):  # pragma: no cover - convenience
    import argparse
    import json
    ap = argparse.ArgumentParser(description="Rigorous enclosures of L^{(k)}(E,1)/k! and certified orders.")
    ap.add_argument("curve", help="Cremona label (e.g. 389a1) or a-invariants like [0,1,1,-2,0]")
    ap.add_argument("--kmax", type=int, default=4)
    ap.add_argument("--prec", type=int, default=DEFAULT_PREC)
    args = ap.parse_args(argv)
    if args.curve.startswith("["):
        ainvs = [int(t) for t in args.curve.strip("[]").split(",")]
        label = None
    else:
        from .ecdata import by_label
        ainvs = by_label(args.curve).ainvs
        label = args.curve
    res = certify_order_rigorous(ainvs, kmax=args.kmax, prec=args.prec, label=label)
    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":  # pragma: no cover
    _cli()
