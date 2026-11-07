"""Certificate that ord_{s=1} L(E,s) >= 3 by the Gross-Zagier formula.

This is the method of J. Buhler, B. Gross and D. Zagier, "On the conjecture of
Birch and Swinnerton-Dyer for an elliptic curve of rank 3", Math. Comp. 44
(1985) 473-481, who proved ord_{s=1} L(E,s) = 3 for 5077a1 (y^2 + y = x^3 - 7x + 6).
Nothing here is new mathematics; the module mechanises their argument and makes
every numerical step rigorous with ball (interval) arithmetic (Arb, via
python-flint).  It proves nothing about BSD beyond the order of vanishing.

Logical chain, for an optimal curve E/Q of conductor N with root number w(E) = -1
---------------------------------------------------------------------------------
 (1) [exact]  w(E) = -1, so Lambda(E,1+t) is odd in t: L(E,1) = 0 and ord_{s=1} L(E,s) is odd.
 (2) [exact]  Pick an odd fundamental discriminant D < -4, coprime to N, such that every
     prime p | N splits in K = Q(sqrt D) (equivalently D is a square mod 4N).  Then
     w(E^D) = +1 and L(E/K, s) = L(E, s) L(E^D, s).
 (3) [exact]  S_D := sum_{0<=a<|D|} (D|a) x^-({oo -> a/|D|}) in Q, where x^- is the
     minus modular symbol of E (PARI msfromell(E,-1)).  By the Birch-Manin twisted
     modular-symbol formula  S_D = L(E^D,1) sqrt|D| / Omega^-  (we only use that S_D is
     a non-zero multiple of L(E^D,1)), so S_D != 0  =>  L(E^D,1) != 0.
 (4) [theorem, Gross-Zagier 1986]  L'(E/K,1) = L'(E,1) L(E^D,1) (as L(E,1) = 0) and
         L'(E/K,1) = 2 Vol(E) / (c_E^2 u_K^2 sqrt|D|) * hhat_K(y_K),
     where y_K in E(K) is the Heegner point (trace to K of phi(x_D)), Vol(E) is the area of
     the period lattice of the Neron differential, c_E the Manin constant, u_K = 1 for D < -4,
     and hhat_K = 2 hhat_Q on E(Q) (hhat_Q = canonical height in the normalisation of
     Cremona / PARI ellheight / Sage).  The normalisation was checked numerically on 120
     (E, D) pairs of rank one, where hhat_K(y_K)/(2 Reg) must be a perfect square.
 (5) [theorem]  y_K - tau(y_K) is torsion (tau = complex conjugation; GZ / Gross: tau(y_K) =
     -w(E) y_K + torsion; independently, rank E^D(Q) = 0 by Kolyvagin or Kato since
     L(E^D,1) != 0).  By the parallelogram law hhat_K(y_K) = hhat_K(y_K + tau y_K)/4 =
     hhat_Q(Q)/2 with Q = y_K + tau y_K in E(Q).  Hence either hhat_K(y_K) = 0 or
     hhat_K(y_K) >= lambda/2, where lambda is a lower bound for hhat_Q on non-torsion points.
 (6) [interval arithmetic]  Rigorous enclosures give |L'(E,1)| <= eps, |L(E^D,1)| <= B_D,
     Vol(E) >= V.  If H := sqrt|D| c_E^2 eps B_D / (2 V) < lambda/2 then hhat_K(y_K) = 0, so
     L'(E,1) L(E^D,1) = 0, so L'(E,1) = 0 and (being odd) ord_{s=1} L(E,s) >= 3.
 (7) [interval arithmetic]  An enclosure of L'''(E,1)/3! that excludes 0 gives ord = 3.

What is rigorous: all of (1)-(3) are exact computations (PARI, rational arithmetic);
(6)-(7) are ball-arithmetic evaluations of the standard rapidly convergent series with
proven tail bounds (|a_n| <= d(n) sqrt n <= 2n), exact integer coefficients a_n (PARI
ellan) and rigorous root isolation + AGM for the periods.  The height lower bound
lambda comes from Sage's implementation of Cremona-Siksek / Thongjunthug (documented
as rigorous) or, as a cross-check, from the Mordell-Weil lattice.  Trusted software:
PARI/GP, FLINT/Arb, Sage.  Trusted theorems: modularity, Gross-Zagier, the action of
complex conjugation on y_K (or Kolyvagin/Kato), Cesnavicius's c_E = 1 for semistable
optimal curves, the twisted modular-symbol formula.
"""
from __future__ import annotations

import math
import time
from fractions import Fraction
from typing import Any, Iterable, Sequence

from .pari_session import get_pari, time_limit

try:  # python-flint provides Arb ball arithmetic
    from flint import arb, ctx as _flint_ctx, fmpz_poly
except ImportError:  # pragma: no cover
    arb = None
    _flint_ctx = None
    fmpz_poly = None

LN2 = math.log(2.0)
LN10 = math.log(10.0)

REFERENCES = {
    "BGZ85": "J. Buhler, B. H. Gross, D. B. Zagier, On the conjecture of Birch and Swinnerton-Dyer for an "
             "elliptic curve of rank 3, Math. Comp. 44 (1985), no. 170, 473-481.",
    "GZ86": "B. H. Gross, D. B. Zagier, Heegner points and derivatives of L-series, Invent. Math. 84 (1986), "
            "225-320 (Thm. I.6.3 and its elliptic-curve form in Ch. V).",
    "Ces18": "K. Cesnavicius, The Manin constant in the semistable case, Compos. Math. 154 (2018), 1889-1920.",
    "CNS24": "K. Cesnavicius, M. Neururer, A. Saha, The Manin constant and the modular degree, "
             "J. Eur. Math. Soc. 26 (2024), 573-637 (c | 6 deg(phi); c | deg(phi) if N is cube-free).",
    "Kol88": "V. A. Kolyvagin, Finiteness of E(Q) and Sha(E,Q) for a subclass of Weil curves, Izv. Akad. "
             "Nauk SSSR 52 (1988) (with GZ86 and Bump-Friedberg-Hoffstein / Murty-Murty).",
    "Kato04": "K. Kato, p-adic Hodge theory and values of zeta functions of modular forms, Asterisque 295 (2004).",
    "Gross91": "B. H. Gross, Kolyvagin's work on modular elliptic curves, in: L-functions and arithmetic "
               "(Durham 1989), LMS Lecture Note Ser. 153 (1991), 235-256.",
    "CS06": "J. E. Cremona, S. Siksek, Computing a lower bound for the canonical height on elliptic curves "
            "over Q, ANTS VII, LNCS 4076 (2006), 275-286.",
    "Tho10": "T. Thongjunthug, Heights on elliptic curves over number fields, period lattices, and complex "
             "elliptic logarithms, PhD thesis, Warwick (2010).",
    "BCDT01": "C. Breuil, B. Conrad, F. Diamond, R. Taylor, On the modularity of elliptic curves over Q, "
              "J. Amer. Math. Soc. 14 (2001), 843-939 (semistable case: Wiles 1995, Taylor-Wiles 1995).",
}


class GZCertificateError(ValueError):
    """Raised when the certificate does not apply to the given curve."""


# ----------------------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------------------
class _prec:
    """Context manager setting the Arb working precision (bits)."""

    def __init__(self, bits: int):
        self.bits = int(bits)

    def __enter__(self):
        self.old = _flint_ctx.prec
        _flint_ctx.prec = self.bits
        return self

    def __exit__(self, *exc):
        _flint_ctx.prec = self.old
        return False


def _need_flint():
    if arb is None:  # pragma: no cover
        raise RuntimeError("python-flint (Arb) is required for the rigorous enclosures")


def _log10_upper(x) -> float:
    """A float >= log10(x) for a positive arb ball x (uses the upper endpoint, then rounds up)."""
    u = x.upper()
    m, e = u.man_exp()
    m = int(m)
    e = int(e)
    if m <= 0:
        return float("-inf")
    return (math.log10(m) + e * math.log10(2.0)) + 1e-12


def _log10_lower(x) -> float:
    """A float <= log10(x) for a positive arb ball x (uses the lower endpoint, then rounds down)."""
    lo = x.lower()
    m, e = lo.man_exp()
    m = int(m)
    e = int(e)
    if m <= 0:
        return float("-inf")
    return (math.log10(m) + e * math.log10(2.0)) - 1e-12


def ball_dict(x, digits: int = 25) -> dict:
    """JSON-friendly description of an arb ball."""
    return {
        "ball": x.str(digits),
        "mid": x.mid().str(digits, radius=False),
        "rad": x.rad().str(5, radius=False),
        "contains_zero": bool(x.contains(0)),
    }


def _upper_float(x) -> float:
    """A float that is >= every point of the ball x (rounded outward)."""
    v = float(x.upper())
    return math.nextafter(v, math.inf)


def _lower_float(x) -> float:
    v = float(x.lower())
    return math.nextafter(v, -math.inf)


# ----------------------------------------------------------------------------------------
# exact arithmetic data of the curve
# ----------------------------------------------------------------------------------------
def curve_data(ainvs: Sequence[int]) -> dict:
    """Exact invariants needed by the certificate (minimal model, conductor, root number, ...)."""
    pari = get_pari()
    E0 = pari.ellinit([int(a) for a in ainvs])
    Emin = pari.ellminimalmodel(E0)
    amin = [int(a) for a in Emin[:5]]
    E = pari.ellinit(amin)
    N = int(pari.ellglobalred(E)[0])
    fa = pari.factor(N)
    fac = [[int(p), int(e)] for p, e in zip(fa[0], fa[1])]
    w = int(pari.ellrootno(E))
    disc = int(pari("(E)->E.disc")(E))
    tors = [int(t) for t in pari.elltors(E)[1]]
    try:
        isog = pari.ellisomat(E)
        class_size = len(isog[0])
        isog_degrees = sorted({int(x) for col in pari.Vec(isog[1]) for x in col})
    except Exception as exc:  # noqa: BLE001
        class_size, isog_degrees = None, f"ellisomat failed: {exc}"
    return dict(ainvs_input=[int(a) for a in ainvs], ainvs=amin, conductor=N, factorization=fac,
                root_number=w, discriminant=disc, torsion_structure=tors, isogeny_class_size=class_size,
                isogeny_degrees_from_E=isog_degrees, semistable=all(e == 1 for _, e in fac), _E=E)


def is_heegner_discriminant(D: int, N: int, odd_only: bool = True) -> bool:
    """D < -4 fundamental, (odd), coprime to N and a square mod 4N (all p | N split in Q(sqrt D))."""
    pari = get_pari()
    if D >= -4:
        return False
    if odd_only and D % 2 == 0:
        return False
    if int(pari.isfundamental(D)) != 1 or math.gcd(D, N) != 1:
        return False
    return int(pari.issquare(pari.Mod(D, 4 * N))) == 1


def heegner_discriminants(N: int, count: int = 10, odd_only: bool = True, max_abs: int = 10**6) -> list[int]:
    """The first `count` discriminants D < -4 (odd if `odd_only`) satisfying the Heegner hypothesis."""
    out: list[int] = []
    D = -4
    while len(out) < count and -D < max_abs:
        D -= 1
        if is_heegner_discriminant(D, N, odd_only):
            out.append(D)
    return out


def minus_modular_symbol(E, timeout: float | None = 3600):
    """[M, x^-] = msfromell(E, -1): the minus modular symbol of E (exact, over Q)."""
    pari = get_pari()
    with time_limit(timeout):
        M, xm = pari.msfromell(E, -1)
    return M, xm


def twisted_minus_symbol_sum(M, xm, D: int) -> Fraction:
    """S_D = sum_{0<=a<|D|} (D|a) x^-({oo -> a/|D|}), an exact rational.

    Numerically S_D = L(E^D,1) sqrt|D| / Omega^- (Omega^- as in PARI's msfromell documentation;
    note the factor sqrt|D|, which PARI 2.17.2's help text omits).  Only S_D != 0 is used.
    """
    pari = get_pari()
    m = -D
    total = pari(0)
    for a in range(1, m):
        k = int(pari.kronecker(D, a))
        if k:
            total += k * pari.mseval(M, xm, pari(f"[oo, {a}/{m}]"))
    return Fraction(str(total))


def an_list(E, M: int) -> list[int]:
    """a_1..a_M of L(E,s) (exact integers, PARI ellan)."""
    pari = get_pari()
    v = pari.ellan(E, M)
    return [int(x) for x in v]


# ----------------------------------------------------------------------------------------
# rigorous enclosures (Arb)
# ----------------------------------------------------------------------------------------
def _tail_bound(alpha, M: int, k: int):
    """Upper bound (arb) for sum_{n>M} 2 e^{-n alpha} / (n alpha)^k  (k = 0, 1, 3, ...).

    Used with |a_n / n| <= 2 (from |a_n| <= d(n) sqrt(n) and d(n) <= 2 sqrt(n)) and
    0 <= G_k(x) <= e^{-x}/x^k for the kernels used here.  For n > M,
    (n alpha)^{-k} <= ((M+1) alpha)^{-k}, and sum_{n>M} e^{-n alpha} = e^{-(M+1) alpha}/(1-e^{-alpha}).
    """
    x0 = alpha * (M + 1)
    geo = (-x0).exp() / (1 - (-alpha).exp())
    return 2 * geo / (x0 ** k) if k else 2 * geo


def _choose_terms(alpha_float: float, log10_target: float, k: int) -> int:
    """Smallest M (approximately) with tail bound <= 10^log10_target; verified later in Arb."""
    a = alpha_float
    # 2 e^{-x0} / (x0^k (1-e^{-a})) <= T,  x0 = (M+1) a
    lnT = log10_target * LN10
    one_minus = -math.expm1(-a)
    x0 = max(1.0, -lnT)
    for _ in range(60):
        x0_new = math.log(2.0 / one_minus) - lnT - k * math.log(max(x0, 1e-300))
        if abs(x0_new - x0) < 1e-9:
            break
        x0 = max(x0_new, 1e-6)
    return max(1, int(math.ceil(x0 / a)) + 2)


def lprime_ball(an: Sequence[int], N: int, log10_target: float = -60.0, prec: int | None = None,
                M: int | None = None, progress=None) -> dict:
    """Rigorous enclosure of L'(E,1), assuming root number -1.

        L'(E,1) = 2 sum_{n>=1} (a_n/n) E_1(2 pi n / sqrt N)          [w(E) = -1]

    The truncation error after M terms is at most 2 * tail_bound(alpha, M, 1) (|a_n/n| <= 2,
    E_1(x) <= e^{-x}/x).  Returns the ball together with the bookkeeping.
    """
    _need_flint()
    alpha_f = 2 * math.pi / math.sqrt(N)
    if M is None:
        M = _choose_terms(alpha_f, log10_target, 1)
    if len(an) < M:
        raise ValueError(f"need {M} coefficients a_n, got {len(an)}")
    if prec is None:
        # Arb's expint enclosures can be ~60-110 bits looser than the working precision (measured:
        # absolute radius up to ~2^(60-prec) for 1 < x < 150), hence the generous guard bits.
        prec = int(math.ceil(-log10_target * LN10 / LN2)) + 110 + int(math.log2(M + 1))
    t0 = time.time()
    with _prec(prec):
        alpha = 2 * arb.pi() / arb(N).sqrt()
        S = arb(0)
        for n in range(1, M + 1):
            a = an[n - 1]
            if a:
                S += arb(a) * (alpha * n).expint(1) / n
            if progress is not None and n % 200000 == 0:
                progress(n, S)
        tail = _tail_bound(alpha, M, 1)
        value = 2 * (S + arb(0, 1) * tail.upper())
        out = dict(value=value, ball=ball_dict(value), M=M, prec_bits=prec,
                   tail_bound_log10=_log10_upper(2 * tail), upper_log10=_log10_upper(abs(value)),
                   formula="L'(E,1) = 2 sum a_n/n E_1(2 pi n/sqrt N)", seconds=round(time.time() - t0, 2))
    return out


def twist_value_ball(an: Sequence[int], N: int, D: int, log10_target: float = -20.0, prec: int | None = None,
                     M: int | None = None) -> dict:
    """Rigorous enclosure of L(E^D,1) for D coprime to N with w(E^D) = +1.

        L(E^D,1) = 2 sum_{n>=1} chi_D(n) a_n / n * exp(-2 pi n / (|D| sqrt N))

    (conductor of E^D is N D^2 for odd fundamental D coprime to N; a_n(E^D) = chi_D(n) a_n(E)).
    """
    _need_flint()
    pari = get_pari()
    m = abs(D)
    chi = [int(pari.kronecker(D, a)) for a in range(m)]          # chi_D is periodic mod |D| (D = 1 mod 4)
    beta_f = 2 * math.pi / (m * math.sqrt(N))
    if M is None:
        M = _choose_terms(beta_f, log10_target, 0)
    if len(an) < M:
        raise ValueError(f"need {M} coefficients a_n, got {len(an)}")
    if prec is None:
        prec = int(math.ceil(-log10_target * LN10 / LN2)) + 40 + 2 * int(math.log2(M + 1))
    t0 = time.time()
    with _prec(prec):
        beta = 2 * arb.pi() / (arb(N).sqrt() * m)
        q = (-beta).exp()
        qn = arb(1)
        S = arb(0)
        for n in range(1, M + 1):
            qn *= q
            c = chi[n % m] * an[n - 1]
            if c:
                S += qn * c / n
        tail = _tail_bound(beta, M, 0)
        value = 2 * (S + arb(0, 1) * tail.upper())
        out = dict(value=value, ball=ball_dict(value), M=M, prec_bits=prec, conductor_twist=N * m * m,
                   tail_bound_log10=_log10_upper(2 * tail),
                   formula="L(E^D,1) = 2 sum chi_D(n) a_n/n exp(-2 pi n/(|D| sqrt N))",
                   seconds=round(time.time() - t0, 2))
    return out


def _G3(x, prec: int):
    """Rigorous ball for G_3(x) = (1/2) int_1^oo e^{-xu} (log u)^2 du/u, x > 0 (an arb).

    Series: G_3(x) = t^3/6 + pi^2 t/12 - zeta(3)/3 - sum_{k>=1} (-x)^k / (k^3 k!),  t = -log x - gamma,
    with the alternating tail after K terms bounded by |term_{K+1}| / (1 - x/(K+2)) (K + 2 > x).
    """
    t = -x.log() - arb.const_euler()
    P = t ** 3 / 6 + arb.pi() ** 2 * t / 12 - arb(3).zeta() / 3
    S = arb(0)
    term = arb(1)
    xf = float(x.mid())
    k = 0
    eps = arb(2) ** (-prec)
    while True:
        k += 1
        term = term * (-x) / k                      # (-x)^k / k!
        S += term / (k ** 3)
        if k > 2 * xf + 10:
            bound = abs(term) * x / (k + 1) / (1 - x / (k + 2)) / ((k + 1) ** 3)
            if bound < eps:
                S += arb(0, 1) * bound.upper()
                break
    return P - S


def third_derivative_ball(an: Sequence[int], N: int, log10_target: float = -30.0, M: int | None = None) -> dict:
    """Rigorous enclosure of L'''(E,1)/3!, valid once L(E,1) = L'(E,1) = L''(E,1) = 0 is known.

        L'''(E,1)/3! = 2 sum_{n>=1} (a_n/n) G_3(2 pi n / sqrt N),  0 <= G_3(x) <= e^{-x}/x^3.

    (In general this sum is Lambda'''(1)/(3! * sqrt N/(2 pi)); it equals L'''(1)/3! because the lower
    derivatives of L vanish.)
    """
    _need_flint()
    alpha_f = 2 * math.pi / math.sqrt(N)
    if M is None:
        M = _choose_terms(alpha_f, log10_target, 3)
    if len(an) < M:
        raise ValueError(f"need {M} coefficients a_n, got {len(an)}")
    xmax = alpha_f * M
    prec = int(math.ceil(-log10_target * LN10 / LN2 + xmax / LN2)) + 60 + int(math.log2(M + 1))
    t0 = time.time()
    with _prec(prec):
        alpha = 2 * arb.pi() / arb(N).sqrt()
        S = arb(0)
        for n in range(1, M + 1):
            a = an[n - 1]
            if a:
                S += arb(a) * _G3(alpha * n, prec) / n
        tail = _tail_bound(alpha, M, 3)
        value = 2 * (S + arb(0, 1) * tail.upper())
        out = dict(value=value, ball=ball_dict(value), M=M, prec_bits=prec,
                   tail_bound_log10=_log10_upper(2 * tail),
                   formula="L'''(E,1)/3! = 2 sum a_n/n G_3(2 pi n/sqrt N)", seconds=round(time.time() - t0, 2))
    return out


def period_area_ball(ainvs_min: Sequence[int], prec: int = 128) -> dict:
    """Rigorous enclosure of Vol(E) = area of C/Lambda, Lambda the period lattice of the Neron
    differential dx/(2y + a1 x + a3) of the given (minimal) model.

    Roots e_i of 4x^3 + b2 x^2 + 2 b4 x + b6 are isolated rigorously (FLINT/Arb), then
      Delta > 0:  w1 = pi/AGM(sqrt(e1-e3), sqrt(e1-e2)),  w2/i = pi/AGM(sqrt(e1-e3), sqrt(e2-e3)), Vol = w1 |w2|;
      Delta < 0:  b = sqrt(3e1^2 + b2 e1/2 + b4/2), a = 3e1 + b2/4,
                  w1 = 2 pi/AGM(2 sqrt b, sqrt(2b + a)),  Im w2 = pi/AGM(2 sqrt b, sqrt(2b - a)),  Vol = w1 Im w2.
    (Cohen, A Course in Computational Algebraic Number Theory, Alg. 7.4.7; all AGMs are of positive reals.)
    """
    _need_flint()
    a1, a2, a3, a4, a6 = (int(a) for a in ainvs_min)
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    disc = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    with _prec(prec):
        f = fmpz_poly([b6, 2 * b4, b2, 4])
        roots = [r for r, mult in f.complex_roots()]
        pi = arb.pi()
        if disc > 0:
            if len(roots) != 3:
                raise RuntimeError("root isolation failed")
            e = sorted((r.real for r in roots), key=lambda z: float(z.mid()), reverse=True)
            e1, e2, e3 = e
            if not (e1 > e2 and e2 > e3):
                raise RuntimeError("real roots not separated; raise prec")
            w1 = pi / arb.agm((e1 - e3).sqrt(), (e1 - e2).sqrt())
            w2 = pi / arb.agm((e1 - e3).sqrt(), (e2 - e3).sqrt())
            vol = w1 * w2
            components = 2
        else:
            real = [r for r in roots if r.imag.contains(0)]
            if len(real) != 1:
                raise RuntimeError("could not identify the real root; raise prec")
            e1 = real[0].real
            bb = (3 * e1 * e1 + arb(b2) / 2 * e1 + arb(b4) / 2).sqrt()
            aa = 3 * e1 + arb(b2) / 4
            w1 = 2 * pi / arb.agm(2 * bb.sqrt(), (2 * bb + aa).sqrt())
            y2 = pi / arb.agm(2 * bb.sqrt(), (2 * bb - aa).sqrt())
            vol = w1 * y2
            components = 1
        return dict(value=vol, ball=ball_dict(vol), real_period=ball_dict(w1), real_components=components,
                    discriminant_sign=1 if disc > 0 else -1)


# ----------------------------------------------------------------------------------------
# lower bound for canonical heights of non-torsion points of E(Q)
# ----------------------------------------------------------------------------------------
def height_lower_bound_sage(ainvs: Sequence[int], tol: float = 1e-3, n_max: int = 20) -> dict:
    """Sage's E.height_function().min(tol, n_max): documented as a rigorous lower bound for the
    canonical height of all points of infinite order on E(Q) (Cremona-Siksek 2006 / Thongjunthug
    2010, interval arithmetic in the intersection tests, divided by the square of the lcm of the
    component-group exponents incl. the archimedean one).  Normalisation: Cremona's (= PARI ellheight)."""
    import sage.all as S  # noqa: WPS433  (lazy, ~10 s)
    E = S.EllipticCurve(S.QQ, [int(a) for a in ainvs])
    t0 = time.time()
    mu = float(E.height_function().min(tol, n_max))
    return dict(value=mu * (1 - 1e-12), raw=mu, method="sage height_function().min (Cremona-Siksek / Thongjunthug)",
                tol=tol, n_max=n_max, seconds=round(time.time() - t0, 2))


def height_lattice_minimum(E, gens: Iterable[Sequence], ) -> dict:
    """Minimum of the canonical height on the lattice spanned by `gens` (Fincke-Pohst, PARI qfminim).

    This equals min{hhat(P): P in E(Q) non-torsion} only if `gens` is a basis of E(Q)/torsion (needs a
    rank proof and a saturation proof).  Floating point, no interval arithmetic: used as a cross-check."""
    pari = get_pari()
    pts = [pari([pari(str(x)), pari(str(y))]) for x, y in gens]
    G = pari.ellheightmatrix(E, pari(pts))
    r = len(pts)
    if r == 0:
        return dict(value=None, note="rank 0")
    res = pari.qfminim(G, None, None, 2)
    mn = float(res[1])
    vec = [int(t) for t in res[2][0]] if int(res[0]) > 0 else None
    return dict(value=mn, minimal_vector=vec, gram_det=float(pari.matdet(G)),
                method="min of the height Gram matrix of the given generators (PARI qfminim, flag 2)")


# ----------------------------------------------------------------------------------------
# Manin constant
# ----------------------------------------------------------------------------------------
def manin_constant_status(cd: dict, manin_bound: float | None = None) -> dict:
    if manin_bound is not None:
        return dict(c_max=float(manin_bound), status="ASSUMED (supplied by caller)", rigorous=False)
    if cd["semistable"] and cd["isogeny_class_size"] == 1:
        return dict(c_max=1.0, rigorous=True,
                    status="c_E = 1: E is the only curve in its isogeny class (hence the optimal quotient of "
                           "J_0(N)) and N is squarefree (semistable), so Cesnavicius (Compos. Math. 2018) applies.")
    if cd["semistable"]:
        return dict(c_max=None, rigorous=False,
                    status="semistable, but E is not known to be optimal here (isogeny class has %s curves); "
                           "c_E = 1 holds for the optimal curve only" % cd["isogeny_class_size"])
    return dict(c_max=None, rigorous=False,
                status="N is not squarefree: no explicit unconditional bound for c_E is used here "
                       "(Manin's conjecture predicts c_E = 1 for the optimal curve; Cesnavicius-Neururer-Saha "
                       "give c | 6 deg(phi), not explicit without the modular degree)")


# ----------------------------------------------------------------------------------------
# modular-symbol fallback through Sage/eclib (cross-check only)
# ----------------------------------------------------------------------------------------
def twisted_minus_symbols_eclib(ainvs: Sequence[int], Ds: Sequence[int]) -> dict:
    """S_D for each D via Sage's eclib minus modular symbol (Cremona's C++ code).

    eclib does its linear algebra modulo a large prime and reconstructs rationals; we use its output
    only to choose D and as a cross-check, never as the sole proof of L(E^D,1) != 0.
    """
    import sage.all as S  # noqa: WPS433
    E = S.EllipticCurve(S.QQ, [int(a) for a in ainvs])
    t0 = time.time()
    ms = E.modular_symbol(sign=-1, implementation="eclib")
    out = {}
    for D in Ds:
        m = -int(D)
        out[int(D)] = str(sum(S.kronecker(D, a) * ms(S.QQ(a) / m) for a in range(1, m)))
    return dict(values=out, seconds=round(time.time() - t0, 2),
                method="Sage E.modular_symbol(sign=-1, implementation='eclib')")


# ----------------------------------------------------------------------------------------
# the certificate
# ----------------------------------------------------------------------------------------
def certify_rank3(ainvs: Sequence[int], D: int | None = None, *, candidate_D: int = 12,
                  log10_target_lprime: float = -60.0, log10_target_twist: float = -20.0,
                  log10_target_third: float = -30.0, height_bound: Any = "sage", gens=None,
                  manin_bound: float | None = None, safety: float = 1e6, want_third: bool = True,
                  modsym: bool = True, modsym_timeout: float | None = 3600, modsym_fallback: str | None = "eclib",
                  label: str | None = None) -> dict:
    """Try to prove ord_{s=1} L(E,s) >= 3 (and = 3) for E = EllipticCurve(ainvs) with root number -1.

    Parameters
    ----------
    D : Heegner discriminant to use (default: the smallest odd D < -4 satisfying the Heegner hypothesis
        for which L(E^D,1) != 0 is established, among the first `candidate_D`).
    height_bound : "sage" (Cremona-Siksek lower bound via Sage), "lattice" (minimum of the MW lattice of
        `gens`; rigorous only if `gens` is a basis), or a float lambda supplied by the caller.
    modsym : "eclib" -> skip PARI and use Sage/eclib symbols only to choose D (non-vanishing then rests on
        the ball enclosure).  True -> try to certify L(E^D,1) != 0 exactly with PARI's modular symbols (feasible for N up to
        roughly 10^4; PARI's msinit can overflow a 2 GB stack at some levels).  Whatever happens, the
        rigorous ball enclosure of L(E^D,1) is computed and proves L(E^D,1) != 0 whenever it excludes 0.
    modsym_fallback : "eclib" -> if PARI fails, compute S_D with Sage/eclib (used to choose D and as a
        cross-check only); None -> no fallback.
    safety : the Gross-Zagier inequality must hold with this extra factor (absorbs any slip in the
        normalisation of the Gross-Zagier constant; the normalisation itself was checked numerically).

    Returns a JSON-serialisable dict; ``conclusion["certified"]`` is True iff every step is rigorous.
    """
    _need_flint()
    pari = get_pari()
    T: dict[str, float] = {}
    t0 = time.time()
    cd = curve_data(ainvs)
    E = cd.pop("_E")
    N = cd["conductor"]
    T["curve_data"] = round(time.time() - t0, 2)
    cert: dict[str, Any] = dict(label=label, curve=cd)
    if cd["root_number"] != -1:
        raise GZCertificateError(f"root number of E is {cd['root_number']}: this certificate needs w(E) = -1")
    if D is not None and not is_heegner_discriminant(D, N):
        raise GZCertificateError(f"D = {D} is not an odd fundamental discriminant < -4 satisfying the Heegner hypothesis")
    cands = [int(D)] if D is not None else heegner_discriminants(N, candidate_D)

    # ---- step 3: modular symbols (exact) -- optional
    t = time.time()
    sym_exact: dict[int, Fraction] = {}
    ms_status = "not attempted"
    if modsym and modsym != "eclib":
        try:
            M_ms, xm = minus_modular_symbol(E, timeout=modsym_timeout)
            for d in cands:
                sym_exact[d] = twisted_minus_symbol_sum(M_ms, xm, d)
            ms_status = "PARI msfromell(E,-1) + mseval (exact rational arithmetic)"
        except Exception as exc:  # noqa: BLE001  (PARI stack overflow, timeout)
            ms_status = f"PARI modular symbols failed: {type(exc).__name__}: {str(exc)[:120]}"
    T["modular_symbols"] = round(time.time() - t, 2)
    sym_eclib = None
    if modsym == "eclib":
        ms_status = "PARI modular symbols skipped (caller); eclib used to choose D"
    if not sym_exact and modsym and (modsym == "eclib" or modsym_fallback == "eclib"):
        try:
            sym_eclib = twisted_minus_symbols_eclib(cd["ainvs"], cands)
        except Exception as exc:  # noqa: BLE001
            sym_eclib = dict(error=f"{type(exc).__name__}: {str(exc)[:160]}")
    cert["modular_symbols"] = dict(status=ms_status,
                                   exact=[dict(D=d, S_D=str(s)) for d, s in sym_exact.items()],
                                   eclib_cross_check=sym_eclib)
    # choose D: first candidate whose symbol is non-zero (exact, else eclib), else the first candidate
    D_used = None
    for d in cands:
        if d in sym_exact:
            if sym_exact[d] != 0:
                D_used = d
                break
        elif sym_eclib and "values" in sym_eclib and sym_eclib["values"].get(d, "0") != "0":
            D_used = d
            break
    if D_used is None:
        if sym_exact and all(v == 0 for v in sym_exact.values()):
            cert["conclusion"] = dict(certified=False, reason="all candidate D have S_D = 0 (L(E^D,1) = 0)")
            return cert
        D_used = cands[0]
    S_D = sym_exact.get(D_used)
    w_twist = int(pari.ellrootno(pari.ellinit(pari.elltwist(E, D_used))))
    cert["heegner"] = dict(D=D_used, class_number=int(pari.qfbclassno(D_used)), u_K=1, odd=True,
                           coprime_to_N=True, square_mod_4N=bool(is_heegner_discriminant(D_used, N)),
                           split_primes=[p for p, _ in cd["factorization"]], twist_root_number=w_twist)
    if w_twist != 1:  # impossible under the Heegner hypothesis (w(E^D) = -w(E)); defensive
        raise GZCertificateError("twist root number is not +1")

    # ---- coefficients
    alpha_f = 2 * math.pi / math.sqrt(N)
    M1 = _choose_terms(alpha_f, log10_target_lprime, 1)
    M3 = _choose_terms(alpha_f, log10_target_third, 3) if want_third else 0
    MD = _choose_terms(2 * math.pi / (abs(D_used) * math.sqrt(N)), log10_target_twist, 0)
    t = time.time()
    an = an_list(E, max(M1, M3, MD))
    T["ellan"] = round(time.time() - t, 2)

    # ---- rigorous enclosures
    lp = lprime_ball(an, N, log10_target_lprime, M=M1)
    T["lprime"] = lp["seconds"]
    tw = twist_value_ball(an, N, D_used, log10_target_twist, M=MD)
    T["twist"] = tw["seconds"]
    vol = period_area_ball(cd["ainvs"])
    om = pari("(E)->E.omega")(E)
    omega_minus = float(-2.0 / vol["real_components"] * float(pari.imag(om[1])))
    twist_nonzero_ball = bool(not tw["value"].contains(0))
    consistency = None
    s_for_check = S_D if S_D is not None else (
        Fraction(sym_eclib["values"][D_used]) if sym_eclib and "values" in sym_eclib else None)
    if s_for_check is not None:
        predicted = float(s_for_check) * omega_minus / math.sqrt(-D_used)
        consistency = dict(symbol=str(s_for_check), predicted_L=predicted,
                           ball_contains_prediction=bool(tw["value"].overlaps(
                               arb(predicted) + arb(0, 1) * arb(abs(predicted) * 1e-12 + 1e-30))),
                           note="S_D * Omega^- / sqrt|D| (PARI periods, floating point) vs the ball for L(E^D,1)")

    # ---- height lower bound
    if isinstance(height_bound, (int, float)):
        hb: dict[str, Any] = dict(value=float(height_bound), method="supplied by caller (caller vouches for it)",
                                  rigorous=None)
    elif height_bound == "sage":
        hb = height_lower_bound_sage(cd["ainvs"])
        hb["rigorous"] = True
    elif height_bound == "lattice":
        if gens is None:
            raise GZCertificateError("height_bound='lattice' needs gens")
        hb = height_lattice_minimum(E, gens)
        hb["rigorous"] = False
        hb["value"] = hb["value"] * 0.999  # floating-point safety
    else:
        raise GZCertificateError(f"unknown height_bound {height_bound!r}")
    if gens is not None and height_bound != "lattice":
        try:
            hb["lattice_cross_check"] = height_lattice_minimum(E, gens)
            lm = hb["lattice_cross_check"]["value"]
            hb["lattice_cross_check"]["consistent"] = bool(lm is None or hb["value"] <= lm * (1 + 1e-9))
        except Exception as exc:  # noqa: BLE001
            hb["lattice_cross_check"] = dict(error=str(exc)[:200])
    lam = hb["value"]

    # ---- Manin constant and the Gross-Zagier inequality
    mc = manin_constant_status(cd, manin_bound)
    with _prec(256):
        Lp_up = abs(lp["value"]).upper()
        LD_up = abs(tw["value"]).upper()
        V_lo = vol["value"].lower()
        if not V_lo > 0:
            raise RuntimeError("volume lower bound not positive")
        c2 = arb(mc["c_max"]) ** 2 if mc["c_max"] is not None else arb(1)
        H = arb(-D_used).sqrt() * c2 * Lp_up * LD_up / (2 * V_lo)       # >= hhat_K(y_K)
        thresh = arb(lam) / 2                                            # y_K non-torsion => hhat_K(y_K) >= lam/2
        holds = bool(H * safety < thresh)
        log10_margin = _log10_lower(thresh) - _log10_upper(H)
        log10_cmax = 0.5 * (log10_margin - math.log10(safety)) if mc["c_max"] is None else None
    gz = dict(formula="hhat_K(y_K) = sqrt|D| c_E^2 u_K^2 |L'(E,1)| |L(E^D,1)| / (2 Vol(E))  (Gross-Zagier)",
              H_upper_log10=_log10_upper(H), threshold=lam / 2, log10_margin=log10_margin, safety=safety,
              inequality_holds=holds, uses_manin_bound=mc["c_max"],
              largest_admissible_manin_constant_log10=log10_cmax)
    cert.update(lprime=_strip(lp), twist=_strip(tw), twist_consistency=consistency, volume=_strip(vol),
                height_lower_bound=hb, manin_constant=mc, gross_zagier=gz)

    # ---- how L(E^D,1) != 0 is known
    if S_D is not None and S_D != 0:
        nonvanishing = f"exact: S_D = {S_D} != 0 (PARI modular symbols); ball also excludes 0: {twist_nonzero_ball}"
        nonvanishing_ok = True
    elif twist_nonzero_ball:
        nonvanishing = ("rigorous ball enclosure of L(E^D,1) excludes 0 (interval arithmetic)"
                        + (f"; eclib cross-check S_D = {sym_eclib['values'].get(D_used)}"
                           if sym_eclib and "values" in sym_eclib else ""))
        nonvanishing_ok = True
    else:
        nonvanishing, nonvanishing_ok = "L(E^D,1) != 0 NOT established", False
    cert["twist_nonvanishing"] = dict(ok=nonvanishing_ok, how=nonvanishing)

    conditional = []
    if mc["c_max"] is None:
        conditional.append(f"Manin constant c_E <= 10^{log10_cmax:.1f} (Manin's conjecture: c_E = 1)")
    if hb.get("rigorous") is False:
        conditional.append(f"height lower bound lambda = {lam:.6g} ({hb['method']})")
    rigorous_ok = bool(holds and nonvanishing_ok and not conditional)

    # ---- third derivative
    third = None
    if want_third:
        th = third_derivative_ball(an, N, log10_target_third, M=M3)
        T["third"] = th["seconds"]
        third = _strip(th)
        third["nonzero"] = bool(not th["value"].contains(0))
        third["pari_ellL1_3_over_6"] = float(pari.ellL1(E, 3)) / 6
    cert["third_derivative"] = third

    if rigorous_ok:
        stmt = "L'(E,1) = 0 exactly, hence ord_(s=1) L(E,s) >= 3"
        if third and third["nonzero"]:
            stmt += "; L'''(E,1)/3! lies in %s, so ord = 3" % third["ball"]["ball"]
    elif holds and nonvanishing_ok:
        stmt = "L'(E,1) = 0 and ord >= 3, CONDITIONAL on: " + "; ".join(conditional)
    elif not holds:
        stmt = "not certified: the Gross-Zagier inequality fails (the Heegner point may be non-torsion)"
    else:
        stmt = "not certified: " + nonvanishing
    cert["conclusion"] = dict(certified=rigorous_ok,
                              ord_lower_bound=3 if (holds and nonvanishing_ok) else 1,
                              ord_exact=3 if (rigorous_ok and third and third["nonzero"]) else None,
                              conditional_on=conditional, statement=stmt)
    cert["assumptions"] = (["height lower bound lambda supplied by the caller"]
                           if isinstance(height_bound, (int, float)) else [])
    cert["exact_steps"] = [
        "root number w(E) = -1 (PARI ellrootno, product of local root numbers)",
        f"D = {D_used}: odd fundamental discriminant, coprime to N, square mod 4N (exact)",
        "L(E^D,1) != 0: " + nonvanishing,
        "a_n (PARI ellan, exact integers); isogeny class size %s (PARI ellisomat)" % cd["isogeny_class_size"],
    ]
    cert["rigorous_numerical_steps"] = [
        "ball enclosures (Arb) of L'(E,1), L(E^D,1), L'''(E,1)/3! with proven tail bounds",
        "ball enclosure of Vol(E) via certified root isolation and real AGMs",
        "height lower bound: " + hb.get("method", "?"),
    ]
    cert["timings"] = T
    cert["seconds_total"] = round(time.time() - t0, 2)
    return cert


def _strip(d: dict) -> dict:
    return {k: v for k, v in d.items() if k != "value"}
