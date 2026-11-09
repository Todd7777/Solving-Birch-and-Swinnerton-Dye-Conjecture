"""Certified bounds for Selmer coranks from p-adic L-functions, using PARI only and exact arithmetic.

For an elliptic curve E/Q and an odd prime p of good ordinary reduction this module computes the
cyclotomic p-adic L-function L_p(E,T) of Mazur, Tate and Teitelbaum as an element of Z_p[[T]]
(up to a non-zero rational factor) by exact Riemann sums of modular symbols. It certifies an
upper bound ord_{T=0} L_p(E,T) <= r and turns it, via Kato's divisibility and Mazur's control
theorem, into the bound corank_{Z_p} Sel_{p^oo}(E/Q) <= r. Combined with r independent
rational points, this proves rank E(Q) = r and the finiteness of Sha(E/Q)[p^oo].

This is the method of Stein and Wuthrich (Math. Comp. 82 (2013), 1757-1792). The module
re-implements its finiteness part with PARI's exact modular symbols (msfromell/mseval) and
writes out every step. It does not depend on PARI's overconvergent-symbol code
(mspadicmoments), which bsd/padic.py uses, so the two give independent certificates.

Mathematics
-----------
1. Integral plus symbol. PARI's msfromell(E, 1) gives a modular symbol x in Hom_Gamma(Delta_0, Q)
   (Gamma = Gamma_0(N), Delta_0 = degree-0 divisors on the cusps) that is a T_l-eigenvector with
   eigenvalue a_l for every prime l not dividing N, and U_l-eigenvector for l | N. Put
   m(r) = x({oo -> r}). By Manin's continued-fraction trick every path between cusps is a sum of
   unimodular paths g{0 -> oo}, g in SL_2(Z), and x(g{0 -> oo}) depends only on the class of the
   bottom row (c : d) in P^1(Z/NZ). So the values of x lie in the Z-module spanned by the
   (N+1)prod(1+1/l) Manin-symbol values. Let g > 0 generate it and put m' = m/g. Then m' is
   Z-valued on all of Q. We compute g from all Manin symbols, not from PARI's internal basis.

2. The measure. Let alpha be the unit root of X^2 - a_p X + p in Z_p (p does not divide a_p). For
   k >= 1 and a in Z_p^x put
        mu(a + p^k Z_p) = alpha^{-k} m'(a/p^k) - alpha^{-k-1} m'(a/p^{k-1}).
   The Hecke relation sum_{j<p} m((r+j)/p) + m(pr) = a_p m(r) together with alpha^2 = a_p alpha - p
   gives finite additivity, so mu is a Z_p-valued measure on Z_p^x. Put
        L(T) = int_{Z_p^x} (1+T)^{log_gamma <x>} dmu(x),     gamma = 1 + p,  <x> = x / omega(x),
   where omega is the Teichmueller character. L(T) is a non-zero rational multiple of the
   Mazur-Swinnerton-Dyer p-adic L-function attached to the trivial Teichmueller component, which
   is the component that controls Sel_{p^oo}(E/Q_oo), Q_oo the cyclotomic Z_p-extension.
   Its constant term is exactly (1 - 1/alpha)^2 m'(0).

3. Riemann sums and their error (proof). Z_p^x is the disjoint union of the balls
   B_{t,k} = omega(t) gamma^k + p^n Z_p, 1 <= t <= p-1, 0 <= k < p^{n-1}. On B_{t,k} one has
   log_gamma <x> = k + p^{n-1} u(x) with u(x) in Z_p. Put
        L^{(n)}(T) = sum_{t,k} mu(B_{t,k}) (1+T)^k.
   The T^j coefficient of L - L^{(n)} is sum_{t,k} int_{B_{t,k}} [binom(k + p^{n-1}u, j) - binom(k, j)] dmu.
   By Vandermonde, binom(k + z, j) - binom(k, j) = sum_{i=1}^{j} binom(z, i) binom(k, j-i), and
   binom(z, i) = (z/i) binom(z-1, i-1) has valuation >= v_p(z) - v_p(i) for z in Z_p. With
   z = p^{n-1}u the integrand lies in p^{e} Z_p, e = n - 1 - floor(log_p j). Integrating a
   p^e Z_p-valued continuous function against a Z_p-valued measure stays in p^e Z_p (closed), so
        coeff_j(L) == coeff_j(L^{(n)})  (mod p^{n-1-floor(log_p j)}),   j >= 1,      (*)
   and the constant terms agree exactly. CERTIFICATE: if v_p(coeff_r(L^{(n)})) < n-1-floor(log_p r),
   then coeff_r(L) != 0, hence ord_{T=0} L <= r.

Theorems turning this into a Selmer bound (p odd, good ordinary, E without CM)
----------------------------------------------------------------------------
Stated as in Stein-Wuthrich (Math. Comp. 82 (2013), sections 6-7), whose running hypotheses are
p != 2 and E without complex multiplication; read there, not in Kato's paper itself:
* Kato [Kat04, Thm 17.4], as SW Thm 7.3: if p is semistable for E and rho_p is surjective, then
  f_E(T) divides L_p(E,T) in Lambda, where f_E is a characteristic series of X_oo = Sel_{p^oo}(E/Q_oo)^dual
  (a torsion Lambda-module, Kato). As SW Thm 7.4: if rhobar_p is not surjective, f_E divides
  p^m L_p(E,T) for some m >= 0. For p >= 5, rhobar_p is surjective iff rho_p is (SW Prop 7.2,
  Serre), so for p >= 5 one of the two always applies. For p = 3 the case "rhobar_3 surjective,
  rho_3 not surjective" is covered by neither statement, so results at p = 3 carry that hypothesis.
* Schneider, Perrin-Riou (SW Thm 6.1): ord_{T=0} f_E >= rank E(Q), with equality iff the cyclotomic
  p-adic height pairing on E(Q) is non-degenerate AND Sha(E/Q)[p^oo] is finite. (The inequality
  corank Sel_{p^oo}(E/Q) <= ord_T f_E is Mazur's control theorem, Invent. Math. 18 (1972).)
Together: if ord_T L_p <= r = rank E(Q), then ord_T f_E = r, so Sha(E/Q)[p^oo] is finite and the
p-adic height pairing is non-degenerate (Schneider's conjecture for E at p); and in general
rank E(Q) <= corank Sel_{p^oo}(E/Q) <= ord_T L_p(E,T).

Trusted software: PARI/GP (msfromell, mseval, ellap, ellglobalred). All p-adic arithmetic is
exact integer arithmetic modulo p^P in Python.
"""
from __future__ import annotations

import math
import time
from fractions import Fraction
from typing import Sequence

from .pari_session import get_pari

REFERENCES = {
    "Kato04": "K. Kato, p-adic Hodge theory and values of zeta functions of modular forms, Asterisque 295 "
              "(2004), 117-290, Theorem 17.4.",
    "Maz72": "B. Mazur, Rational points of abelian varieties with values in towers of number fields, "
             "Invent. Math. 18 (1972), 183-266 (control theorem).",
    "MTT86": "B. Mazur, J. Tate, J. Teitelbaum, On p-adic analogues of the conjectures of Birch and "
             "Swinnerton-Dyer, Invent. Math. 84 (1986), 1-48.",
    "SW13": "W. Stein, C. Wuthrich, Algorithms for the arithmetic of elliptic curves using Iwasawa theory, "
            "Math. Comp. 82 (2013), 1757-1792.",
}

_GP_VALUES = None


def _gp_values():
    """GP closure returning [m(a/q) : a = a0..a1] (exact rationals), m(r) = x({oo -> r})."""
    global _GP_VALUES
    if _GP_VALUES is None:
        _GP_VALUES = get_pari()("(M,xp,q,a0,a1)->vector(a1-a0+1,i,mseval(M,xp,[oo,(a0+i-1)/q]))")
    return _GP_VALUES


# ------------------------------------------------------------------------------------------
# P^1(Z/N) and SL_2(Z) lifts
# ------------------------------------------------------------------------------------------
def _factor(n: int) -> list[tuple[int, int]]:
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            out.append((d, e))
        d += 1 if d == 2 else 2
    if n > 1:
        out.append((n, 1))
    return out


def p1_list(N: int) -> list[tuple[int, int]]:
    """All points (c : d) of P^1(Z/NZ), one representative each, as pairs of residues mod N.

    P^1(Z/l^e) = {(c : 1) : c mod l^e} u {(1 : d) : d = 0 mod l}, and P^1(Z/N) is the product
    over l^e || N (Chinese remainder theorem; the scalars (Z/N)^x act componentwise).
    """
    if N == 1:
        return [(0, 1)]
    pts = [(0, 1)]
    mod = 1
    for l, e in _factor(N):
        q = l ** e
        local = [(c, 1) for c in range(q)] + [(1, d) for d in range(0, q, l)]
        new = []
        inv_mod = pow(mod, -1, q)
        inv_q = pow(q, -1, mod) if mod > 1 else 0
        for c0, d0 in pts:
            for c1, d1 in local:
                c = (c0 * q * inv_q + c1 * mod * inv_mod) % (mod * q)
                d = (d0 * q * inv_q + d1 * mod * inv_mod) % (mod * q)
                new.append((c, d))
        pts, mod = new, mod * q
    assert len(pts) == N * math.prod(1 + Fraction(1, l) for l, _ in _factor(N))
    return pts


def lift_to_sl2z(c: int, d: int, N: int) -> tuple[int, int, int, int]:
    """Integers (a, b, c', d') with a d' - b c' = 1 and (c', d') = (c, d) mod N; needs gcd(c, d, N) = 1."""
    c %= N
    d %= N
    if c == 0:
        c = N
    t = 0
    while math.gcd(c, d + t * N) != 1:
        t += 1
    d = d + t * N
    # extended gcd: x d + y c = 1
    r0, r1, x0, x1, y0, y1 = d, c, 1, 0, 0, 1
    while r1:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        x0, x1 = x1, x0 - qq * x1
        y0, y1 = y1, y0 - qq * y1
    assert r0 == 1 and x0 * d + y0 * c == 1
    return x0, -y0, c, d


# ------------------------------------------------------------------------------------------
# the primitive integral plus symbol
# ------------------------------------------------------------------------------------------
class PlusSymbol:
    """m'(r) = x({oo -> r}) / g, with x = PARI msfromell(E, 1) and g > 0 the generator of the
    Z-module spanned by the values of x on all Manin symbols (so m' is Z-valued and primitive)."""

    def __init__(self, E, N: int | None = None):
        pari = get_pari()
        t0 = time.time()
        self.E = E
        self.N = int(pari.ellglobalred(E)[0]) if N is None else int(N)
        self.M, self.xp = pari.msfromell(E, 1)
        self._m = {}
        vals = self._manin_values()
        nz = [v for v in vals if v != 0]
        if not nz:
            raise RuntimeError("plus symbol vanishes on every Manin symbol")
        num, den = 0, 1
        for v in nz:
            num = math.gcd(num, v.numerator)
            den = den * v.denominator // math.gcd(den, v.denominator)
        self.g = Fraction(num, den)
        self.n_manin = len(vals)
        self.seconds = round(time.time() - t0, 3)

    def raw(self, a: int, q: int) -> Fraction:
        """x({oo -> a/q}) as an exact Fraction (q >= 1)."""
        key = (a % q, q)
        v = self._m.get(key)
        if v is None:
            pari = get_pari()
            w = pari.mseval(self.M, self.xp, pari([pari("oo"), pari(a % q) / q]))
            v = Fraction(int(pari.numerator(w)), int(pari.denominator(w)))
            self._m[key] = v
        return v

    def _manin_values(self) -> list[Fraction]:
        out = []
        for c, d in p1_list(self.N):
            a, b, c1, d1 = lift_to_sl2z(c, d, self.N)
            # g{0 -> oo} = {b/d1 -> a/c1}; its value is m(a/c1) - m(b/d1), with m(oo) = 0
            v1 = self.raw(a, c1) if c1 > 0 else (self.raw(-a, -c1) if c1 < 0 else Fraction(0))
            if d1 == 0:
                v0 = Fraction(0)
            elif d1 > 0:
                v0 = self.raw(b, d1)
            else:
                v0 = self.raw(-b, -d1)
            out.append(v1 - v0)
        return out

    def __call__(self, a: int, q: int) -> int:
        v = self.raw(a, q) / self.g
        if v.denominator != 1:
            raise ArithmeticError(f"non-integral primitive symbol value at {a}/{q}: {v}")
        return v.numerator

    def values_mod(self, q: int) -> list[int]:
        """[m'(a/q) : 0 <= a < q] as Python ints (one vectorised PARI call)."""
        pari = get_pari()
        vec = _gp_values()(self.M, self.xp, q, 0, q - 1)
        out = []
        g_num, g_den = self.g.numerator, self.g.denominator
        for w in vec:
            num, den = int(pari.numerator(w)), int(pari.denominator(w))
            # w / g = num * g_den / (den * g_num) must be an integer
            top, bot = num * g_den, den * g_num
            if top % bot:
                raise ArithmeticError(f"non-integral primitive symbol value at denominator {q}")
            out.append(top // bot)
        return out

    def hecke_check(self, primes: Sequence[int] = (2, 3, 5, 7), cusps: Sequence[tuple[int, int]] =
                    ((1, 3), (2, 7), (5, 12), (3, 11))) -> bool:
        """sum_{j<l} m'((r+j)/l) + m'(l r) = a_l m'(r) (l not dividing N) at a few cusps r = a/q."""
        pari = get_pari()
        for l in primes:
            if self.N % l == 0:
                continue
            al = int(pari.ellap(self.E, l))
            for a, q in cusps:
                lhs = sum(self(a + j * q, l * q) for j in range(l))
                g = math.gcd(l * a, q)
                lhs += self(l * a // g, q // g)
                if lhs != al * self(a, q):
                    return False
        return True


# ------------------------------------------------------------------------------------------
# p-adic helpers (exact integers modulo p^P)
# ------------------------------------------------------------------------------------------
def unit_root(ap: int, p: int, P: int) -> int:
    """The unit root of X^2 - a_p X + p modulo p^P (Hensel), for p not dividing a_p."""
    if ap % p == 0:
        raise ValueError("p divides a_p: no unit root (supersingular)")
    pP = p ** P
    x = ap % p
    for _ in range(2 * P + 5):
        f = (x * x - ap * x + p) % pP
        fp = (2 * x - ap) % pP
        x = (x - f * pow(fp, -1, pP)) % pP
    assert (x * x - ap * x + p) % pP == 0 and x % p != 0
    return x


def vp(x: int, p: int) -> int | None:
    """p-adic valuation of an integer (None for 0)."""
    if x == 0:
        return None
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def coefficient_precision(p: int, n: int, j: int) -> int | None:
    """Guaranteed precision e of coeff_j(L^{(n)}) as an approximation of coeff_j(L), from (*).

    None means exact (j = 0)."""
    if j == 0:
        return None
    return n - 1 - int(math.floor(math.log(j, p) + 1e-12))


def riemann_coefficients(sym: PlusSymbol, p: int, ap: int, n: int, jmax: int, P: int | None = None) -> dict:
    """T^j coefficients (j <= jmax) of the level-n Riemann sum L^{(n)}(T), modulo p^P (default P = n + 1).

    Every coefficient c_j is returned as an integer in [0, p^P) with the precision e_j of (*): the
    true coefficient of L(T) is congruent to c_j modulo p^{e_j} (e_0 = exact).
    """
    if n < 1:
        raise ValueError("level n must be >= 1")
    if P is None:
        P = n + 1
    pP = p ** P
    alpha = unit_root(ap, p, P)
    ainv = pow(alpha, -1, pP)
    z_n = pow(ainv, n, pP)                 # alpha^{-n}
    z_n1 = z_n * ainv % pP                 # alpha^{-n-1}
    pn, pn1 = p ** n, p ** (n - 1)
    Vn = sym.values_mod(pn)
    Vn1 = sym.values_mod(pn1) if n > 1 else [sym(0, 1)]
    teich = [pow(t, pn1, pn) for t in range(1, p)]      # omega(t) mod p^n
    S = [0] * pn1
    gk = 1
    gamma = 1 + p
    for k in range(pn1):
        s = 0
        for w in teich:
            b = w * gk % pn
            s += z_n * Vn[b] - z_n1 * Vn1[b % pn1]
        S[k] = s % pP
        gk = gk * gamma % pn
    coeffs = []
    for j in range(jmax + 1):
        c = 0
        for k in range(j, pn1):
            if S[k]:
                c += S[k] * math.comb(k, j)
        coeffs.append(c % pP)
    prec = [coefficient_precision(p, n, j) for j in range(jmax + 1)]
    exact_constant = (1 - ainv) ** 2 * sym(0, 1) % pP
    if coeffs[0] != exact_constant:
        raise ArithmeticError("constant term disagrees with (1 - 1/alpha)^2 m'(0): bug or wrong symbol")
    return dict(p=p, n=n, P=P, coeffs=coeffs, precision=prec, alpha_mod=alpha, symbol_evaluations=pn + pn1)


def multiplicative_bound(E, p: int, *, jmax: int = 4, n: int = 2, sym: PlusSymbol | None = None) -> dict:
    """Certified upper bound for ord_{T=0} L_p(E,T) at an odd prime p of multiplicative reduction (p || N).

    Here f is already a U_p-eigenform with eigenvalue a_p = +1 (split) or -1 (non-split), and
    mu(a + p^k Z_p) = a_p^{-k} m'(a/p^k) is a Z_p-valued measure: finite additivity is
    U_p m = sum_{j<p} m((r+j)/p) = a_p m. The Riemann sums and the error bound (*) are as in the
    ordinary case, and the constant term is (1 - 1/a_p) m'(0) exactly. That constant vanishes in the
    split case, which is the exceptional zero of Mazur-Tate-Teitelbaum.

    Stein-Wuthrich Thm 7.3 (p semistable, rho_p surjective): L_p = f_E d(T), and T | d(T) when the
    reduction is split. Hence ord_T f_E <= ord_T L_p - 1 (split) or <= ord_T L_p (non-split). When E has
    multiplicative reduction at p >= 5, rhobar_p can fail to be surjective only if E has a rational
    p-isogeny (Serre, as quoted in SW section 7), and then rho_p is surjective by SW Prop 7.2.
    """
    pari = get_pari()
    N = int(pari.ellglobalred(E)[0])
    if p < 3 or N % p or (N // p) % p == 0:
        raise ValueError("needs an odd prime p with p || N")
    ap = int(pari.ellap(E, p))
    if ap not in (1, -1):
        raise ValueError("a_p must be +-1 at a multiplicative prime")
    if sym is None:
        sym = PlusSymbol(E, N)
    P = n + 1
    pP, q, q1 = p ** P, p ** n, p ** (n - 1)
    V = sym.values_mod(q)
    apinv_n = pow(ap, -n, pP)
    teich = [pow(t, q1, q) for t in range(1, p)]
    S, gk = [0] * q1, 1
    for k in range(q1):
        S[k] = sum(apinv_n * V[w * gk % q] for w in teich) % pP
        gk = gk * (1 + p) % q
    coeffs = [sum(S[k] * math.comb(k, j) for k in range(j, q1)) % pP for j in range(jmax + 1)]
    if coeffs[0] != (1 - pow(ap, -1, pP)) * sym(0, 1) % pP:
        raise ArithmeticError("constant term disagrees with (1 - 1/a_p) m'(0)")
    rc = dict(p=p, n=n, P=P, coeffs=coeffs, precision=[coefficient_precision(p, n, j) for j in range(jmax + 1)])
    cert = _read_certificate(rc)
    split = ap == 1
    r = cert["ord_upper_bound"]
    fE = None if r is None else (r - 1 if split else r)
    degrees = {int(d) for row in pari.Vec(pari.ellisomat(E)[1]) for d in row}
    no_p_isogeny = all(d % p for d in degrees)
    return dict(p=p, a_p=ap, reduction="split multiplicative" if split else "nonsplit multiplicative",
                level=n, modulus=f"{p}^{P}", coeffs=coeffs, facts=cert["facts"], ord_upper_bound=r,
                char_series_ord_upper_bound=fE, isogeny_degrees=sorted(degrees),
                rho_p_surjective_by_isogeny_criterion=bool(p >= 5 and no_p_isogeny and int(pari.elliscm(E)) == 0),
                note="ord_T f_E <= ord_T L_p - 1 (split, SW Thm 7.3) needs rho_p surjective; "
                     "for p >= 5 that holds when E has no rational p-isogeny")


def _read_certificate(rc: dict) -> dict:
    """From Riemann-sum coefficients, the certified facts about L(T)."""
    p = rc["p"]
    facts = []
    first_nonzero = None
    for j, (c, e) in enumerate(zip(rc["coeffs"], rc["precision"])):
        v = vp(c, p)
        if e is None:                                         # exact constant term (mod p^P)
            nonzero = c != 0
            facts.append(dict(j=j, value_mod_pP=c, valuation=v, precision="exact", certified_nonzero=nonzero,
                              known_zero_mod=None if nonzero else rc["P"]))
        else:
            nonzero = v is not None and v < e
            facts.append(dict(j=j, value_mod_pP=c, valuation=v, precision=e, certified_nonzero=nonzero,
                              known_zero_mod=None if nonzero else e))
        if nonzero and first_nonzero is None:
            first_nonzero = j
    return dict(facts=facts, ord_upper_bound=first_nonzero)


def ordinary_bound(E, p: int, *, jmax: int = 3, n_max: int | None = None, max_evaluations: int = 4 * 10 ** 6,
                   target: int | None = None, sym: PlusSymbol | None = None) -> dict:
    """Certified upper bound for ord_{T=0} L_p(E,T), p odd of good ordinary reduction.

    Raises the level n from 2 until some coefficient c_j, j <= jmax, is certified non-zero (or, with
    `target`, until c_target is), or until the next level would need more than `max_evaluations`
    modular-symbol values. Returns the certificate as a dict ('ord_upper_bound' is None if nothing
    was certified).
    """
    pari = get_pari()
    N = int(pari.ellglobalred(E)[0])
    if p < 3 or N % p == 0:
        raise ValueError("needs an odd prime of good reduction")
    ap = int(pari.ellap(E, p))
    if ap % p == 0:
        raise ValueError(f"p = {p} is supersingular (a_p = {ap})")
    t0 = time.time()
    if sym is None:
        sym = PlusSymbol(E, N)
    levels = []
    n = 2
    best = None
    while True:
        if p ** n + p ** (n - 1) > max_evaluations or (n_max is not None and n > n_max):
            break
        rc = riemann_coefficients(sym, p, ap, n, jmax)
        cert = _read_certificate(rc)
        levels.append(dict(n=n, P=rc["P"], coeffs=rc["coeffs"], ord_upper_bound=cert["ord_upper_bound"]))
        best = (rc, cert)
        r = cert["ord_upper_bound"]
        if r is not None and (target is None or r <= target):
            break
        n += 1
    out = dict(p=p, a_p=ap, conductor=N, reduction="good ordinary", symbol_content=str(sym.g),
               manin_symbols=sym.n_manin, levels=levels, seconds=round(time.time() - t0, 3))
    if best is None:
        out.update(ord_upper_bound=None, facts=[], note="budget too small for level 2")
        return out
    rc, cert = best
    out.update(level=rc["n"], modulus=f"{p}^{rc['P']}", ord_upper_bound=cert["ord_upper_bound"],
               facts=cert["facts"])
    r = cert["ord_upper_bound"]
    if r is not None:
        f = cert["facts"][r]
        out["leading"] = dict(j=r, valuation=f["valuation"], precision=f["precision"])
    return out


def selmer_bound(ainvs: Sequence[int], p: int, **kw) -> dict:
    """corank Sel_{p^oo}(E/Q) <= r with the theorem chain spelled out (good ordinary odd p)."""
    pari = get_pari()
    E = pari.ellinit([int(a) for a in ainvs])
    E = pari.ellinit(pari.ellminimalmodel(E)[:5])
    res = ordinary_bound(E, p, **kw)
    r = res["ord_upper_bound"]
    res["theorems"] = [
        "Kato 2004 Thm 17.4 as stated in Stein-Wuthrich 2013 Thm 7.3/7.4: f_E | p^m L_p(E,T) "
        "(p odd good ordinary, E without CM; for p = 3 needs rho_3 surjective or rhobar_3 not surjective)",
        "Mazur 1972 control theorem; Schneider/Perrin-Riou (SW Thm 6.1): rank <= corank Sel <= ord_T f_E",
        "exact Riemann sums of the integral plus symbol with the error bound (*) of bsd/kato_bound.py",
    ]
    cm = int(pari.elliscm(E)) != 0
    res["hypotheses_unverified"] = ([f"p = 3: rho_3 surjective or rhobar_3 not surjective"] if p == 3 else []) + \
                                   (["E has CM: outside Stein-Wuthrich's running hypothesis"] if cm else [])
    res["selmer_corank_upper_bound"] = r
    res["statement"] = (f"corank_Zp Sel_{{p^oo}}(E/Q) <= ord_T L_p(E,T) <= {r} (p = {p})" if r is not None
                        else f"no bound certified at p = {p} within the budget")
    return res


# ------------------------------------------------------------------------------------------
# independence of rational points, exactly (reduction modulo primes)
# ------------------------------------------------------------------------------------------
def _is_origin(P) -> bool:
    return len(P) == 1                                    # PARI represents the origin as [0]


def independence_certificate(E, points: Sequence[Sequence], torsion_order: int, *, ells: Sequence[int] =
                             (2, 3, 5, 7, 11, 13), max_q: int = 20000, max_columns: int = 80) -> dict:
    """Exact certificate that the points P_1..P_r of E(Q) are independent modulo torsion.

    Proof. Suppose sum n_i P_i is torsion with (n_i) != 0. Dividing by gcd(n_i) keeps the sum torsion,
    so take gcd(n_i) = 1. For a prime ell not dividing #E(Q)_tors the torsion point is ell times a
    torsion point, so for every prime q of good reduction (q != ell) the reductions satisfy
    sum (n_i mod ell) P_i = 0 in E(F_q)/ell E(F_q), a non-trivial F_ell-relation. When the
    ell-Sylow subgroup S of E(F_q) is cyclic of order ell^v, write #E(F_q) = ell^v m. Then
    P -> ell^{v-1} m P is a surjection E(F_q) -> S[ell] = Z/ell with kernel ell E(F_q), read off as
    k(P) with ell^{v-1} m P = k G, G a generator of S[ell]. If the r x (#q) matrix (k(P_i)) over
    F_ell has rank r, no such relation exists: the points are independent. Exact arithmetic only.
    """
    pari = get_pari()
    r = len(points)
    if r == 0:
        return dict(ok=True, ell=None, primes=[], method="empty set")
    N = int(pari.ellglobalred(E)[0])
    disc = int(pari("(E)->E.disc")(E))
    pts = [[Fraction(x), Fraction(y)] for x, y in points]
    ainvs = [int(a) for a in E[:5]]
    for ell in ells:
        if torsion_order % ell == 0:
            continue
        cols, used = [], []
        q = 2
        while q < max_q and len(cols) < max_columns:
            q = int(pari.nextprime(q + 1))
            if disc % q == 0 or q == ell:
                continue
            Eq = pari.ellinit(ainvs, q)
            nq = int(pari.ellcard(Eq))
            if nq % ell:
                continue
            cyc = [int(t) for t in pari.ellgroup(Eq)]
            if sum(1 for t in cyc if t % ell == 0) != 1:
                continue                                  # ell-Sylow not cyclic
            v, m = 0, nq
            while m % ell == 0:
                m //= ell
                v += 1
            c = ell ** (v - 1) * m
            G = None
            for _ in range(60):                           # c R generates S[ell] unless it is 0
                cand = pari.ellmul(Eq, pari.random(Eq), c)
                if not _is_origin(cand):
                    G = cand
                    break
            if G is None:
                continue
            mult = [pari("[0]")]
            for _ in range(1, ell):
                mult.append(pari.elladd(Eq, mult[-1], G))
            keys = [str(T) for T in mult]
            col = []
            for x, y in pts:
                if x.denominator % q == 0:
                    Pq = pari("[0]")                      # reduces to the origin
                else:
                    Pq = pari([pari.Mod(x.numerator * pow(x.denominator, -1, q), q),
                               pari.Mod(y.numerator * pow(y.denominator, -1, q), q)])
                    if not pari.ellisoncurve(Eq, Pq):
                        raise ArithmeticError(f"reduction of {x, y} mod {q} is not on the curve")
                col.append(keys.index(str(pari.ellmul(Eq, Pq, c))))
            cols.append(col)
            used.append(q)
            mat = pari.matrix(r, len(cols), [cols[j][i] for i in range(r) for j in range(len(cols))])
            if int(pari.matrank(pari.Mod(1, ell) * mat)) == r:
                return dict(ok=True, ell=ell, primes=used,
                            method="F_ell-rank of the images in prod_q E(F_q)/ell E(F_q) (exact)")
    return dict(ok=False, ell=None, primes=[], method="failed (no full-rank matrix found)")
