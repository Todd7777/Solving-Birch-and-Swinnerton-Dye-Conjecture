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


# ------------------------------------------------------------------------------------------
# supersingular primes with a_p = 0
# ------------------------------------------------------------------------------------------
SUPERSINGULAR_PROOF = r"""
Supersingular p with a_p = 0 (p odd, p not dividing N). Notation: Gamma_n = (1+pZ_p)/(1+p^{n+1}Z_p), gamma = 1+p,
for a unit b put l(b) in Z/p^n with <b> = gamma^{l(b)} mod p^{n+1}, omega_n(T) = (1+T)^{p^n} - 1, and let
    theta_n(T) = sum_{b in (Z/p^{n+1})^x} m'(b/p^{n+1}) (1+T)^{l(b)}       (degree < p^n, integer coefficients)
be the Mazur-Tate element of level p^{n+1}.

(1) Three-term relation. Grouping b by c = b mod p^n and using the Hecke relation
    sum_{j<p} m'((r+j)/p) + m'(pr) = a_p m'(r) = 0  at r = c/p^n
gives, for n >= 2, theta_n == -Phi_{p^{n-1}}(1+T) theta_{n-2}  (mod omega_{n-1}), and theta_1 == -(p-1) m'(0)
(mod omega_0 = T). (The p lifts of d mod p^{n-1} to c mod p^n have l(c) = l(d) + i p^{n-2}, i < p.)
Evaluating at zeta - 1, zeta of order p^t: theta_n(zeta-1) = 0 if t = n-1, and theta_n(zeta-1) = -p theta_{n-2}(zeta-1)
if t <= n-2. By induction theta_{2k+1}(zeta-1) = 0 for every even t with 2 <= t <= 2k. So theta_{2k+1} is
divisible in Z[T] by the monic, separable polynomial omega^+_{2k} = prod_{even t in [2, 2k]} Phi_{p^t}(1+T).

(2) The p-adic L-function. alpha = sqrt(-p) and beta = -alpha are the roots of X^2 + p. Stein-Wuthrich (Prop 3.1,
section 3.5) define L_alpha = lim P_n coefficientwise, P_n = sum mu_alpha(omega(a)(1+p)^j + p^n)(1+T)^j. (Their
convergence proof bounds the denominators c_k of P_k uniformly, Lemma 3.6; here c_k <= (k+1)/2 instead, and their
Prop 3.5 still gives ord_p(a_{k+1,j} - a_{k,j}) >= k/2 - 2 - floor(log_p j) -> oo, so the limit exists.) The
D_p-valued L_p = G omega_E - p H phi(omega_E) where L_alpha = G + alpha H, G, H in Q_p[[T]]. Since omega_E, phi(omega_E)
is a basis, ord_T L_p = min(ord_T G, ord_T H). From mu_alpha one gets P_n = alpha^{-n} theta_{n-1}
- alpha^{-n-1} Phi_{p^{n-1}}(1+T) theta_{n-2}, so for n = 2k+2 the Q_p-part is G_n = (-p)^{-(k+1)} theta_{2k+1}.
Write log^+_{<=2k} = (1/p) prod_{even t in [2, 2k]} Phi_{p^t}(1+T)/p, which tends coefficientwise to Pollack's
log^+ = (1/p) prod_{t even} Phi_{p^t}(1+T)/p (constant term 1/p). By (1),
    G_{2k+2} = log^+_{<=2k} * L^-_k,    L^-_k := (-1)^{k+1} theta_{2k+1} / omega^+_{2k}  in Z[T]
(because omega^+_{2k} = p^{k+1} log^+_{<=2k}). Hence L^- := G / log^+ = lim_k L^-_k exists in Z_p[[T]] (integral),
and ord_T G = ord_T L^-. Likewise, for n = 2k+1 the alpha-part is H_n = (-p)^{-(k+1)} theta_{2k}, theta_{2k} is
divisible by omega^-_{2k-1} = prod_{odd t <= 2k-1} Phi_{p^t}(1+T) = p^{k+1} log^-_{<=2k-1}, and
L^+_k := (-1)^{k+1} theta_{2k} / omega^-_{2k-1} tends to L^+ := H / log^- in Z_p[[T]].

(3) Congruence. G_n and G_{n'} agree at every root of omega_{n-1} (SW Cor 3.3, componentwise), and
log^+_{<=2k}(zeta) = log^+_{<=2k'}(zeta) != 0 when zeta = 1 or zeta has odd order exponent t <= 2k+1. So
L^-_{k'} - L^-_k vanishes at the roots of T omega^-_{2k+1}, omega^-_{2k+1} = prod_{odd t <= 2k+1} Phi_{p^t}(1+T). Letting
k' -> oo and dividing (Weierstrass division in Z_p[[T]] by a distinguished separable polynomial):
    L^- == L^-_k  (mod T omega^-_{2k+1} Z_p[[T]]).          For k = 0:  L^- == -theta_1  (mod omega_1 Z_p[[T]]).
The same argument for H gives L^+ == L^+_k (mod T omega^+_{2k} Z_p[[T]]); for k = 1, L^+ == theta_2 / Phi_p(1+T).
Each of the k+1 (resp. k) cyclotomic factors has all T^i coefficients, i < p-1, divisible by p, so for j <= p-1:
    [T^j] L^- == [T^j] L^-_k  (mod p^{k+1}),        [T^j] L^+ == [T^j] L^+_k  (mod p^k).
In particular [T^j] L^- == -[T^j] theta_1 (mod p). If a coefficient [T^j], j <= r, of L^- or L^+ is certified
non-zero this way, then ord_T L_p = min(ord_T L^-, ord_T L^+) <= r.

(4) Selmer bound. Stein-Wuthrich Thm 7.3/7.4 (Kato; for p >= 5 one applies, by Prop 7.2) hold for semistable p,
including good supersingular p, with L_p and f_E D_p-valued (their section 6.2), so ord_T f_E <= ord_T L_p.
Their Thm 6.2 (Perrin-Riou) gives ord_T f_E >= rank E(Q), with equality iff the D_p-valued regulator is nonzero
and Sha(E/Q)[p^oo] is finite. So ord_T L_p <= rank E(Q) proves both.
"""

_GP_THETA1 = None


def _gp_theta1():
    """GP loop for theta_1. The p-1 units with a given l(b) are omega(i) gamma^k, i < p, and omega(p-i) = -omega(i)
    mod p^2; the plus symbol is even and 1-periodic, m'(1 - a/p^2) = m'(-a/p^2) = m'(a/p^2). So the inner sum is
    twice the sum over i <= (p-1)/2. The full sum is recomputed at k = 0 and must agree (error otherwise)."""
    global _GP_THETA1
    if _GP_THETA1 is None:
        _GP_THETA1 = get_pari()("""(M,xp,p,gnum,gden,J)->
{
  my(q=p^2, gk=Mod(1,q), tch=vector(p-1,t,lift(Mod(t,q)^p)), s, c=vector(J+1), h=(p-1)/2, full);
  full = sum(i=1,p-1, mseval(M,xp,[oo,tch[i]/q]));
  for(k=0,p-1,
    s=0;
    for(i=1,h, s += mseval(M,xp,[oo,lift(tch[i]*gk)/q]));
    s = 2*s;
    if(k==0 && s!=full, error("theta_1: plus symbol not even at level p^2"));
    s = s*gden/gnum;
    for(j=0,J, c[j+1] += s*binomial(k,j));
    gk *= (1+p));
  c
}""")
    return _GP_THETA1


def theta1_coefficients(sym: PlusSymbol, p: int, jmax: int = 3) -> list[int]:
    """Exact integer T^j coefficients (j <= jmax) of theta_1 (level p^2, degree < p), one GP loop."""
    c = _gp_theta1()(sym.M, sym.xp, p, sym.g.numerator, sym.g.denominator, jmax)
    return [int(x) for x in c]


_GP_THETA_X = None


def theta_X(sym: PlusSymbol, p: int, n: int) -> list[int]:
    """theta_n in the variable X = 1+T: exact integer coefficients S_k of X^k, k < p^n (level p^{n+1}).

    One GP loop over the p^n (p-1) units mod p^{n+1}; no table of all symbol values is stored."""
    global _GP_THETA_X
    if _GP_THETA_X is None:
        _GP_THETA_X = get_pari()("""(M,xp,p,n,gnum,gden)->
{
  my(q=p^(n+1), gk=Mod(1,q), tch=vector(p-1,t,lift(Mod(t,q)^(p^n))), S=vector(p^n), s);
  for(k=0,p^n-1,
    s=0;
    for(i=1,p-1, s += mseval(M,xp,[oo,lift(tch[i]*gk)/q]));
    S[k+1] = s*gden/gnum;
    gk *= (1+p));
  S
}""")
    return [int(x) for x in _GP_THETA_X(sym.M, sym.xp, p, n, sym.g.numerator, sym.g.denominator)]


def _divide_X(c: list[int], d: list[int]) -> tuple[list[int], list[int]]:
    """Exact division of polynomials in X with integer coefficients by a monic d: (quotient, remainder)."""
    c = list(c)
    dd = len(d) - 1
    if len(c) <= dd:
        return [0], c
    Q = [0] * (len(c) - dd)
    for i in range(len(c) - 1, dd - 1, -1):
        coef = c[i]
        if coef:
            Q[i - dd] = coef
            for j in range(dd + 1):
                c[i - dd + j] -= coef * d[j]
    return Q, c[:dd]


def _cyclotomic_X(p: int, m: int) -> list[int]:
    """Phi_{p^m}(X) = sum_{i<p} X^{i p^{m-1}} as a coefficient list."""
    out = [0] * ((p - 1) * p ** (m - 1) + 1)
    for i in range(p):
        out[i * p ** (m - 1)] = 1
    return out


def _X_to_T(c: list[int], jmax: int) -> list[int]:
    return [sum(ck * math.comb(k, j) for k, ck in enumerate(c) if ck and k >= j) for j in range(jmax + 1)]


def supersingular_structure_check(sym: PlusSymbol, p: int, levels: int = 3) -> dict:
    """Exact checks of steps (1)-(3) of SUPERSINGULAR_PROOF at small p (cost p^{levels+1} symbol values):
    theta_2 divisible by Phi_p(1+T) and theta_3 by Phi_{p^2}(1+T), and L^-_1 == L^-_0 (mod omega_1)."""
    out = {}
    th = {n: theta_X(sym, p, n) for n in range(1, levels + 1)}
    if levels >= 2:
        Q2, R2 = _divide_X(th[2], _cyclotomic_X(p, 1))
        out["theta_2 divisible by Phi_p(1+T)"] = not any(R2)
    if levels >= 3:
        Q3, R3 = _divide_X(th[3], _cyclotomic_X(p, 2))
        out["theta_3 divisible by Phi_{p^2}(1+T)"] = not any(R3)
        # L^-_1 = (+1) theta_3 / omega^+_2 and L^-_0 = -theta_1 must agree modulo omega_1 = X^p - 1
        red = [0] * p
        for k, ck in enumerate(Q3):
            red[k % p] += ck
        for k, ck in enumerate(th[1]):
            red[k % p] += ck
        out["L^-_1 == L^-_0 (mod omega_1)"] = not any(red)
    return out


def _divide_cyclotomic_X(c: list[int], p: int, m: int) -> tuple[list[int], list[int]]:
    """Exact division by Phi_{p^m}(X) in linear time: multiply by X^{p^{m-1}} - 1, divide by X^{p^m} - 1."""
    s, P = p ** (m - 1), p ** m
    F = [0] * (len(c) + s)
    for k, ck in enumerate(c):
        if ck:
            F[k + s] += ck
            F[k] -= ck
    if len(F) <= P:
        return [0], F
    Q = [0] * (len(F) - P)
    for k in range(len(F) - 1, P - 1, -1):
        if F[k]:
            Q[k - P] = F[k]
            F[k - P] += F[k]
            F[k] = 0
    return Q, F[:P]


def _divide_cyclotomics_X(c: list[int], p: int, ts) -> tuple[list[int], bool]:
    """Divide successively by Phi_{p^t}(X), t in ts; returns (quotient, all remainders zero)."""
    ok = True
    for t in ts:
        c, R = _divide_cyclotomic_X(c, p, t)
        ok = ok and not any(R)
    return c, ok


def _cyclo_product_X(p: int, ts) -> list[int]:
    out = [1]
    for t in ts:
        d = _cyclotomic_X(p, t)
        new = [0] * (len(out) + len(d) - 1)
        for i, a in enumerate(out):
            if a:
                for j, b in enumerate(d):
                    if b:
                        new[i + j] += a * b
        out = new
    return out


def plus_minus_approximation(sym: PlusSymbol, p: int, sign: str, k: int, jmax: int = 3) -> dict:
    """[T^j] L^sign for j <= jmax, modulo p^prec, from one Mazur-Tate element (SUPERSINGULAR_PROOF, (2)-(3)).

    sign '-': L^-_k = (-1)^{k+1} theta_{2k+1} / omega^+_{2k}, precision p^{k+1}, cost p^{2k+2} symbol values.
    sign '+': L^+_k = (-1)^{k+1} theta_{2k} / omega^-_{2k-1} (k >= 1), precision p^k, cost p^{2k+1}.
    The exact divisions are checked to leave remainder 0 (step (1))."""
    if sign == "-":
        n, ts, prec = 2 * k + 1, [t for t in range(2, 2 * k + 1, 2)], k + 1
    elif sign == "+" and k >= 1:
        n, ts, prec = 2 * k, [t for t in range(1, 2 * k, 2)], k
    else:
        raise ValueError("sign must be '-' (k >= 0) or '+' (k >= 1)")
    if jmax > p - 1:
        raise ValueError("the congruence is proved for j <= p-1")
    th = theta_X(sym, p, n) if n >= 2 else None
    if n == 1:
        c = theta1_coefficients(sym, p, jmax)          # theta_1 needs no division (omega^+_0 = 1)
        coeffs = [(-1) ** (k + 1) * x for x in c]
    else:
        Q, ok = _divide_cyclotomics_X(th, p, ts)
        if not ok:
            raise ArithmeticError(f"theta_{n} is not divisible by prod Phi_(p^t), t in {ts}")
        coeffs = [(-1) ** (k + 1) * x for x in _X_to_T(Q, jmax)]
    mod = p ** prec
    return dict(sign=sign, k=k, level=n + 1, precision=prec, coeffs_mod=[x % mod for x in coeffs],
                valuations=[vp(x % mod, p) for x in coeffs], symbol_evaluations=(p - 1) * p ** n)


def supersingular_bound(E, p: int, *, jmax: int = 3, sym: PlusSymbol | None = None, cross_check_plus: bool = False,
                        target: int | None = None, max_evaluations: int = 2 * 10 ** 8) -> dict:
    """Certified upper bound for ord_{T=0} L_p(E,T) (D_p-valued) at a supersingular prime p >= 5 with a_p = 0.

    Uses [T^j] L^- == -[T^j] theta_1 (mod p), proved in SUPERSINGULAR_PROOF; cost p^2 symbol values. If no
    coefficient up to `target` (default jmax) is certified, tries L^+ from theta_2 (mod p) and then L^- from
    theta_3 (mod p^2), within max_evaluations. With cross_check_plus, the L^+ check is always run.
    """
    pari = get_pari()
    N = int(pari.ellglobalred(E)[0])
    if p < 5 or N % p == 0:
        raise ValueError("needs a prime p >= 5 of good reduction")
    ap = int(pari.ellap(E, p))
    if ap != 0:
        raise ValueError(f"a_p = {ap} != 0: not the supersingular case treated here")
    t0 = time.time()
    if sym is None:
        sym = PlusSymbol(E, N)
    if target is None:
        target = jmax
    c = theta1_coefficients(sym, p, jmax)
    if c[0] != -(p - 1) * sym(0, 1):
        raise ArithmeticError("theta_1(0) != -(p-1) m'(0): bug or wrong symbol")
    Lminus = [(-x) % p for x in c]                     # [T^j] L^- mod p, j <= jmax < p
    attempts = [dict(sign="-", k=0, level=2, precision=1, coeffs_mod=Lminus)]

    def first_nonzero(a):
        return next((j for j, v in enumerate(a["coeffs_mod"]) if v), None)

    r = first_nonzero(attempts[0])
    plan = [("+", 1), ("-", 1)]
    for sign, k in plan:
        if not (cross_check_plus and sign == "+") and r is not None and r <= target:
            break
        cost = (p - 1) * p ** (2 * k if sign == "+" else 2 * k + 1)
        if cost > max_evaluations:
            continue
        a = plus_minus_approximation(sym, p, sign, k, jmax)
        attempts.append(a)
        ra = first_nonzero(a)
        if ra is not None and (r is None or ra < r):
            r = ra
    out = dict(p=p, a_p=ap, reduction="good supersingular", theta1_coeffs=c, Lminus_mod_p=Lminus,
               attempts=attempts, ord_upper_bound=r, seconds=round(time.time() - t0, 2),
               method="[T^j] of L^-, L^+ from Mazur-Tate elements (SUPERSINGULAR_PROOF)")
    plus = [a for a in attempts if a["sign"] == "+"]
    if plus:
        out["Lplus_mod_p"] = plus[0]["coeffs_mod"]
    return out


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
    """corank Sel_{p^oo}(E/Q) <= r with the theorem chain spelled out (good ordinary odd p, or supersingular
    p >= 5 with a_p = 0, which is dispatched to supersingular_bound)."""
    pari = get_pari()
    E = pari.ellinit([int(a) for a in ainvs])
    E = pari.ellinit(pari.ellminimalmodel(E)[:5])
    N = int(pari.ellglobalred(E)[0])
    if p >= 5 and N % p and int(pari.ellap(E, p)) == 0:
        res = supersingular_bound(E, p, **{k: v for k, v in kw.items() if k in ("jmax", "sym", "target")})
        r = res["ord_upper_bound"]
        res["theorems"] = [
            "Stein-Wuthrich 2013 Thm 7.3/7.4 (Kato) for semistable p, D_p-valued in the supersingular case: "
            "ord_T f_E <= ord_T L_p (p >= 5, E without CM)",
            "Perrin-Riou (SW Thm 6.2): ord_T f_E >= rank E(Q), equality iff D_p-regulator != 0 and Sha[p^oo] finite",
            "ord_T L_p = min(ord_T L^-, ord_T L^+) and [T^j] L^- == -[T^j] theta_1 mod p (SUPERSINGULAR_PROOF)",
        ]
        cm = int(pari.elliscm(E)) != 0
        res["hypotheses_unverified"] = ["E has CM: outside Stein-Wuthrich's running hypothesis"] if cm else []
        res["selmer_corank_upper_bound"] = r
        res["statement"] = (f"corank_Zp Sel_{{p^oo}}(E/Q) <= ord_T L_p(E,T) <= {r} (p = {p}, supersingular)"
                            if r is not None else f"no bound certified at p = {p}")
        return res
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
