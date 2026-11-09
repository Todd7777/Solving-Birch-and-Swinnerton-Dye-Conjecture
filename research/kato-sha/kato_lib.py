"""Certified upper bounds for ord_{T=0} of cyclotomic p-adic L-functions of elliptic curves
over Q, computed by EXACT Riemann sums of modular symbols, and the resulting finiteness of
Sha(E/Q)[p^oo] for curves of rank >= 2 (method of Stein-Wuthrich, Math. Comp. 82 (2013)).

Everything p-adic below is exact integer arithmetic modulo p^M; the only external inputs are

  * eclib's plus modular symbol attached to the isogeny class (via Sage), which we rescale to the
    *primitive integer-valued* symbol  m'(r) = m(r)/g, where g > 0 generates the Z-module spanned
    by the values of m on all Manin symbols of Gamma_0(N) (all values m(r), r in Q, lie in that
    Z-module by Manin's continued-fraction trick, so m' is Z-valued and primitive);
  * the a_p of E (Sage / PARI, exact);
  * the Mordell-Weil generators of Cremona's tables, whose independence is re-certified exactly
    (reduction modulo primes, see ``independence_certificate``).

Ordinary p (p odd, p not dividing N, a_p not divisible by p).  Let alpha be the unit root of
x^2 - a_p x + p, mu(a + p^k Z_p) = alpha^{-k} m'(a/p^k) - alpha^{-k-1} m'(a/p^{k-1}) (a Z_p-valued
measure on Z_p^x because m' is Z-valued and T_p m' = a_p m'), and
        L(T) = int_{Z_p^x} (1+T)^{log_gamma <x>} dmu(x),     gamma = 1+p,
which is a nonzero rational multiple of the Mazur-Swinnerton-Dyer p-adic L-function (trivial
Teichmueller component).  The level-n Riemann sum
        L^{(n)}(T) = sum_{t=1}^{p-1} sum_{k=0}^{p^{n-1}-1} mu(w(t) gamma^k + p^n Z_p) (1+T)^k
(w = Teichmueller lift) satisfies, for 1 <= j < p^{n-1},
        v_p( coeff_j(L) - coeff_j(L^{(n)}) ) >= n - 1 - floor(log_p j)                     (*)
because for x in the ball around w(t)gamma^k one has log_gamma<x> = k + p^{n-1}u, u in Z_p, and
v_p binom(k + p^{n-1}u, j) - binom(k, j)) >= min_{1<=i<=j} (n-1-v_p(i)).  The constant term is
exact (finite additivity).  CERTIFICATE: v_p(coeff_r(L^{(n)})) < n-1-floor(log_p r)  =>  coeff_r(L)
is nonzero  =>  ord_{T=0} L <= r.

Supersingular p with a_p = 0 (p odd).  Mazur-Tate elements
        theta_n = sum_{a in (Z/p^{n+1})^x} m'(a/p^{n+1}) [<a> mod p^{n+1}]  in Z[Gamma_n] = Z[X]/(X^{p^n}-1),
X = 1+T.  Pollack's functions L^-, L^+ in Lambda (here: L^- is the coefficient of log^+, L^+ the
coefficient of alpha*log^-, as in Skinner's AWS-2018 notes, Sec. 4.1.4) satisfy
        L^- == -theta_1                 mod  omega_1 = (1+T)^p - 1
        L^+ ==  theta_2 / Phi_p(1+T)    mod  T*Phi_{p^2}(1+T)
        L^- ==  theta_3 / Phi_{p^2}(1+T) mod T*Phi_p(1+T)*Phi_{p^3}(1+T)
(the quotients are exact polynomial divisions, checked to have zero remainder), which give the
T^j coefficients of L^-, L^+, L^- modulo p, p, p^2 respectively for 1 <= j <= p-1.  See report.md
for the derivation (including integrality of L^+-) and for the theorems used afterwards.
"""
from __future__ import annotations

import math
import time

from sage.all import EllipticCurve, Integer, QQ, Rational, ZZ, gcd as sage_gcd
from sage.modular.modsym.p1list import P1List

DATA = '/home/user/Solving-Birch-and-Swinnerton-Dye-Conjecture/data/ecdata'


# ----------------------------------------------------------------------------------------------
# data loading
# ----------------------------------------------------------------------------------------------
def load_classes(rank_wanted=2, Nmax=9999, extra=('5077a1',)):
    """Curve #1 of every isogeny class of rank `rank_wanted` with N <= Nmax (Cremona allcurves),
    plus the labels in `extra`.  Returns list of dicts with label, N, ainvs, rank, torsion, gens."""
    files = ['allcurves.00000-09999']
    rows = {}
    for fn in files:
        for line in open(f'{DATA}/{fn}'):
            f = line.split()
            N, iso, num = int(f[0]), f[1], int(f[2])
            lab = f'{N}{iso}{num}'
            if N > Nmax:
                continue
            if (num == 1 and int(f[4]) == rank_wanted) or lab in extra:
                rows[lab] = dict(label=lab, N=N, iso=f'{N}{iso}', ainvs=eval(f[3]), rank=int(f[4]),
                                 torsion=int(f[5]))
    for line in open(f'{DATA}/allgens.00000-09999'):
        f = line.split()
        lab = f'{f[0]}{f[1]}{f[2]}'
        if lab in rows:
            r = int(f[4])
            tors_struct = f[5]
            gens = []
            for tok in f[6:6 + r]:
                X, Y, Z = [int(t) for t in tok.strip('[]').split(':')]
                gens.append((X, Y, Z))
            rows[lab]['gens'] = gens
            rows[lab]['torsion_structure'] = tors_struct
    out = sorted(rows.values(), key=lambda d: (d['N'], d['iso']))
    return out


def point_from_proj(E, XYZ):
    X, Y, Z = XYZ
    return E(QQ(X) / Z, QQ(Y) / Z)


# ----------------------------------------------------------------------------------------------
# exact certificate that the given points are Z-independent
# ----------------------------------------------------------------------------------------------
def independence_certificate(E, pts, tors_order, max_q=5000, max_cols=60):
    """Exact certificate that the points `pts` are Z-independent in E(Q).

    If n_1 P_1 + ... + n_r P_r is torsion for a primitive vector (n_i), then for a prime ell not
    dividing #E(Q)_tors the images of the P_i in prod_q E(F_q)/ell E(F_q) (q of good reduction)
    satisfy the nontrivial F_ell-relation (n_i mod ell).  For q with cyclic ell-Sylow subgroup S_q
    (order ell^v), E(F_q)/ell = S_q/ell S_q = Z/ell via  P -> k  with  ell^{v-1} m P = k G',
    m = #E(F_q)/ell^v, G' a generator of S_q[ell].  We collect such columns until the r x (#q)
    matrix over F_ell has rank r.  Pure exact arithmetic over finite fields."""
    from sage.all import GF, matrix, prime_range
    r = len(pts)
    N = int(E.conductor())
    for ell in [2, 3, 5, 7, 11, 13]:
        if tors_order % ell == 0:
            continue
        cols, used = [], []
        for q in prime_range(3, max_q):
            if N % q == 0 or q == ell:
                continue
            Eq = E.change_ring(GF(q))
            nq = int(Eq.cardinality())
            if nq % ell:
                continue
            v, m = 0, nq
            while m % ell == 0:
                m //= ell
                v += 1
            inv = [int(d) for d in Eq.abelian_group().invariants()]
            if sum(1 for d in inv if d % ell == 0) != 1:
                continue            # ell-Sylow not cyclic: skip this q
            c = ell**(v - 1) * m
            Gp = None
            for _ in range(50):
                R = Eq.random_point()
                if c * R != Eq(0):
                    Gp = c * R      # generator of S_q[ell]
                    break
            if Gp is None:
                continue
            mult = [Eq(0)]
            for k in range(1, ell):
                mult.append(mult[-1] + Gp)
            col = []
            for P in pts:
                x, y = P.xy()
                if x.denominator() % q == 0:
                    Pq = Eq(0)
                else:
                    Pq = Eq(x, y)
                Q = c * Pq
                col.append(mult.index(Q))
            cols.append(col)
            used.append(int(q))
            Mx = matrix(GF(ell), r, len(cols), lambda i, j: cols[j][i])
            if Mx.rank() == r:
                return dict(ok=True, ell=ell, primes=used,
                            method='reduction: F_ell-rank of images in prod_q E(F_q)/ell E(F_q)')
            if len(cols) >= max_cols:
                break
    return dict(ok=False, ell=None, primes=[], method='failed')


# ----------------------------------------------------------------------------------------------
# primitive modular symbol
# ----------------------------------------------------------------------------------------------
class PrimitiveModSym:
    """m'(r) = m(r)/g with m = eclib plus modular symbol, g = generator of the Z-span of the values
    of m on all Manin symbols (P^1(Z/NZ)).  Values are exact Python ints."""

    def __init__(self, E):
        t0 = time.time()
        self.E = E
        self.N = int(E.conductor())
        self.m = E.modular_symbol(sign=1, implementation='eclib')
        self._ms = self.m._modsym
        self._scal = QQ(self.m._scaling)
        self.t_build = time.time() - t0
        t0 = time.time()
        self.cache = {}
        self.g, self.n_manin = self._content()
        self.t_content = time.time() - t0

    def raw(self, a, q):
        """eclib value (times Sage's scaling) at the cusp a/q, a, q ints, q >= 1."""
        key = (a % q, q)
        v = self.cache.get(key)
        if v is None:
            r = Rational((a % q, q))
            v = self._ms(r, sign=1, base_at_infinity=True) * self._scal
            self.cache[key] = v
        return v

    def _content(self):
        L = P1List(self.N)
        vals = []
        for i in range(len(L)):
            a, b, c, d = L.lift_to_sl2z(i)
            # psi = m(g.oo) - m(g.0) = m(a/c) - m(b/d); m(oo) = 0
            v1 = self.raw(a, c) if c != 0 else QQ(0)
            v0 = self.raw(b, d) if d != 0 else QQ(0)
            vals.append(QQ(v1 - v0))
        nz = [v for v in vals if v != 0]
        if not nz:
            raise RuntimeError('modular symbol vanishes on all Manin symbols')
        num = 0
        den = 1
        for v in nz:
            num = math.gcd(num, int(v.numerator()))
            den = den * int(v.denominator()) // math.gcd(den, int(v.denominator()))
        g = QQ(num) / den
        self._manin_vals = vals
        return g, len(L)

    def __call__(self, a, q):
        v = self.raw(a, q) / self.g
        if v.denominator() != 1:
            raise ArithmeticError(f'non-integral primitive symbol value at {a}/{q}: {v}')
        return int(v)

    def sanity_checks(self, primes=(2, 3, 5, 7, 11, 13), trials=4, seed=1):
        """Checks: m'(oo)=0 implicit; evenness m'(-r)=m'(r); Hecke relation
        sum_{j<l} m'((r+j)/l) + m'(l r) = a_l m'(r) for l not dividing N; Gamma_0(N)-invariance."""
        import random
        rnd = random.Random(seed)
        E = self.E
        out = dict(even=True, hecke=True, gamma0=True)
        for _ in range(trials):
            q = rnd.randint(2, 400)
            a = rnd.randint(1, q - 1)
            if math.gcd(a, q) != 1:
                continue
            if self(a, q) != self(-a, q):
                out['even'] = False
            for l in primes:
                if self.N % l == 0:
                    continue
                al = int(E.ap(l))
                # r = a/q ; (r+j)/l = (a + j q)/(l q) ; l r = l a / q
                lhs = sum(self(a + j * q, l * q) for j in range(l)) + self(l * a, q)
                if lhs != al * self(a, q):
                    out['hecke'] = False
            # Gamma_0(N) invariance: gamma = [[x, y],[N z, w]] det 1: m(gamma r) - m(gamma oo) = m(r) - m(oo)
            for _ in range(2):
                z = rnd.randint(1, 5)
                w = rnd.randint(1, 50)
                while math.gcd(w, self.N * z) != 1:
                    w += 1
                # find x,y with x w - y N z = 1
                g_, x, y = _xgcd(w, self.N * z)
                y = -y
                # gamma r = (x a + y q)/(N z a + w q), gamma oo = x/(N z)
                num, den = x * a + y * q, self.N * z * a + w * q
                if den < 0:
                    num, den = -num, -den
                lhs = self(num, den) - self(x, self.N * z)
                if lhs != self(a, q):
                    out['gamma0'] = False
        return out


def _xgcd(a, b):
    """returns g, x, y with a x + b y = g"""
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        qq, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - qq * x1
        y0, y1 = y1, y0 - qq * y1
    return a, x0, y0


# ----------------------------------------------------------------------------------------------
# p-adic helpers (exact integers mod p^M)
# ----------------------------------------------------------------------------------------------
def vp_int(x, p, cap=10**9):
    if x == 0:
        return cap
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def unit_root(ap, p, M):
    """unit root of x^2 - ap x + p modulo p^M (Newton/Hensel), requires p not dividing ap."""
    pM = p**M
    x = ap % p
    for _ in range(2 * M + 5):
        f = (x * x - ap * x + p) % pM
        fp = (2 * x - ap) % pM
        x = (x - f * pow(fp, -1, pM)) % pM
    assert (x * x - ap * x + p) % pM == 0 and x % p != 0
    return x


def teichmuller_reps(p, M):
    pM = p**M
    return [pow(t, p**(M - 1), pM) for t in range(1, p)]


def binom_coeffs(S, jmax, mod):
    """given S[k] (coefficient of (1+T)^k), return T^j coefficients for j<=jmax, mod `mod`."""
    out = []
    for j in range(jmax + 1):
        s = 0
        for k, sk in enumerate(S):
            if sk and k >= j:
                s += sk * math.comb(k, j)
        out.append(s % mod)
    return out


def ordinary_riemann(ms, p, ap, n, jmax, M):
    """T^j coefficients (j <= jmax) of the level-n Riemann sum L^{(n)}, mod p^M."""
    pM = p**M
    alpha = unit_root(ap, p, M)
    ainv = pow(alpha, -1, pM)
    z = pow(ainv, n, pM)
    z1 = z * ainv % pM
    pn = p**n
    pn1 = p**(n - 1)
    teich = teichmuller_reps(p, max(M, n))
    S = [0] * pn1
    gk = 1
    for k in range(pn1):
        s = 0
        for w in teich:
            b = (w * gk) % pn
            s += z * ms(b, pn) - z1 * ms(b % pn1, pn1)
        S[k] = s % pM
        gk = gk * (1 + p) % pn
    return binom_coeffs(S, jmax, pM)


def ordinary_precision(p, n, j):
    if j == 0:
        return 10**9   # exact
    return n - 1 - int(math.floor(math.log(j, p) + 1e-12))


def theta_X_coeffs(ms, p, n):
    """theta_n as exact integer coefficients c_k of X^k, X = 1+T, k < p^n."""
    q = p**(n + 1)
    teich = teichmuller_reps(p, n + 1)
    c = [0] * (p**n)
    gk = 1
    for k in range(p**n):
        s = 0
        for w in teich:
            s += ms((w * gk) % q, q)
        c[k] = s
        gk = gk * (1 + p) % q
    return c


def divide_by_cyclotomic(cX, p, m):
    """exact division of sum c_k X^k by Phi_{p^m}(X) = (X^{p^m}-1)/(X^{p^{m-1}}-1).
    Returns quotient coefficients (in X) and the remainder (must be zero)."""
    # multiply by (X^{p^{m-1}} - 1), then divide by (X^{p^m} - 1)
    s = p**(m - 1)
    P = p**m
    F = [0] * (len(cX) + s)
    for k, c in enumerate(cX):
        if c:
            F[k + s] += c
            F[k] -= c
    # long division by X^P - 1 from the top
    deg = len(F) - 1
    Q = [0] * max(deg - P + 1, 0)
    for k in range(deg, P - 1, -1):
        if F[k]:
            Q[k - P] = F[k]
            F[k - P] += F[k]
            F[k] = 0
    rem = F[:P]
    return Q, rem


def X_to_T(cX, jmax):
    """T^j coefficients (exact ints) of sum c_k (1+T)^k for j <= jmax."""
    return [sum(c * math.comb(k, j) for k, c in enumerate(cX) if c and k >= j) for j in range(jmax + 1)]
