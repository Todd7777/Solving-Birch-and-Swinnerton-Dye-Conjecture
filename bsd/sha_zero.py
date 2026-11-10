"""Certificates that the p-primary part Sha(E/Q)(p) is trivial, prime by prime.

This follows Stein-Wuthrich, Math. Comp. 82 (2013), Algorithm 11.1 (read for this module), whose output is an
upper bound #Sha(E/Q)(p) <= p^{b_p}. Since #Sha(E/Q)(p) is a square (Cassels), b_p <= 1 proves Sha(E/Q)(p) = 0,
equivalently Sha(E/Q)[p] = 0. Two steps of their algorithm are replaced by short proofs, given below, so that
no p-saturation of the Mordell-Weil basis and no Dickson classification is needed. The p-adic L-function inputs
come from bsd/kato_bound.py.

Theorem chain (Stein-Wuthrich; p odd, E without CM, not additive at p)
------------------------------------------------------------------------
* SW Thm 7.3 (Kato): if E is semistable at p and rho_p is surjective, then L_p(E,T) = f_E(T) d(T) with d in
  Lambda (and T | d in the split multiplicative case; L_p and f_E D_p-valued in the supersingular case).
  SW Prop 7.2: for p >= 5, rho_p is surjective iff rhobar_p is.
* SW Thm 6.1 (Schneider, Perrin-Riou, Jones): if ord_T f_E = r = rank E(Q), the leading coefficient f_E^* has the
  valuation of eps_p * prod_v c_v * #Sha(E/Q)(p) / #E(Q)(p)^2 * Reg_gamma(E/Q); eps_p = (1 - 1/alpha)^2 with
  ord_p(eps_p) = 2 ord_p(N_p), N_p = #E~(F_p) (SW section 3.1); in the split multiplicative case eps_p is replaced
  by L_p / log_p(kappa(gamma)), L_p = log_p(q_E)/ord_p(q_E) the L-invariant. Reg_gamma = Reg_p / log_p(1+p)^r.
* SW Thm 6.2 (Perrin-Riou), supersingular p: (1 - phi)^{-2} f_E^* has the valuation of prod_v c_v * #Sha(E/Q)(p)
  * Reg(E/Q), a D_p-valued statement. SW write Reg_p there and in Algorithm 11.1; with T-coefficients the
  consistent normalisation is Reg_gamma = Reg_p / log_p(1+p)^r (as in their example 12.2). This module uses
  Reg_gamma, which gives the larger, hence weaker, bound: with Reg_p the bound would be smaller by 2r.
  Likewise at split multiplicative p it uses L_p / log_p(kappa(gamma)) (Thm 6.1) rather than L_p (as printed in
  Algorithm 11.1), again the weaker choice.
Combining: ord_p #Sha(E/Q)(p) <= b_p := ord_p L_p^* - ord_p(eps_p) - sum ord_p c_v + 2 ord_p #E(Q)_tors
- ord_p Reg_gamma, because ord_p f_E^* <= ord_p L_p^* (Kato).

Lemma A (transvections; replaces Serre's/Dickson's criterion). Let G <= GL_2(F_l) contain a transvection.
Then G fixes a line (lies in a Borel subgroup) or G contains SL_2(F_l). Proof: if two transvections in G have
different axes, in a suitable basis they are (1 a; 0 1) and (1 0; b 1) with a, b != 0; their powers give all
upper and all lower unipotent matrices, which generate SL_2(F_l) (row reduction). Otherwise all transvections in
G have one axis L; conjugating a transvection by g in G gives a transvection with axis g(L), so g(L) = L.
Application: let E have multiplicative reduction at q with l not dividing v_q(q_E) = -v_q(j) = ord_q(Delta_min).
By Tate's uniformisation, over Q_q (after an unramified quadratic twist, invisible to inertia) E[l] is an
extension of Z/l by mu_l, and an inertia element sigma moves q_E^{1/l} to zeta^c q_E^{1/l}. Over
K = Q_q^ur(zeta_l), whose valuation is v_q (q != l) or (l-1) v_q (q = l), q_E has valuation prime to l, so it is not
an l-th power and K(q_E^{1/l})/K has degree l: some sigma in Gal(Qbar_q/K), inside inertia, acts as (1 c; 0 1),
c != 0, a transvection. (For q != l this is the ramification criterion quoted in SW section 7 from Ribet-Stein
[RS01, 2.4]; the case q = l is used for 389a1 at l = 389, the only bad prime.) A Frobenius at a good prime q' with
a_{q'}^2 - 4q' a non-square mod l has no eigenvalue in F_l, so the image is not in a Borel subgroup. Hence the
image contains SL_2(F_l); its determinant is the cyclotomic character, so it is GL_2(F_l).

Saturation (needed where PARI's regulator is used for the bound: supersingular and split multiplicative p).
If G = <P_1..P_r> has finite index in E(Q), E(Q)[p] = 0, and p divides the index, there is Q with pQ = sum a_i P_i,
a not 0 mod p. Reducing at a good prime l, sum a_i P_i~ lies in pE(F_l). When ord_p #E(F_l) = 1, P -> (#E(F_l)/p) P
maps E(F_l) onto E(F_l)[p] = Z/p with kernel exactly pE(F_l), so l gives a linear form f_l on F_p^r with f_l(a) = 0.
If the forms found have rank r, then a = 0: contradiction. Then every regulator computed from the P_i (real,
p-adic, D_p-valued) equals that of E(Q) times [E(Q):G]^2, a p-adic unit, so its valuation is that of E(Q).

Lemma B (integrality of p-adic heights; replaces the regulator computation and the saturation check). Let p be
odd of good ordinary reduction, with p not dividing prod_v c_v, and let P be in E(Q). SW (4.1) define
h(P) = h(mP)/m^2 where mP has good reduction everywhere and lies in the formal group, and there
h(mP) = 2 log_p(e/sigma_p(t)), t = -x/y, e^2 = denominator of x(mP), sigma_p in t + t^2 Z_p[[t]] (SW, citing
[MT91]). On a minimal model v_p(t) = v_p(e), so e/sigma_p(t) is a p-adic unit and h(mP) is in pZ_p. One can take
m prime to p unless p | N_p.
 (i) If p does not divide N_p, every height pairing value is in pZ_p, so for any basis of E(Q)/tors,
     ord_p Reg_p >= r and ord_p Reg_gamma >= 0.
 (ii) If ord_p N_p = 1, the reduction map E(Q) -> E~(F_p)[p] = Z/p has a kernel K of index 1 or p. Choosing a
     basis with P_2..P_r in K gives <P_1,P_1> in p^-1 Z_p, <P_1,P_j> = <pP_1,P_j>/p in Z_p, <P_i,P_j> in pZ_p, so
     ord_p Reg_p >= r - 2 and ord_p Reg_gamma >= -2, while ord_p eps_p = 2.
 In both cases b_p <= ord_p L_p^* - sum ord_p c_v + 2 ord_p #E(Q)_tors, for ANY basis of E(Q): no saturation
 proof is needed. The PARI regulator (ellpadicregulator, which reproduces SW (4.2) for 446d1 at p = 5 digit for
 digit) is computed only as a cross-check. If b_p = 0 the bound is attained, and by the remark after SW
 Algorithm 11.1 the Iwasawa main conjecture then holds for (E, p): f_E / L_p is a unit of Lambda.
"""
from __future__ import annotations

import time
from typing import Sequence

from .pari_session import get_pari


def _v(x, p) -> int | None:
    pari = get_pari()
    return None if x == 0 else int(pari.valuation(x, p))


# ------------------------------------------------------------------------------------------
# surjectivity of rhobar_l (Lemma A)
# ------------------------------------------------------------------------------------------
def rhobar_surjective_certificate(E, ell: int, qmax: int = 20000) -> dict:
    """Certificate that rhobar_ell is surjective via Lemma A: a multiplicative prime q with ell not dividing
    ord_q(Delta_min) = -v_q(j) (transvection in inertia at q; q != ell preferred, q = ell allowed), and a good
    prime q' with a_{q'}^2 - 4q' a non-square mod ell."""
    pari = get_pari()
    if ell < 3:
        raise ValueError("ell must be an odd prime")
    j = pari("(E)->E.j")(E)
    N = int(pari.ellglobalred(E)[0])
    trans = None
    fa = pari.factor(N)
    mult = [q for q, e in zip([int(x) for x in fa[0]], [int(x) for x in fa[1]]) if e == 1]
    for q in sorted(mult, key=lambda q: q == ell):              # q != ell first
        vq = -int(pari.valuation(j, q))
        if vq % ell != 0:
            trans = dict(q=q, ord_q_Delta_min=vq, reduction="multiplicative", q_equals_ell=q == ell)
            break
    if trans is None:
        return dict(ok=False, reason="no multiplicative prime q with ell not dividing ord_q(Delta_min)")
    q2 = 2
    while q2 < qmax:
        q2 = int(pari.nextprime(q2 + 1))
        if N % q2 == 0 or q2 == ell:
            continue
        aq = int(pari.ellap(E, q2))
        d = (aq * aq - 4 * q2) % ell
        if d and int(pari.kronecker(d, ell)) == -1:
            return dict(ok=True, ell=ell, transvection=trans, frobenius=dict(q=q2, a_q=aq, disc_mod_ell=d),
                        method="Lemma A: transvection at a multiplicative prime + Frobenius with irreducible char poly")
    return dict(ok=False, reason=f"no Frobenius with non-square discriminant mod {ell} below {qmax}")


# ------------------------------------------------------------------------------------------
# p-saturation of a set of generators
# ------------------------------------------------------------------------------------------
_CARD_CACHE: dict = {}


def _card(E, ell: int) -> int:
    key = (tuple(int(a) for a in E[:5]), ell)
    if key not in _CARD_CACHE:
        _CARD_CACHE[key] = ell + 1 - int(get_pari().ellap(E, ell))
    return _CARD_CACHE[key]


def _rank_mod_p(rows: list[list[int]], p: int) -> int:
    m = [list(r) for r in rows]
    rank, ncols = 0, len(m[0]) if m else 0
    for c in range(ncols):
        piv = next((i for i in range(rank, len(m)) if m[i][c] % p), None)
        if piv is None:
            continue
        m[rank], m[piv] = m[piv], m[rank]
        inv = pow(m[rank][c], -1, p)
        m[rank] = [x * inv % p for x in m[rank]]
        for i in range(len(m)):
            if i != rank and m[i][c] % p:
                f = m[i][c]
                m[i] = [(x - f * y) % p for x, y in zip(m[i], m[rank])]
        rank += 1
    return rank


def saturation_certificate(E, gens: Sequence, p: int, lmax: int = 10 ** 7, max_forms: int = 30) -> dict:
    """Proof that <gens> has index prime to p in E(Q), by reduction maps (see the module docstring); the gens
    must be independent points generating a subgroup of finite index, which the rank-r conclusion also shows
    up to torsion (a relation sum a_i P_i = 0 with a primitive would make every form vanish on a).

    If max_forms useful primes l leave the forms of rank < r, the common kernel points at the p-divisible
    combinations: each kernel vector a is tested with ellisdivisible(E, sum a_i P_i, p), and a hit is returned as
    a witness that p divides the index."""
    pari = get_pari()
    if int(pari.elltors(E)[0]) % p == 0:
        return dict(ok=False, reason="p divides #E(Q)_tors")
    N = int(pari.ellglobalred(E)[0])
    ainvs = [int(a) for a in E[:5]]
    pts = [[pari(c) for c in P] for P in gens]
    den = 1
    for P in pts:
        for c in P:
            den = den * int(pari.denominator(c))
    r = len(pts)
    forms, used = [], []
    ell = 2
    while ell < lmax:
        ell = int(pari.nextprime(ell + 1))
        if ell == p or N % ell == 0 or den % ell == 0:
            continue
        n = _card(E, ell)
        if n % p or n % (p * p) == 0:
            continue
        El = pari.ellinit(ainvs, ell)
        R = [pari.ellmul(El, [pari.Mod(c, ell) for c in P], n // p) for P in pts]
        nz = [i for i, x in enumerate(R) if len(x) == 2]
        if not nz:
            continue
        g = R[nz[0]]
        f = [int(pari.elllog(El, x, g, p)) if len(x) == 2 else 0 for x in R]
        forms.append(f)
        used.append(dict(l=ell, card=n, form=f))
        if _rank_mod_p(forms, p) == r:
            return dict(ok=True, p=p, primes=used, method="reduction maps E(Q)/p -> E(F_l)/p with ord_p #E(F_l) = 1")
        if len(forms) >= max_forms:
            break
    out = dict(ok=False, p=p, primes=used[:5], n_forms=len(forms), reason=f"forms of rank < {r}")
    if forms:
        K = pari.matker(pari.matrix(len(forms), r, [x for f in forms for x in f]) * pari.Mod(1, p))
        witnesses = []
        for k in range(int(pari.matsize(K)[1])):
            a = [int(pari.centerlift(K[i, k])) for i in range(r)]
            Q = [0]
            for ai, P in zip(a, pts):
                Q = pari.elladd(E, Q, pari.ellmul(E, P, ai))
            if int(pari.ellisdivisible(E, Q, p)):
                witnesses.append(a)
        if witnesses:
            out.update(not_saturated=True, witness=witnesses,
                       reason=f"sum a_i P_i is divisible by {p} in E(Q) for a = {witnesses[0]}: p divides the index")
    return out


# ------------------------------------------------------------------------------------------
# the bound b_p
# ------------------------------------------------------------------------------------------
def _local_data(E, p: int) -> dict:
    pari = get_pari()
    N = int(pari.ellglobalred(E)[0])
    cs = []
    for q in [int(x) for x in pari.factor(N)[0]]:
        cs.append(dict(q=q, c=int(pari.elllocalred(E, q)[3])))
    tors = int(pari.elltors(E)[0])
    return dict(N=N, tamagawa=cs, sum_vp_c=sum(int(pari.valuation(c["c"], p)) for c in cs),
                torsion=tors, vp_tors=int(pari.valuation(tors, p)))


def ordinary_bound(E, p: int, vp_Lstar: int, rank: int, gens: Sequence | None = None, regulator_digits: int = 6) -> dict:
    """b_p at a good ordinary odd prime p via Lemma B; with `gens`, also the exact b_p from PARI's regulator."""
    pari = get_pari()
    loc = _local_data(E, p)
    if loc["N"] % p == 0:
        raise ValueError("p must be a prime of good reduction")
    ap = int(pari.ellap(E, p))
    if ap % p == 0:
        raise ValueError("p is supersingular")
    Np = p + 1 - ap
    vN = int(pari.valuation(Np, p))
    if vN > 1:
        return dict(ok=False, reason=f"ord_p(N_p) = {vN} > 1: Lemma B not applicable")
    if loc["sum_vp_c"]:
        return dict(ok=False, reason="p divides a Tamagawa number: Lemma B not applicable")
    v_eps = 2 * vN
    reg_gamma_lb = 0 if vN == 0 else -2
    b = vp_Lstar - v_eps - loc["sum_vp_c"] + 2 * loc["vp_tors"] - reg_gamma_lb
    out = dict(p=p, a_p=ap, N_p=Np, vp_N_p=vN, vp_eps=v_eps, vp_Lstar=vp_Lstar, reg_gamma_lower_bound=reg_gamma_lb,
               tamagawa=loc["tamagawa"], torsion=loc["torsion"], b_p=b, sha_p_trivial=b <= 1)
    if gens is not None:
        R = pari.ellpadicregulator(E, p, regulator_digits, pari(gens))
        vR = _v(R, p)
        out["pari_vp_Reg_p"] = vR
        if vR is not None:
            exact = vp_Lstar - v_eps - loc["sum_vp_c"] + 2 * loc["vp_tors"] - (vR - rank)
            out["b_p_from_pari_regulator"] = exact          # valid if gens is a basis; >= 0 if consistent
            out["consistent"] = exact >= 0 and (vR - rank) >= reg_gamma_lb
    return out


def supersingular_bound(E, p: int, Lminus_T2_mod_p: int, rank: int, gens: Sequence, digits: int = 6,
                        Lplus_T2_mod_p: int | None = None) -> dict:
    """b_p at a supersingular p >= 5 (a_p = 0), Reg_gamma convention (the weaker one, see the module docstring).
    The bound is valid when `gens` is p-saturated (saturation_certificate).

    L^- = G / log^+ and L^+ = H / log^- are named as in bsd.kato_bound.SUPERSINGULAR_PROOF (there L^- is the part
    congruent to -theta_1; log^+-(0) = 1/p). The leading coefficient of L_p = G omega - p H phi(omega) is
    G_2 omega - p H_2 phi(omega) with G_2 = [T^2]L^-/p and H_2 = [T^2]L^+/p (ord_T L^+- >= 2). With phi^2 = -1/p,
    (1-phi)^{-2} L_p^* has omega-coordinate p((p-1)G_2 + 2pH_2)/(p+1)^2 and phi(omega)-coordinate
    p^2(2G_2 - (p-1)H_2)/(p+1)^2 (at T = 0 these give SW's interpolation (1-phi)^{-2} L_p(0) = (L(E,1)/Omega) omega,
    since L^-(0) = (p-1) L(E,1)/Omega and L^+(0) = 2 L(E,1)/Omega). If [T^2]L^- is a unit, the first coordinate has
    valuation exactly 0 and the second >= 1 (L^+ is integral), so the D_p-valuation is 0. PARI's D_p-valued
    regulator (basis omega, F omega, F = p phi) gives the valuation of Reg_p in the lattice
    Z_p omega + Z_p phi(omega); that lattice equals Z_p omega + Z_p eta when the eta-coordinate of phi(omega) is a
    unit, which is checked with PARI's Frobenius matrix.

    Direction check (when [T^2]L^+ mod p is known): Kato's L_p = f_E d with d in Lambda, and Perrin-Riou's formula
    for (1-phi)^{-2} f_E^* as a multiple of Reg_p, make (1-phi)^{-2} L_p^* a scalar multiple of Reg_p. So the ratio
    (phi(omega)-coordinate) / (p * omega-coordinate), which is -2 - [T^2]L^+/[T^2]L^- mod p, must equal the same
    ratio for PARI's regulator, R[2]/R[1] mod p. This is not used in the bound.
    """
    pari = get_pari()
    loc = _local_data(E, p)
    ap = int(pari.ellap(E, p))
    if ap != 0 or p < 5:
        raise ValueError("needs p >= 5 with a_p = 0")
    if Lminus_T2_mod_p % p == 0:
        return dict(ok=False, reason="[T^2]L^- not a unit: valuation of the leading term not determined here")
    vL = 0
    R = pari.ellpadicregulator(E, p, digits, pari(gens))
    r1, r2 = R[0], R[1]
    v_reg = min(x for x in (_v(r1, p), None if r2 == 0 else _v(r2, p) + 1) if x is not None)
    Fr = pari.ellpadicfrobenius(E, p, digits)
    # PARI's matrix of F on (omega, eta): column 1 = F(omega); phi(omega) = F(omega)/p
    eta_coord = Fr[1, 0] / p
    lattices_agree = _v(eta_coord, p) == 0 and _v(Fr[0, 0] / p, p) is not None and _v(Fr[0, 0] / p, p) >= 0
    v_reg_gamma = v_reg - rank
    b = vL - v_reg_gamma - loc["sum_vp_c"] + 2 * loc["vp_tors"]
    out = dict(p=p, a_p=ap, vp_leading=vL, pari_regulator=[str(r1), str(r2)], vp_Reg_p=v_reg,
               vp_Reg_gamma=v_reg_gamma, lattices_omega_eta_equal=lattices_agree, b_p=b,
               b_p_with_Reg_p_as_printed_in_SW=b - rank, sha_p_trivial=b <= 1 and lattices_agree,
               consistent=b >= 0)
    if Lplus_T2_mod_p is not None and _v(r1, p) == v_reg and r2 != 0 and _v(r2, p) >= v_reg:
        scale = pari(p) ** v_reg
        reg_ratio = int(pari.lift(pari.Mod(r2 / scale, p))) * pow(int(pari.lift(pari.Mod(r1 / scale, p))), -1, p) % p
        L_ratio = (-2 - Lplus_T2_mod_p * pow(Lminus_T2_mod_p, -1, p)) % p
        out["direction_check"] = dict(reg_ratio_mod_p=reg_ratio, L_ratio_mod_p=L_ratio, agree=reg_ratio == L_ratio)
    return out


def split_multiplicative_bound(E, p: int, vp_Lstar_T_r_plus_1: int, rank: int, gens: Sequence, digits: int = 6) -> dict:
    """b_p at a split multiplicative p, with eps_p replaced by L-invariant / log_p(kappa(gamma)) (SW Thm 6.1),
    the weaker of the two printed conventions; the modified height is PARI's (SW section 4.2)."""
    pari = get_pari()
    loc = _local_data(E, p)
    if int(pari.ellap(E, p)) != 1:
        raise ValueError("not split multiplicative")
    Ep = pari.ellinit([int(a) for a in E[:5]], pari(f"O({p}^{digits})"))
    u2, u, q = pari("(E)->E.tate")(Ep)[:3]
    ordq = int(pari.valuation(q, p))
    Linv = pari.log(q) / ordq
    s2 = pari.ellpadics2(E, p, digits)
    H = pari.ellpadicheightmatrix(E, p, digits, pari(gens))
    M = H[0] + H[1] * (-s2 + 1 / pari.log(q) / u2)
    R = pari.matdet(M)
    v_reg_gamma = _v(R, p) - rank
    v_eps_replacement = _v(Linv, p) - 1
    b = vp_Lstar_T_r_plus_1 - v_eps_replacement - loc["sum_vp_c"] + 2 * loc["vp_tors"] - v_reg_gamma
    return dict(p=p, ord_p_q=ordq, vp_L_invariant=_v(Linv, p), vp_Reg_gamma=v_reg_gamma, vp_Lstar=vp_Lstar_T_r_plus_1,
                tamagawa=loc["tamagawa"], b_p=b, b_p_with_L_invariant_as_printed_in_algorithm=b - 1,
                sha_p_trivial=b <= 1, consistent=b >= 0)
