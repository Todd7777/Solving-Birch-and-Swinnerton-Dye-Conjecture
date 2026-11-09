#!/usr/bin/env python3
"""Certify that the 3-adic Galois representation of 389a1 is surjective onto GL_2(Z_3).

Stein-Wuthrich Thm 7.3 needs rho_3 surjective at p = 3. The argument:

 (1) rhobar_3 is surjective. The 3-division polynomial psi_3 is irreducible with Galois group S_4, so
     the projective image is PGL_2(F_3). A subgroup of GL_2(F_3) with that projective image and
     surjective determinant is GL_2(F_3), because SL_2(F_3) has no subgroup of index 2 (checked below).
 (2) If the image H_9 of rho_9 were proper, then K = [H_9, H_9] would map onto [GL_2(F_3), GL_2(F_3)] =
     SL_2(F_3) (checked). K meets V = ker(SL_2(Z/9) -> SL_2(Z/3)) in an SL_2(F_3)-submodule of
     V = sl_2(F_3), which is irreducible (checked), so in 1 or V. With V, K = SL_2(Z/9) and H_9 = GL_2(Z/9).
     With 1, K is a lift of SL_2(F_3), hence conjugate to Elkies' group G (all lifts are conjugate:
     checked, as Elkies states). The same argument for H_9 n SL_2 gives H_9 n SL_2 = G, so H_9 lies in
     N(G), which is Elkies' G' (checked: N(G) n SL_2 = G and det(N(G)) = (Z/9)^x). Then j(E) = f(x)
     for some x in P^1(Q) (moduli interpretation of X_9 = X(9)/G', Elkies 2006; j(E) != 0, 1728).
     The script checks that j(389a1) = 1404928/389 is not a value of f on P^1(Q) (exact factorisation).
     So H_9 = GL_2(Z/9).
 (3) Lifting lemma: a closed subgroup of GL_2(Z_3) with surjective determinant that maps onto GL_2(Z/9)
     is GL_2(Z_3). Proof: the closure C of its commutator subgroup maps onto [GL_2(Z/9), GL_2(Z/9)] =
     SL_2(Z/9) (checked). For n >= 1 and X of trace 0, (I + 3^n X)^3 = I + 3^{n+1} X mod 3^{n+2}, so
     C contains every level of the congruence filtration of SL_2(Z_3); C is closed, so C = SL_2(Z_3).

Trusted input not re-derived here: the moduli interpretation in (2), Elkies' rational function f
(transcribed from arXiv:math/0612734 and checked against both of his forms and his table of values),
and PARI's polgalois and factor.
"""
from __future__ import annotations

import itertools
import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))


def mul(A, B, m):
    a, b, c, d = A
    e, f, g, h = B
    return ((a * e + b * g) % m, (a * f + b * h) % m, (c * e + d * g) % m, (c * f + d * h) % m)


def det(A, m):
    a, b, c, d = A
    return (a * d - b * c) % m


def inv(A, m):
    a, b, c, d = A
    di = pow(det(A, m), -1, m)
    return ((d * di) % m, (-b * di) % m, (-c * di) % m, (a * di) % m)


def gl2(m):
    return [A for A in itertools.product(range(m), repeat=4) if det(A, m) % 3 != 0]


def closure(gens, m):
    I = (1, 0, 0, 1)
    seen = {I}
    frontier = [I]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = mul(x, g, m)
                if y not in seen:
                    seen.add(y)
                    new.append(y)
        frontier = new
    return seen


def commutator_subgroup(Hs, m):
    comms = {mul(mul(inv(a, m), inv(b, m), m), mul(a, b, m), m) for a in Hs for b in Hs}
    return closure(list(comms), m)


def red(A, m):
    return tuple(x % m for x in A)


def main():
    GL3, GL9 = gl2(3), gl2(9)
    SL3 = {A for A in GL3 if det(A, 3) == 1}
    SL9 = {A for A in GL9 if det(A, 9) == 1}
    assert len(GL3) == 48 and len(SL3) == 24 and len(GL9) == 3888 and len(SL9) == 648
    checks = {}
    # SL_2(F_3) has no subgroup of index 2: such a subgroup would contain all squares
    squares = closure([mul(A, A, 3) for A in SL3], 3)
    checks["SL2(F3) has no index-2 subgroup"] = len(squares) == 24
    checks["[GL2(F3),GL2(F3)] = SL2(F3)"] = commutator_subgroup(GL3, 3) == SL3
    checks["[GL2(Z/9),GL2(Z/9)] = SL2(Z/9)"] = commutator_subgroup(GL9, 9) == SL9
    # V = kernel of reduction on SL_2(Z/9); irreducible under conjugation by SL_2(Z/9)
    V = [A for A in SL9 if red(A, 3) == (1, 0, 0, 1)]
    assert len(V) == 27
    invariant = []
    for g in V:
        if g == (1, 0, 0, 1):
            continue
        sub = closure([g], 9)
        for h in V:                                      # subgroups of order 9 = <g, h>
            sub2 = closure([g, h], 9)
            for S_ in {frozenset(sub), frozenset(sub2)}:
                if 1 < len(S_) < 27 and all(mul(mul(x, s, 9), inv(x, 9), 9) in S_ for x in SL9 for s in S_):
                    invariant.append(S_)
    checks["V = sl2(F3) irreducible (no proper non-trivial invariant subgroup)"] = not invariant
    # lifts of SL_2(Z/3) = <S, T | S^2 = (ST)^3 = -I, T^3 = I> into SL_2(Z/9), standard S = [[0,-1],[1,0]]
    # (Elkies' printed S-tilde = [[0,2],[4,0]] reduces to this S mod 3)
    S3, T3 = (0, 2, 1, 0), (1, 1, 0, 1)
    I9, mI9 = (1, 0, 0, 1), (8, 0, 0, 8)
    liftsS = [A for A in SL9 if red(A, 3) == S3 and mul(A, A, 9) == mI9]
    liftsT = [A for A in SL9 if red(A, 3) == T3 and mul(mul(A, A, 9), A, 9) == I9]
    pairs, lifts = [], set()
    for s in liftsS:
        for t in liftsT:
            st = mul(s, t, 9)
            if mul(mul(st, st, 9), st, 9) == mI9:
                H = frozenset(closure([s, t], 9))
                if len(H) == 24 and len({red(A, 3) for A in H}) == 24:
                    pairs.append((s, t))
                    lifts.add(H)
    checks["number of lift pairs (S~, T~) (Elkies: 27)"] = len(pairs)
    checks["Elkies' S~ = [[0,2],[4,0]] is one of the lifts of S"] = any(s == (0, 2, 4, 0) for s, _ in pairs)
    G = set(next(iter(lifts)))
    conjG = {frozenset(mul(mul(x, g, 9), inv(x, 9), 9) for g in G) for x in SL9}
    checks["every lift of SL2(Z/3) into SL2(Z/9) is conjugate to G (Elkies: all 27 conjugate)"] = (
        bool(lifts) and lifts <= conjG)
    checks["number of distinct lifted subgroups"] = len(lifts)
    Gl = list(G)
    N = {x for x in GL9 if all(mul(mul(x, g, 9), inv(x, 9), 9) in G for g in Gl)}
    checks["N(G) n SL2 = G and det(N(G)) = (Z/9)^x, so N(G) is the extension G' of Elkies"] = (
        {x for x in N if det(x, 9) == 1} == G and {det(x, 9) for x in N} == {1, 2, 4, 5, 7, 8})
    # j(389a1) is not a value of Elkies' f on P^1(Q)
    from bsd.pari_session import get_pari
    pari = get_pari()
    num = pari("3^7*(x^2-1)^3*(x^6+3*x^5+6*x^4+x^3-3*x^2+12*x+16)^3*(2*x^3+3*x^2-3*x-5)")
    den = pari("(x^3-3*x-1)^9")
    A_ = pari("x^6+6*x^5+4*x^3+12*x^2-18*x-23")
    B_ = pari("7*x^6+24*x^5+18*x^4-26*x^3-33*x^2+18*x+28")
    checks["Elkies' two forms of f agree"] = num == 1728 * den + 27 * A_ ** 2 * B_ ** 2 * pari("2*x^3-3*x^2+4")
    table = {"-2": 419904, "0": -44789760, "-1/2": 15786448344, "2": 24992518538304, "-3/2": -92515041526500,
             "-1/3": -70043919611288518656, "1": 0, "-1": 0}
    checks["Elkies' table of values reproduced"] = all(
        pari.subst(num, "x", pari(x)) / pari.subst(den, "x", pari(x)) == v for x, v in table.items()) and \
        pari.polcoef(num, 27) / pari.polcoef(den, 27) == 4374
    E = pari.ellinit([0, 1, 1, -2, 0])
    j = pari("(E)->E.j")(E)
    fa = pari.factor(num - j * den)[0]
    roots = [g for g in fa if int(pari.poldegree(g)) == 1]
    checks["j(389a1) = 1404928/389 has no rational preimage under f (and f(oo) = 4374 != j)"] = not roots and j != 4374
    g3 = pari.polgalois(pari.elldivpol(E, 3))
    checks["psi_3 irreducible with Galois group S4"] = bool(pari.polisirreducible(pari.elldivpol(E, 3))) and str(g3[3]) == "S4"
    for k, v in checks.items():
        print(f"{'OK  ' if v is True or (isinstance(v, int) and not isinstance(v, bool)) else 'FAIL'} {k}: {v}")
    ok = all(v for k, v in checks.items() if isinstance(v, bool))
    print("CONCLUSION:", "rho_3 of 389a1 is surjective onto GL_2(Z_3)" if ok else "NOT established")
    import json
    with open(os.path.join(HERE, "rho3_certificate.json"), "w") as fh:
        json.dump(dict(curve="389a1", rho3_surjective=ok, checks={k: v for k, v in checks.items()},
                       trusted=["moduli interpretation of X_9 = X(9)/G' (Elkies 2006)",
                                "Elkies' rational function f (transcribed, checked against his two forms and table)",
                                "PARI polgalois, factor"]), fh, indent=1)
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
