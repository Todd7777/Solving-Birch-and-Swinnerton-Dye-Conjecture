"""Tests for bsd.kato_bound: exact p-adic L-function Riemann sums, Selmer bounds, point independence.

The Riemann sums are compared with PARI's overconvergent-symbol p-adic L-series (mspadicseries), an
independent algorithm, modulo the precision that the error bound (*) guarantees.
"""
from fractions import Fraction

import pytest

from bsd import kato_bound as kb
from bsd.pari_session import get_pari

pari = get_pari()

C11 = [0, -1, 1, -10, -20]
C37 = [0, 0, 1, -1, 0]
C389 = [0, 1, 1, -2, 0]
C5077 = [0, 0, 1, -7, 6]


@pytest.mark.parametrize("N", [1, 2, 11, 12, 37, 360, 389, 5077])
def test_p1_list_and_lifts(N):
    pts = kb.p1_list(N)
    assert len(set(pts)) == len(pts)
    for c, d in pts:
        a, b, c1, d1 = kb.lift_to_sl2z(c, d, N) if N > 1 else (1, 0, 0, 1)
        assert a * d1 - b * c1 == 1
        if N > 1:
            assert (c1 - c) % N == 0 and (d1 - d) % N == 0


def test_plus_symbol_content_and_hecke():
    s11 = kb.PlusSymbol(pari.ellinit(C11))
    assert s11.g == Fraction(1, 10) and s11(0, 1) == 2 and s11.n_manin == 12
    assert s11.hecke_check()
    s389 = kb.PlusSymbol(pari.ellinit(C389))
    assert s389(0, 1) == 0 and s389.hecke_check()          # L(E,1) = 0 for 389a1


def _pari_coeffs(ainvs, p, n=8, jmax=3):
    E = pari.ellinit(ainvs)
    M, phi = pari.msfromell(E, 1)
    mu = pari.mspadicmoments(pari.mspadicinit(M, p, n, 0), phi, 1)
    ser = pari.truncate(pari.mspadicseries(mu))
    return [pari.polcoef(ser, j) for j in range(jmax + 1)]


@pytest.mark.parametrize("ainvs,p", [(C11, 3), (C11, 7), (C37, 5), (C389, 5), (C389, 7)])
def test_riemann_sums_match_pari_overconvergent(ainvs, p):
    E = pari.ellinit(ainvs)
    sym = kb.PlusSymbol(E)
    rc = kb.riemann_coefficients(sym, p, int(pari.ellap(E, p)), 4, 3, P=6)
    for j, (c, e) in enumerate(zip(rc["coeffs"], rc["precision"])):
        e = 6 if e is None else e
        ref = _pari_coeffs(ainvs, p)[j]
        if ref == 0:
            assert c % p ** e == 0
            continue
        assert int(pari.padicprec(ref, p)) >= e
        # ours is PARI's divided by the content g (same orientation for these curves)
        ref_int = int(pari.lift(pari.Mod(pari.truncate(ref) / pari(str(sym.g)), p ** e)))
        assert c % p ** e == ref_int, (j, c, ref)


@pytest.mark.parametrize("ainvs,p,r", [(C11, 3, 0), (C37, 5, 1), (C389, 3, 2), (C389, 5, 2), (C389, 7, 2),
                                       (C5077, 5, 3)])
def test_selmer_bounds(ainvs, p, r):
    res = kb.selmer_bound(ainvs, p)
    assert res["ord_upper_bound"] == r
    # the coefficients below r must not be certified non-zero (Kato: ord_T L_p >= rank)
    assert not any(f["certified_nonzero"] for f in res["facts"][:r])


def test_supersingular_rejected():
    with pytest.raises(ValueError):
        kb.ordinary_bound(pari.ellinit(C11), 19)            # a_19(11a1) = 0


def test_independence_certificate_controls():
    E = pari.ellinit(C389)
    assert kb.independence_certificate(E, [(-1, 1), (0, 0)], 1)["ok"]
    assert not kb.independence_certificate(E, [(-1, 1), (-1, 1)], 1, max_columns=20)["ok"]
    assert not kb.independence_certificate(E, [(-1, 1), (0, 0), (1, 0)], 1, max_columns=20)["ok"]  # P3 = P1+P2
    assert kb.independence_certificate(pari.ellinit(C5077), [(2, 0), (-1, 3), (4, 6)], 1)["ok"]
    assert not kb.independence_certificate(pari.ellinit(C11), [(5, 5)], 5, max_columns=20)["ok"]  # torsion


@pytest.mark.parametrize("ainvs,p,ord_Lp,ord_f", [(C389, 389, 3, 2),     # split: exceptional zero
                                                   (C11, 11, 1, 0),       # split, rank 0 (Greenberg-Stevens)
                                                   (C37, 37, 1, 1)])      # non-split, rank 1
def test_multiplicative_primes(ainvs, p, ord_Lp, ord_f):
    res = kb.multiplicative_bound(pari.ellinit(ainvs), p)
    assert res["ord_upper_bound"] == ord_Lp and res["char_series_ord_upper_bound"] == ord_f
    assert res["rho_p_surjective_by_isogeny_criterion"] == (p >= 5)


# ---------------------------------------------------------------------------------------------
# supersingular primes (a_p = 0): L^-, L^+ from Mazur-Tate elements
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize("ainvs,p", [(C37, 17), (C11, 19)])
def test_supersingular_structure(ainvs, p):
    """Exact checks of the three-term relation and the divisibilities behind L^- == -theta_1 mod omega_1."""
    checks = kb.supersingular_structure_check(kb.PlusSymbol(pari.ellinit(ainvs)), p, 3)
    assert checks and all(checks.values()), checks


@pytest.mark.parametrize("ainvs,p,r", [(C11, 19, 0), (C37, 17, 1), (C37, 19, 1), (C389, 107, 2)])
def test_supersingular_bounds(ainvs, p, r):
    res = kb.supersingular_bound(pari.ellinit(ainvs), p)
    assert res["ord_upper_bound"] == r
    assert all(v == 0 for v in res["Lminus_mod_p"][:r])        # Kato + rank r force these to vanish


def test_supersingular_interpolation_constants():
    """L^-(0) = (p-1) m'(0) and L^+(0) = 2 m'(0) (Pollack's interpolation, as quoted by Skinner) for 11a1, m'(0) = 2."""
    res = kb.supersingular_bound(pari.ellinit(C11), 19, cross_check_plus=True)
    assert res["Lminus_mod_p"][0] == (18 * 2) % 19 and res["Lplus_mod_p"][0] == 4


def test_supersingular_389a1_two_certificates():
    res = kb.supersingular_bound(pari.ellinit(C389), 107, cross_check_plus=True)
    assert res["Lminus_mod_p"][:3] == [0, 0, 21] and res["Lplus_mod_p"][:3] == [0, 0, 76]


def test_stein_wuthrich_1483a1_p5():
    """Stein-Wuthrich, Math. Comp. 82 (2013), section 12.2: 1483a1 (rank 2) at the supersingular prime 5 has
    ord_T L_p = 2, and the displayed series (1-phi)^{-2} L_p has T^2, T^3 coefficients 1 + O(5) (omega part)
    and 4*5 + O(5^2) (phi(omega) part), with their period normalisation (= ours / 2, as for 446d1)."""
    import math
    from fractions import Fraction
    E = pari.ellinit([0, 1, 1, 2, 2])
    sym = kb.PlusSymbol(E)
    assert kb.supersingular_bound(E, 5, sym=sym)["ord_upper_bound"] == 2
    tc = lambda S: [sum(c * math.comb(k, j) for k, c in enumerate(S) if c and k >= j) for j in range(4)]
    G = [Fraction(c, 625) for c in tc(kb.theta_X(sym, 5, 7))]          # G_8 = (-5)^-4 theta_7
    H = [Fraction(c, 625) for c in tc(kb.theta_X(sym, 5, 6))]          # H_7 = (-5)^-4 theta_6
    for j in (2, 3):
        w = Fraction(5, 36) * (4 * G[j] + 10 * H[j]) / 2
        f = Fraction(5, 36) * (10 * G[j] - 20 * H[j]) / 2
        assert w.denominator % 5 and (w.numerator * pow(w.denominator, -1, 5)) % 5 == 1
        f5 = f / 5
        assert f.numerator % 5 == 0 and f5.denominator % 5 and (f5.numerator * pow(f5.denominator, -1, 5)) % 5 == 4


def test_supersingular_rejects_ordinary_prime():
    with pytest.raises(ValueError):
        kb.supersingular_bound(pari.ellinit(C389), 5)               # a_5(389a1) = -3
