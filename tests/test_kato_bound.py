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
