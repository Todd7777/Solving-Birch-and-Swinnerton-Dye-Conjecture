"""Fast tests (< 60 s, PARI + python-flint only) for bsd/gz_certificate.py.

The headline test reproduces Buhler-Gross-Zagier (Math. Comp. 44, 1985) for 5077a1:
L'(E,1) = 0 exactly via the Gross-Zagier formula, and L'''(E,1)/3! = 1.7318499001...
"""
import math

import pytest

flint = pytest.importorskip("flint")

from bsd.gz_certificate import (GZCertificateError, certify_rank3, heegner_discriminants,  # noqa: E402
                                is_heegner_discriminant, lprime_ball, minus_modular_symbol,
                                period_area_ball, twisted_minus_symbol_sum, an_list)
from bsd.pari_session import get_pari  # noqa: E402

pari = get_pari()


def test_heegner_discriminants_5077():
    Ds = heegner_discriminants(5077, 3)
    assert Ds == [-7, -19, -47]
    for D in Ds:
        assert int(pari.kronecker(D, 5077)) == 1 and D % 2 != 0 and int(pari.isfundamental(D)) == 1
    assert not is_heegner_discriminant(-8, 5077)          # even
    assert not is_heegner_discriminant(-11, 5077)         # 5077 inert in Q(sqrt -11)


@pytest.mark.parametrize("ainvs,D,expected", [
    ([0, -1, 1, -10, -20], -15, 1),   # 11a1, rank 0
    ([0, -1, 1, -10, -20], -23, 1),
    ([0, 0, 1, -1, 0], -7, 2),        # 37a1, rank 1
    ([0, 0, 1, -7, 6], -7, 8),        # 5077a1, rank 3
])
def test_twisted_symbol_convention(ainvs, D, expected):
    """sum_a (D|a) x^-({oo -> a/|D|}) = L(E^D,1) sqrt|D| / Omega^-  (checked against PARI's ellL1 of the twist)."""
    E = pari.ellinit(ainvs)
    M, xm = minus_modular_symbol(E)
    S = twisted_minus_symbol_sum(M, xm, D)
    assert S == expected
    cinf = 2 if float(pari("(E)->E.disc")(E)) > 0 else 1
    omega_minus = -2.0 / cinf * float(pari.imag(pari("(E)->E.omega")(E)[1]))
    LD = float(pari.ellL1(pari.ellinit(pari.elltwist(E, D)), 0))
    assert abs(LD * math.sqrt(-D) / omega_minus - float(S)) < 1e-10


@pytest.mark.parametrize("ainvs", [[0, 0, 1, -7, 6], [0, -1, 1, -10, -20], [0, 0, 1, -79, 342]])
def test_period_area_ball_contains_pari_area(ainvs):
    vol = period_area_ball(ainvs)
    area = float(pari("(E)->E.area")(pari.ellinit(ainvs)))
    v = vol["value"]
    assert v.overlaps(flint.arb(area) + flint.arb(0, 1) * flint.arb(abs(area) * 1e-13))
    assert v.rad() < 1e-20


def test_lprime_ball_rank_one():
    """For 37a1 the enclosure of L'(E,1) must contain PARI's value and exclude 0."""
    E = pari.ellinit([0, 0, 1, -1, 0])
    an = an_list(E, 2000)
    res = lprime_ball(an, 37, log10_target=-30)
    pv = float(pari.ellL1(E, 1))
    assert res["value"].overlaps(flint.arb(pv) + flint.arb(0, 1) * flint.arb(1e-14))
    assert not res["value"].contains(0)


def test_certify_5077a1_reproduces_bgz():
    cert = certify_rank3([0, 0, 1, -7, 6], height_bound=0.17, candidate_D=2, label="5077a1")
    assert cert["heegner"]["D"] == -7
    assert cert["modular_symbols"]["exact"][0] == {"D": -7, "S_D": "8"}
    assert cert["twist_nonvanishing"]["ok"]
    assert cert["twist_consistency"]["ball_contains_prediction"]
    assert cert["gross_zagier"]["inequality_holds"]
    assert cert["gross_zagier"]["log10_margin"] > 40
    assert cert["manin_constant"]["c_max"] == 1.0
    assert cert["conclusion"]["certified"] and cert["conclusion"]["ord_exact"] == 3
    mid = float(cert["third_derivative"]["ball"]["mid"])
    assert abs(mid - 1.7318499001193006898) < 1e-15


def test_negative_control_rank_one_curve_is_not_certified():
    """37a1 has L'(E,1) != 0 and a non-torsion Heegner point: the inequality must fail."""
    cert = certify_rank3([0, 0, 1, -1, 0], height_bound=0.04, candidate_D=2, want_third=False)
    assert not cert["gross_zagier"]["inequality_holds"]
    assert not cert["conclusion"]["certified"]


def test_rejects_root_number_plus_one():
    with pytest.raises(GZCertificateError):
        certify_rank3([0, 1, 1, -2, 0], height_bound=0.1)   # 389a1, w = +1
