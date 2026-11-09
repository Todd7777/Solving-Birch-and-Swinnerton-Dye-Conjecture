"""Tests for bsd.bsd1 (certified rank = analytic rank), including negative controls."""
from fractions import Fraction

import pytest

from bsd import bsd1
from bsd.pari_session import get_pari

pari = get_pari()

CASES = {
    # label: (ainvs, rank, generators, torsion order)
    "11a1": ([0, -1, 1, -10, -20], 0, [], 5),
    "37a1": ([0, 0, 1, -1, 0], 1, [(0, 0)], 1),
    "389a1": ([0, 1, 1, -2, 0], 2, [(-1, 1), (0, 0)], 1),
    "5077a1": ([0, 0, 1, -7, 6], 3, [(-2, 3), (-1, 3), (0, 2)], 1),
}


@pytest.mark.parametrize("label", list(CASES))
def test_certify_known_classes(label):
    ainvs, r, gens, tors = CASES[label]
    gens = [(Fraction(x), Fraction(y)) for x, y in gens]
    res = bsd1.certify_class(label, ainvs, r, gens, tors, kato_max_conductor=6000)
    assert res["bsd1_certified"], res["problems"]
    assert res["analytic"]["order"] == r
    if r >= 2:
        assert res["descent"]["unconditional"] and res["descent"]["rank_upper"] == r
        assert res["independence"]["ok"]
        assert res["kato"]["ok"] and res["kato"]["p"] >= 5 and res["kato"]["ord_upper_bound"] == r
        assert {s["p"] for s in res["sha_p_finite"]} >= {2}


def test_leading_coefficients():
    ainvs, r, gens, tors = CASES["5077a1"]
    res = bsd1.certify_class("5077a1", ainvs, r, [(Fraction(x), Fraction(y)) for x, y in gens], tors)
    mid = float(res["analytic"]["leading_ball"].strip("[").split("+/-")[0])
    assert abs(mid - 1.7318499001193007) < 1e-15                # L'''(E,1)/3!, Buhler-Gross-Zagier


def test_wrong_rank_claims_are_not_certified():
    ainvs, _, gens, tors = CASES["389a1"]
    gens = [(Fraction(x), Fraction(y)) for x, y in gens]
    # rank 0 claimed: parity fits (w = +1), but the ball for L(E,1) contains 0
    res = bsd1.certify_class("389a1", ainvs, 0, [], tors)
    assert not res["bsd1_certified"] and "order of vanishing not certified" in res["problems"]
    # rank 1 claimed: parity is wrong
    res = bsd1.certify_class("389a1", ainvs, 1, gens[:1], tors)
    assert not res["bsd1_certified"]
    # rank 2 claimed with a dependent pair of points
    res = bsd1.certify_class("389a1", ainvs, 2, [gens[0], gens[0]], tors)
    assert not res["bsd1_certified"] and "independence of the generators not certified" in res["problems"]
    # 37a1 with rank 3 claimed: only one generator exists
    res = bsd1.certify_class("37a1", [0, 0, 1, -1, 0], 3, [(Fraction(0), Fraction(0))], 1)
    assert not res["bsd1_certified"]


def test_point_not_on_curve_is_rejected():
    ainvs, r, gens, tors = CASES["389a1"]
    res = bsd1.certify_class("389a1", ainvs, r, [(Fraction(-1), Fraction(1)), (Fraction(1), Fraction(1))], tors)
    assert not res["bsd1_certified"] and "a generator is not on the curve" in res["problems"]


def test_descent_bound_rational_two_torsion_needs_no_class_group():
    # 960d1-like: a curve with full rational 2-torsion uses only Q (no cubic field to certify)
    E = pari.ellinit([0, 0, 0, -1, 0])                       # y^2 = x^3 - x, rank 0
    d = bsd1.descent_bound(E)
    assert d["cubic_fields"] == [] and d["unconditional"] and d["rank_upper"] == 0
