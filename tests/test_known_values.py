"""Regression tests against independently known values (Cremona's tables, LMFDB).

Run:  python -m pytest -q
"""
import math
from fractions import Fraction

import pytest

from bsd.invariants import compute_invariants
from bsd.analytic import analytic_data, exact_L1_ratio, certify_order
from bsd.descent import mordell_weil
from bsd.verify import verify_curve
from bsd.theorems import p_part_status, sha_p_finite_via_kato
from bsd.pari_session import get_pari

CURVES = {
    # label: (ainvs, conductor, rank, torsion, tamagawa product, sha_an, regulator)
    "11a1": ([0, -1, 1, -10, -20], 11, 0, 5, 5, 1, 1.0),
    "37a1": ([0, 0, 1, -1, 0], 37, 1, 1, 1, 1, 0.0511114082399688),
    "389a1": ([0, 1, 1, -2, 0], 389, 2, 1, 1, 1, 0.152460177943144),
    "5077a1": ([0, 0, 1, -7, 6], 5077, 3, 1, 1, 1, 0.417143558758384),
    "571a1": ([0, -1, 1, -929, -10595], 571, 0, 1, 1, 4, 1.0),
    "681b1": ([1,1,0,-1154,-15345], 681, 0, 4, 4, 9, 1.0),
}


@pytest.mark.parametrize("label", list(CURVES))
def test_invariants(label):
    ainvs, N, r, tors, tam, sha, reg = CURVES[label]
    inv = compute_invariants(ainvs)
    assert inv.conductor == N
    assert inv.torsion_order == tors
    assert inv.tamagawa_product == tam
    assert inv.root_number == (-1) ** r


@pytest.mark.parametrize("label", list(CURVES))
def test_full_bsd_quotient(label):
    ainvs, N, r, tors, tam, sha, reg = CURVES[label]
    rec = verify_curve(ainvs, label=label)
    assert rec.analytic["analytic_rank"] == r
    assert rec.descent["rank_lower"] == r and rec.descent["rank_upper"] == r
    assert abs(rec.descent["regulator"] - reg) < 1e-9
    assert rec.sha_an_int == sha
    assert rec.sha_is_square
    assert rec.sha_rel_err < 1e-9
    assert rec.errors == []


def test_exact_L_ratio_11a1():
    pari = get_pari()
    E = pari.ellinit([0, -1, 1, -10, -20])
    assert exact_L1_ratio(E) == Fraction(1, 5)


def test_exact_vanishing_389a1():
    pari = get_pari()
    E = pari.ellinit([0, 1, 1, -2, 0])
    assert exact_L1_ratio(E) == 0
    lower, exact, _ = certify_order(2, +1, Fraction(0))
    assert lower == 2 and exact


def test_rank3_not_certifiable():
    # root number -1, L(E,1)=0: only ord >= 1 is rigorous
    lower, exact, msg = certify_order(3, -1, Fraction(0))
    assert lower == 1 and not exact


def test_period_convention_37a1():
    inv = compute_invariants([0, 0, 1, -1, 0])
    assert inv.real_components == 2
    assert abs(inv.omega - 5.98691729246392) < 1e-10
    assert abs(inv.omega_plus - 2.99345864623196) < 1e-10


def test_saturation_fixes_index():
    # ellrank returns an index-3 sublattice for 5077a1; saturation must repair it
    pari = get_pari()
    E = pari.ellinit([0, 0, 1, -7, 6])
    d = mordell_weil(E, known_gens=None, saturation_bound=100)
    assert abs(d.regulator - 0.417143558758384) < 1e-9


def test_p_part_ledger_skinner_urban_11a1():
    inv = compute_invariants([0, -1, 1, -10, -20])
    # 11a1 has a rational 5-isogeny, so E[5] is reducible: p = 5 not covered
    assert p_part_status(inv, 0, 5).status == "not covered"
    # p = 7: good ordinary (a_7 = -2), irreducible, ramified at 11 (ord_11 Δ = 5, 7 ∤ 5)
    assert p_part_status(inv, 0, 7).status == "proven"


def test_kato_finiteness_logic():
    inv = compute_invariants([0, 1, 1, -2, 0])   # 389a1
    assert sha_p_finite_via_kato(inv, rank_lower=2, p=5, r_p=2)["sha_p_finite"]
    assert not sha_p_finite_via_kato(inv, rank_lower=2, p=5, r_p=3)["sha_p_finite"]
    assert not sha_p_finite_via_kato(inv, rank_lower=2, p=389, r_p=2)["sha_p_finite"]
