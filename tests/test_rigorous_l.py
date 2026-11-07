"""Tests for bsd.rigorous_l: rigorous Arb enclosures of L^{(k)}(E,1)/k! and certified orders.

PARI's values (lfun at 320 bits) are floating-point references, not rigorous.  The tests
check that they lie inside our balls and that the balls are tiny.  Run:
python -m pytest -q tests/test_rigorous_l.py   (a few seconds)
"""
import math

import pytest
from flint import acb, arb, ctx

from bsd.pari_session import get_pari
from bsd.rigorous_l import (C_series, C_series_rgamma, K_values, _dec_upper, _ntail_bound,
                            ball_to_json, certify_order_rigorous, lambda_by_integration,
                            taylor_enclosures)

CURVES = {
    "11a1": [0, -1, 1, -10, -20],
    "37a1": [0, 0, 1, -1, 0],
    "389a1": [0, 1, 1, -2, 0],
    "5077a1": [0, 0, 1, -7, 6],
}
PARI_BITS = 320


def _pari_c(E, k):
    pari = get_pari()
    v = pari.lfun(E, 1, k, precision=PARI_BITS)
    with ctx.workprec(PARI_BITS):
        return arb(str(v).replace(" E", "e")) / math.factorial(k)


def _with_old_precision(pari, digits):
    old = int(pari.default("realprecision"))
    pari.set_real_precision(digits)
    return old


@pytest.fixture(scope="module")
def pari():
    p = get_pari()
    old = _with_old_precision(p, 110)      # so that str() of a PARI real shows ~110 digits
    yield p
    p.set_real_precision(old)


@pytest.mark.parametrize("label", list(CURVES))
def test_enclosure_contains_pari(pari, label):
    ainvs = CURVES[label]
    enc = taylor_enclosures(ainvs, kmax=4, prec=128)
    E = pari.ellinit(ainvs)
    for k in range(5):
        pv = _pari_c(E, k)
        ck = enc.c[k]
        if ck.is_zero():                        # exact zero (root number -1, k = 0)
            assert enc.root_number == -1 and k == 0
            assert abs(float(pv.mid())) < 1e-60
            continue
        assert ck.contains(pv), (label, k, ck.str(40), pv.str(40))
        assert float(ck.rad()) < 1e-30


def test_lambda_matches_lfunlambda(pari):
    # PARI's Lambda carries the gamma factor Gamma_R(s)Gamma_R(s+1) = 2 (2 pi)^{-s} Gamma(s): twice ours.
    for label in ("11a1", "389a1"):
        enc = taylor_enclosures(CURVES[label], kmax=4, prec=128)
        E = pari.ellinit(CURVES[label])
        for k in range(5):
            v = pari.lfunlambda(E, 1, k, precision=PARI_BITS)
            with ctx.workprec(PARI_BITS):
                lv = arb(str(v).replace(" E", "e")) / math.factorial(k) / 2
            if enc.lam[k].is_zero():
                assert (-1) ** k == -enc.root_number and abs(float(lv.mid())) < 1e-60
            else:
                assert enc.lam[k].contains(lv)


def test_certified_orders():
    expected = {"11a1": 0, "37a1": 1, "389a1": 2}
    for label, r in expected.items():
        res = certify_order_rigorous(CURVES[label], kmax=4, label=label)
        assert res["order_certified"], res["conclusion"]
        assert res["order"] == r
        assert r in res["c_nonzero_certified"]
        assert [int(k) for k in res["c_zero_exact"]] == list(range(r))
    res389 = certify_order_rigorous(CURVES["389a1"], kmax=4)
    assert res389["modular_symbols"]["zero"] is True          # L(E,1) = 0 exactly, from modular symbols
    assert "parity" in res389["c_zero_exact"]["1"]


def test_rank3_gives_bound_not_vanishing():
    res = certify_order_rigorous(CURVES["5077a1"], kmax=4)
    assert not res["order_certified"]
    assert res["order_candidates"] == [1, 3]
    assert 3 in res["c_nonzero_certified"]
    assert float(res["L_prime_abs_upper"]) < 1e-30          # rigorous |L'(E,1)| bound, not a proof of 0
    assert res["c"][0]["exact_zero"]                          # c_0 = 0 by parity


def test_external_certificate_hook():
    # With an externally proven lower bound ord >= 3, for instance from Gross-Zagier-Kolyvagin
    # plus a proof that rank >= 2, the rigorous c_3 != 0 closes the gap: ord = 3.
    res = certify_order_rigorous(CURVES["5077a1"], kmax=4,
                                 external_lower_bound=(3, "test input: assumed proven elsewhere"))
    assert res["order_certified"] and res["order"] == 3
    assert "external" in res["c_zero_exact"]["1"] and res["external_inputs"]["lower_bound"]["ord_at_least"] == 3
    # An external claim contradicted by an enclosure is rejected: L'(E,1) != 0 for 37a1.
    with pytest.raises(RuntimeError):
        certify_order_rigorous(CURVES["37a1"], kmax=2, external_zeros={1: "bogus"})


def test_without_modular_symbols_rank2_is_not_certified():
    res = certify_order_rigorous(CURVES["389a1"], kmax=4, use_modular_symbols=False)
    assert not res["order_certified"] and res["order_candidates"] == [0, 2]


def test_K_values_against_closed_forms_and_integration():
    with ctx.workprec(200):
        for xs in ("0.04", "1.3", "25"):
            x = arb(xs)
            vals = K_values(x, [0, 1, 2, 3, 4], 200)
            assert vals[0].overlaps((-x).exp() / x)                 # J_0 = e^{-x}/x
            assert vals[1].overlaps(x.expint(1) / x)                # J_1 = E_1(x)/x
    with ctx.workprec(100):
        x = arb(2)
        vals = K_values(x, [2, 3, 4], 100)
        for k in (2, 3, 4):
            integral = acb.integral(lambda u, an, k=k: (-2 * u).exp() * u.log(analytic=an) ** k, 1, 50).real
            tail = (x * 50).gamma_upper(k + 1) / x ** (k + 1)       # >= int_50^oo e^{-2u}(log u)^k du
            assert vals[k].overlaps((integral + arb(0, tail)) / math.factorial(k))


def test_two_independent_rigorous_methods_agree():
    # K_k series machinery vs Arb's certified integrator on int_1^oo g(u)(log u)^k du.
    for label in ("11a1", "37a1", "389a1"):
        enc = taylor_enclosures(CURVES[label], kmax=4, prec=128)
        lam_int = lambda_by_integration(CURVES[label], kmax=4, prec=128)
        for k in range(5):
            assert lam_int[k].overlaps(enc.lam[k]), (label, k)
            if not enc.lam[k].is_zero():
                assert float(lam_int[k].rad()) < 1e-30


def test_C_series_two_ways():
    for N in (11, 389, 19999):
        A, B = C_series(N, 6, 200), C_series_rgamma(N, 6, 200)
        assert all(a.overlaps(b) for a, b in zip(A, B))


def test_ntail_bound_dominates_partial_tail(pari):
    # (B3) bounds sum_{n>M} |a_n| K_k(x_n); compare against the next 400 terms.
    ainvs, N, M = CURVES["37a1"], 37, 20
    an = [int(v) for v in pari.ellan(pari.ellinit(ainvs), M + 400)]
    with ctx.workprec(128):
        a = 2 * arb.pi() / arb(N).sqrt()
        for k in range(5):
            partial = sum((abs(an[n - 1]) * K_values(a * n, [k], 128)[k] for n in range(M + 1, M + 401)), arb(0))
            assert partial < _ntail_bound(a, M, k)


def test_model_invariance():
    ai = CURVES["389a1"]
    u = 3
    ai2 = [ai[0] * u, ai[1] * u ** 2, ai[2] * u ** 3, ai[3] * u ** 4, ai[4] * u ** 6]
    e1, e2 = taylor_enclosures(ai, 3), taylor_enclosures(ai2, 3)
    assert e2.conductor == 389 and all(x.overlaps(y) for x, y in zip(e1.c, e2.c))


def test_json_bounds_are_upper_bounds():
    with ctx.workprec(128):
        x = arb("-0.000123456789123456789", "3e-30")
        j = ball_to_json(x)
        assert arb(j["abs_upper"]) >= abs(x).upper()
        assert arb(j["rad"]) >= x.rad()
        assert not j["exact_zero"] and not j["contains_zero"]
        assert _dec_upper(arb(0)) == "0"
