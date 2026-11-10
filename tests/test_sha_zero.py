"""Tests for bsd.sha_zero: Lemma A by brute force over GL_2(F_l), the surjectivity and saturation certificates with
negative controls, and the bounds b_p against the worked examples of Stein-Wuthrich, Math. Comp. 82 (2013):
446d1 at p = 5 (anomalous ordinary, their (4.2) and section 7.1), at p = 223 (split multiplicative, section 12.1),
and 1483a1 at p = 5 (supersingular, section 12.2). Their period normalisation is ours divided by 2."""
import random

import pytest

from bsd import kato_bound as kb
from bsd import sha_zero as sz
from bsd.pari_session import get_pari

pari = get_pari()

C389, G389 = [0, 1, 1, -2, 0], [[-1, 1], [0, 0]]
E0, G0 = [1, -1, 0, -4, 4], [[2, 0], [1, -1]]          # 446d1 with SW's generators (their section 3.3)
C1483, G1483 = [0, 1, 1, 2, 2], [[-1, 0], [0, 1]]      # 1483a1 with SW's generators (their section 12.2)


def _same(x, ref: str, n: int, p: int) -> bool:
    """x == ref modulo p^n, for p-adic numbers given as PARI objects / strings."""
    return int(pari.valuation(x - pari(ref), p)) >= n


# ---------------------------------------------------------------------------------------------
# Lemma A
# ---------------------------------------------------------------------------------------------
def _mul(a, b, l):
    return ((a[0] * b[0] + a[1] * b[2]) % l, (a[0] * b[1] + a[1] * b[3]) % l,
            (a[2] * b[0] + a[3] * b[2]) % l, (a[2] * b[1] + a[3] * b[3]) % l)


def _closure(gens, l):
    seen, todo = set(gens), list(gens)
    while todo:
        x = todo.pop()
        for g in gens:
            y = _mul(x, g, l)
            if y not in seen:
                seen.add(y)
                todo.append(y)
    return seen


def _fixes_a_line(gens, l):
    lines = [(1, a) for a in range(l)] + [(0, 1)]
    for v in lines:
        if all((m[0] * v[0] + m[1] * v[1]) * v[1] % l == (m[2] * v[0] + m[3] * v[1]) * v[0] % l for m in gens):
            return True
    return False


@pytest.mark.parametrize("l", [3, 5])
def test_lemma_A_dichotomy_bruteforce(l):
    """Every subgroup <t, g> of GL_2(F_l) with t a transvection fixes a line or contains SL_2(F_l)."""
    t = (1, 1, 0, 1)
    n_sl2 = l * (l * l - 1)
    for a in range(l ** 4):
        g = (a % l, a // l % l, a // l ** 2 % l, a // l ** 3)
        if (g[0] * g[3] - g[1] * g[2]) % l == 0:
            continue
        G = _closure([t, g], l)
        sl2 = sum(1 for m in G if (m[0] * m[3] - m[1] * m[2]) % l == 1)
        assert _fixes_a_line([t, g], l) or sl2 == n_sl2


def test_lemma_A_irreducible_frobenius_l7():
    """With g of irreducible characteristic polynomial (the certificate's Frobenius), <t, g> contains SL_2(F_7)."""
    l, t = 7, (1, 1, 0, 1)
    irr = [g for g in ((a, b, c, d) for a in range(l) for b in range(l) for c in range(l) for d in range(l))
           if (g[0] * g[3] - g[1] * g[2]) % l and pow(((g[0] + g[3]) ** 2 - 4 * (g[0] * g[3] - g[1] * g[2])) % l,
                                                      (l - 1) // 2, l) == l - 1]
    assert len(irr) == 21 * 42                                   # 21 irreducible quadratics, 42 matrices each
    for g in random.Random(0).sample(irr, 40):
        G = _closure([t, g], l)
        assert sum(1 for m in G if (m[0] * m[3] - m[1] * m[2]) % l == 1) == l * (l * l - 1)


# ---------------------------------------------------------------------------------------------
# surjectivity and saturation certificates
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize("ell", [3, 5, 7, 107, 2957, 4999])
def test_rhobar_certificate_389a1(ell):
    c = sz.rhobar_surjective_certificate(pari.ellinit(C389), ell)
    assert c["ok"] and c["transvection"]["q"] == 389 and not c["transvection"]["q_equals_ell"]
    q, aq = c["frobenius"]["q"], c["frobenius"]["a_q"]
    assert pow((aq * aq - 4 * q) % ell, (ell - 1) // 2, ell) == ell - 1


def test_rhobar_certificate_389a1_at_389_uses_inertia_at_389():
    c = sz.rhobar_surjective_certificate(pari.ellinit(C389), 389)
    assert c["ok"] and c["transvection"] == dict(q=389, ord_q_Delta_min=1, reduction="multiplicative",
                                                 q_equals_ell=True)


def test_rhobar_certificate_negative_controls():
    # 11a1: Delta = -11^5, so 5 | ord_11(Delta): no transvection is available (rhobar_5 is reducible)
    assert not sz.rhobar_surjective_certificate(pari.ellinit([0, -1, 1, -10, -20]), 5)["ok"]
    # 11a3: ord_11(Delta) = 1 gives a transvection, but E(Q)[5] != 0 puts the image in a Borel subgroup, so every
    # Frobenius has a split characteristic polynomial mod 5 and the search must fail
    c = sz.rhobar_surjective_certificate(pari.ellinit([0, -1, 1, 0, 0]), 5, qmax=3000)
    assert not c["ok"] and "Frobenius" in c["reason"]


@pytest.mark.parametrize("p", [3, 5, 107, 389, 2957])
def test_saturation_389a1(p):
    c = sz.saturation_certificate(pari.ellinit(C389), G389, p)
    assert c["ok"] and len(c["primes"]) >= 2
    for u in c["primes"]:
        assert u["card"] % p == 0 and u["card"] % (p * p) != 0


def test_saturation_negative_controls():
    E = pari.ellinit(C389)
    P1 = pari.ellmul(E, G389[0], 3)
    bad = sz.saturation_certificate(E, [P1, G389[1]], 3)
    assert not bad["ok"] and bad["not_saturated"] and bad["witness"][0] in ([1, 0], [-1, 0])
    assert sz.saturation_certificate(E, [P1, G389[1]], 5)["ok"]
    # ellrank's points for 446d1 have index 3 in E(Q); SW's generators are a basis
    e0 = pari.ellinit(E0)
    assert sz.saturation_certificate(e0, [[-2, 2], [-1, 3]], 3).get("not_saturated")
    assert sz.saturation_certificate(e0, G0, 3)["ok"]


# ---------------------------------------------------------------------------------------------
# good ordinary p (Lemma B)
# ---------------------------------------------------------------------------------------------
@pytest.mark.parametrize("p,vN,vR", [(5, 0, 2), (859, 1, 0)])
def test_ordinary_bound_389a1(p, vN, vR):
    b = sz.ordinary_bound(pari.ellinit(C389), p, 0, 2, G389)
    assert (b["vp_N_p"], b["pari_vp_Reg_p"]) == (vN, vR)
    assert b["b_p"] == 0 and b["b_p_from_pari_regulator"] == 0 and b["sha_p_trivial"] and b["consistent"]


def test_ordinary_bound_stein_wuthrich_446d1_p5():
    """SW section 7.1: 5 is anomalous for 446d1 (N_5 = 10), [T^2]L_5 has valuation 1, Sha(5) = 0. Lemma B(ii)
    gives b_5 = 1 without the regulator; PARI's regulator is SW's (4.2) digit for digit and gives b_5 = 0."""
    E = pari.ellinit(E0)
    rc = kb.riemann_coefficients(kb.PlusSymbol(E), 5, int(pari.ellap(E, 5)), 4, 3, P=4)
    assert rc["precision"][2] >= 2 and kb.vp(rc["coeffs"][2], 5) == 1
    assert _same(pari.ellpadicregulator(E, 5, 10, pari(G0)), "2*5 + 2*5^2 + 5^4 + 4*5^5 + 2*5^7", 8, 5)
    b = sz.ordinary_bound(E, 5, 1, 2, G0)
    assert b["vp_N_p"] == 1 and b["vp_eps"] == 2 and b["b_p"] == 1 and b["sha_p_trivial"]
    assert b["pari_vp_Reg_p"] == 1 and b["b_p_from_pari_regulator"] == 0 and b["consistent"]


def test_ordinary_bound_rejects():
    with pytest.raises(ValueError):
        sz.ordinary_bound(pari.ellinit(C389), 107, 0, 2)          # supersingular
    with pytest.raises(ValueError):
        sz.ordinary_bound(pari.ellinit(C389), 389, 0, 2)          # bad


# ---------------------------------------------------------------------------------------------
# supersingular p (D_p-valued)
# ---------------------------------------------------------------------------------------------
def test_supersingular_bound_389a1_p107():
    E = pari.ellinit(C389)
    b = sz.supersingular_bound(E, 107, 21, 2, G389, Lplus_T2_mod_p=76)    # [T^2]L^-, [T^2]L^+ mod 107 (T2 data)
    assert b["b_p"] == 0 and b["vp_Reg_p"] == 2 and b["lattices_omega_eta_equal"] and b["sha_p_trivial"]
    assert b["direction_check"] == dict(reg_ratio_mod_p=81, L_ratio_mod_p=81, agree=True)
    # swapping L^- and L^+ breaks the proportionality: the check is not vacuous
    assert not sz.supersingular_bound(E, 107, 76, 2, G389, Lplus_T2_mod_p=21)["direction_check"]["agree"]


def test_supersingular_bound_stein_wuthrich_1483a1_p5():
    """SW section 12.2: the normalised D_p-valued regulator of 1483a1 at 5 is
    (1 + 2*5 + 3*5^2 + 5^3 + O(5^5)) omega + (4*5 + 4*5^2 + 4*5^3 + 5^4 + 2*5^5 + O(5^6)) phi(omega), and
    Sha(5) = 0. PARI (basis omega, F omega, F = 5 phi) agrees to O(5^6) in the phi(omega)-coordinate and to O(5^4)
    in the omega-coordinate; there PARI has 4*5^4 at every precision tried, where SW print O(5^5)."""
    E = pari.ellinit(C1483)
    R = pari.ellpadicregulator(E, 5, 10, pari(G1483))
    lg2 = pari.log(pari("6 + O(5^12)")) ** 2
    assert _same(R[0] / lg2, "1 + 2*5 + 3*5^2 + 5^3", 4, 5)
    assert _same(5 * R[1] / lg2, "4*5 + 4*5^2 + 4*5^3 + 5^4 + 2*5^5", 6, 5)
    s = kb.supersingular_bound(E, 5, cross_check_plus=True)
    b = sz.supersingular_bound(E, 5, s["Lminus_mod_p"][2], 2, G1483, Lplus_T2_mod_p=s["Lplus_mod_p"][2])
    assert b["b_p"] == 0 and b["sha_p_trivial"] and b["direction_check"]["agree"]
    assert sz.saturation_certificate(E, G1483, 5)["ok"] and sz.rhobar_surjective_certificate(E, 5)["ok"]


def test_supersingular_bound_rejects():
    with pytest.raises(ValueError):
        sz.supersingular_bound(pari.ellinit(C389), 5, 1, 2, G389)
    assert not sz.supersingular_bound(pari.ellinit(C389), 107, 0, 2, G389).get("sha_p_trivial")


# ---------------------------------------------------------------------------------------------
# split multiplicative p
# ---------------------------------------------------------------------------------------------
def test_split_multiplicative_stein_wuthrich_446d1_p223():
    """SW section 12.1: L-invariant 179*223 + 85*223^2 + 30*223^3 + O(223^4), regulator 153*223^2 + 125*223^3 +
    124*223^4 + O(223^5), [T^3]L_p = 139 + O(223) (ours is twice theirs), Sha(223) = 0."""
    E = pari.ellinit(E0)
    m = kb.multiplicative_bound(E, 223)
    assert m["ord_upper_bound"] == 3 and m["coeffs"][3] % 223 == 2 * 139 % 223
    p, n = 223, 6
    Ep = pari.ellinit(E0, pari(f"O({p}^{n})"))
    u2, u, q = pari("(E)->E.tate")(Ep)[:3]
    assert _same(pari.log(q) / pari.valuation(q, p), "179*223 + 85*223^2 + 30*223^3", 4, p)
    H = pari.ellpadicheightmatrix(E, p, n, pari(G0))
    R = pari.matdet(H[0] + H[1] * (-pari.ellpadics2(E, p, n) + 1 / pari.log(q) / u2))
    assert _same(R, "153*223^2 + 125*223^3 + 124*223^4", 5, p)
    b = sz.split_multiplicative_bound(E, p, 0, 2, G0)
    assert b["b_p"] == 0 and b["vp_L_invariant"] == 1 and b["vp_Reg_gamma"] == 0 and b["sha_p_trivial"]
    assert b["b_p_with_L_invariant_as_printed_in_algorithm"] == -1      # the literal Algorithm 11.1 is inconsistent


def test_split_multiplicative_389a1_p389():
    E = pari.ellinit(C389)
    m = kb.multiplicative_bound(E, 389)
    assert m["ord_upper_bound"] == 3 and m["coeffs"][3] % 389 == 321
    b = sz.split_multiplicative_bound(E, 389, 0, 2, G389)
    assert b["b_p"] == 0 and b["vp_L_invariant"] == 1 and b["vp_Reg_gamma"] == 0 and b["consistent"]
