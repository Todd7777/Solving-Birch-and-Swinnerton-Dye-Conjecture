"""What is *proven* for a given curve — a mechanical theorem ledger.

For each curve we record (i) the status of the rank conjecture and of the
finiteness of Ш, which is a theorem when the analytic rank is ≤ 1
(Gross–Zagier + Kolyvagin, with Kato for rank 0), and open otherwise; and
(ii) for analytic rank ≤ 1, a prime-by-prime ledger of the p-part of the
BSD *formula*, obtained by mechanically checking the hypotheses of:

  * Skinner–Urban (Invent. Math. 195, 2014), Thm 2, together with Kato
    (Astérisque 295, 2004): analytic rank 0, p ≥ 3 of good ordinary reduction,
    E[p] irreducible, and some prime q ≠ p with q ‖ N at which E[p] is
    ramified.  [We take the statement from the paper; the hypotheses are
    checked exactly as listed there — see docs/02-what-is-proven.md.]
  * Jetchev–Skinner–Wan (Camb. J. Math. 5, 2017): analytic rank 1, E
    semistable, p ≥ 5 of good reduction with E[p] irreducible (also p = 3
    when E is ordinary at 3, or supersingular with a_3 = 0).
  * Castella (Camb. J. Math. 6, 2018): analytic rank 1, E semistable,
    p > 3 a prime of multiplicative reduction with E[p] irreducible.

Everything else is reported as "not covered by the theorems mechanised
here" — which is NOT the same as "unknown in the literature": many further
cases are covered by Kolyvagin's and Kato's explicit bounds (as implemented
in Sage's prove_BSD), by Wan / Sprung / Castella–Çiperiani–Skinner–Sprung
(supersingular primes), Skinner (multiplicative p, rank 0),
Castella–Grossi–Lee–Skinner (Eisenstein primes), and Rubin (CM curves).

Facts used in the hypothesis checks:
  * E[p] is reducible as a Galois module iff E admits a rational p-isogeny.
  * For q ‖ N (multiplicative reduction), E[p] is unramified at q iff
    p | ord_q(Δ_min) (Tate curve), so "ramified at q" ⇔ p ∤ ord_q(Δ_min).
  * Good ordinary at p ⇔ p ∤ N and a_p ≢ 0 (mod p).
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

from .invariants import CurveInvariants, ap as _ap

SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]


@dataclass
class PPart:
    p: int
    status: str            # "proven" | "not covered"
    theorem: str
    checks: dict


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def primes_dividing(n: int) -> list[int]:
    n = abs(n)
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def ledger_primes(inv: CurveInvariants, sha_an_int: int | None) -> list[int]:
    ps = set(SMALL_PRIMES)
    ps.update(inv.bad_primes)
    ps.update(primes_dividing(inv.tamagawa_product))
    ps.update(primes_dividing(inv.torsion_order))
    if sha_an_int:
        ps.update(primes_dividing(sha_an_int))
    return sorted(ps)


def p_part_status(inv: CurveInvariants, analytic_rank: int, p: int) -> PPart:
    good = p not in inv.bad_primes
    a_p = _ap(inv.ainvs, p) if good else None
    ordinary = bool(good and a_p % p != 0)
    supersingular_ap0 = bool(good and a_p == 0)
    irreducible = p not in inv.reducible_primes
    mult_primes = [ld for ld in inv.local_data if ld.conductor_exponent == 1]
    ramified_q = [ld.p for ld in mult_primes if ld.p != p and ld.ord_disc % p != 0]
    multiplicative_at_p = any(ld.p == p for ld in mult_primes)
    checks = dict(good=good, a_p=a_p, ordinary=ordinary, irreducible=irreducible,
                  ramified_at_some_q_exactly_dividing_N=ramified_q[:3], semistable=inv.semistable,
                  multiplicative_at_p=multiplicative_at_p)
    if analytic_rank == 0:
        if p >= 3 and good and ordinary and irreducible and ramified_q:
            return PPart(p, "proven", "Skinner–Urban 2014 Thm 2 + Kato 2004", checks)
        return PPart(p, "not covered", "", checks)
    if analytic_rank == 1:
        if inv.semistable and good and irreducible and (p >= 5 or (p == 3 and (ordinary or supersingular_ap0))):
            return PPart(p, "proven", "Jetchev–Skinner–Wan 2017", checks)
        if inv.semistable and multiplicative_at_p and p > 3 and irreducible:
            return PPart(p, "proven", "Castella 2018 (multiplicative p)", checks)
        return PPart(p, "not covered", "", checks)
    return PPart(p, "open", "analytic rank ≥ 2: no theorem applies", checks)


def p_part_ledger(inv: CurveInvariants, analytic_rank: int, sha_an_int: int | None) -> list[dict]:
    return [asdict(p_part_status(inv, analytic_rank, p)) for p in ledger_primes(inv, sha_an_int)]


def rank_and_sha_status(analytic_rank: int, rigorous_lower: int, rigorous_exact: bool,
                        rank_lower: int, rank_upper: int, n_gens: int,
                        l1_exact_nonzero: bool | None = None) -> tuple[str, str]:
    """Human-readable status of (rank conjecture, finiteness of Ш) for this curve.

    `l1_exact_nonzero` says whether L(E,1) ≠ 0 was certified exactly (modular symbols); otherwise the
    rank-0 conclusion rests on the numerical non-vanishing of L(E,1) and is labelled accordingly.
    """
    if analytic_rank == 0 and rigorous_exact:
        if l1_exact_nonzero:
            rk = ("THEOREM: rank E(Q) = 0 and Ш(E/Q) is finite (L(E,1) ≠ 0 certified exactly by modular symbols; "
                  "Kolyvagin 1988 with Gross–Zagier 1986, independently Kato 2004).")
        else:
            rk = ("THEOREM (modulo the numerical non-vanishing of L(E,1)): rank E(Q) = 0 and Ш(E/Q) is finite "
                  "(Kolyvagin 1988 with Gross–Zagier 1986, independently Kato 2004).")
        if rank_upper == 0:
            rk += " 2-descent independently gives rank 0."
        return rk, "finite (theorem)"
    if analytic_rank == 1 and rigorous_exact:
        rk = ("THEOREM (modulo the numerical non-vanishing of L'(E,1)): rank E(Q) = 1 and Ш(E/Q) finite "
              "(Gross–Zagier 1986 + Kolyvagin 1990, with Bump–Friedberg–Hoffstein / Murty–Murty 1990-91 "
              "for the auxiliary quadratic field).")
        if rank_lower == rank_upper == 1:
            rk += " 2-descent independently gives rank 1."
        return rk, "finite (theorem, given L'(E,1) ≠ 0)"
    # analytic rank >= 2 (or not rigorously pinned)
    sha_open = "OPEN: finiteness of Ш is not known for any curve of rank ≥ 2 (only Ш[p^∞] for individual primes p)"
    if rigorous_exact and rank_lower == rank_upper == analytic_rank:
        rk = (f"VERIFIED FOR THIS CURVE: rank E(Q) = {rank_lower} (2-descent upper bound, {n_gens} independent points) and "
              f"ord_{{s=1}} L(E,s) = {analytic_rank} (lower bound proven as in the analytic certificate; non-vanishing of "
              f"L^({analytic_rank})(E,1) numerical), so the rank part of BSD holds for this curve. No general theorem covers "
              f"analytic rank ≥ 2: this is a curve-by-curve verification.")
        return rk, sha_open
    consistent = rank_lower <= analytic_rank <= rank_upper
    rk = (f"OPEN for this curve: analytic rank {analytic_rank} numerically, but only ord ≥ {rigorous_lower} is proven "
          f"(no elliptic curve over Q is known to have provable analytic rank ≥ 4). 2-descent gives "
          f"{rank_lower} ≤ rank ≤ {rank_upper}" + (f" and {n_gens} independent points were found" if n_gens else "") +
          (". Consistent with the rank conjecture." if consistent else ". INCONSISTENT — investigate (likely a computational failure)."))
    return rk, sha_open


def sha_p_finite_via_kato(inv: CurveInvariants, rank_lower: int, p: int, r_p: int) -> dict:
    """Finiteness of Ш(E/Q)[p^∞] from Kato's theorem and a p-adic L-function computation.

    Ingredients (all theorems):
      * Kato (2004, Thm 17.4): for p ≥ 3 of good ordinary reduction, X_∞ = Sel_{p^∞}(E/Q_∞)^∨ is
        Λ-torsion and char_Λ(X_∞) divides the Mazur–Swinnerton-Dyer p-adic L-function L_p(E,T)
        in Λ ⊗ Q_p.  Hence ord_{T=0} char(X_∞) ≤ ord_{T=0} L_p(E,T).
      * Mazur's control theorem (good ordinary p): corank_{Z_p} Sel_{p^∞}(E/Q) = rank_{Z_p}(X_∞)_Γ,
        and by the structure theorem rank_{Z_p}(X_∞)_Γ ≤ ord_{T=0} char(X_∞).
      * PARI's ellpadicbsd certifies that the r_p-th derivative of L_p at T = 0 is non-zero,
        so ord_{T=0} L_p(E,T) ≤ r_p.
    Therefore corank Sel_{p^∞}(E/Q) ≤ r_p.  If r_p ≤ (number of independent rational points
    found) ≤ rank E(Q) ≤ corank Sel_{p^∞}(E/Q), all are equal and Ш(E/Q)[p^∞] is finite.
    (This is the mechanism of Stein–Wuthrich, Math. Comp. 82 (2013), who also obtain upper
    bounds for #Ш[p^∞] using the integral form of Kato's divisibility.)
    """
    good = p not in inv.bad_primes
    a_p = _ap(inv.ainvs, p) if good else None
    ordinary = bool(good and a_p % p != 0)
    ok = p >= 3 and good and ordinary and r_p >= 0 and r_p <= rank_lower
    return dict(p=p, good_ordinary=ordinary, r_p=r_p, rank_lower=rank_lower,
                sha_p_finite=ok,
                reason=("corank Sel_{p^∞} ≤ ord_{T=0} L_p ≤ r_p = rank: Ш[p^∞] finite (Kato 2004 + control theorem + p-adic L computation)"
                        if ok else "hypotheses not met (need p ≥ 3 good ordinary and r_p ≤ number of independent points)"))
