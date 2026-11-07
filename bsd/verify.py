"""Verify the full BSD formula for one curve and assemble the record."""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field, asdict
from fractions import Fraction

from .pari_session import get_pari, time_limit
from .invariants import compute_invariants, CurveInvariants
from .analytic import analytic_data, AnalyticData, certify_order
from .descent import mordell_weil, DescentData
from .theorems import p_part_ledger, rank_and_sha_status, sha_p_finite_via_kato, SMALL_PRIMES
from .padic import padic_bsd


@dataclass
class BSDRecord:
    label: str
    invariants: dict
    analytic: dict
    descent: dict
    sha_an: float | None
    sha_an_int: int | None
    sha_is_square: bool | None
    sha_rel_err: float | None
    rank_status: str
    sha_status: str
    p_part: list[dict]
    padic: list[dict]
    sage: dict | None
    cremona: dict | None
    timings: dict
    errors: list[str]

    def to_dict(self) -> dict:
        d = asdict(self)
        return _jsonable(d)


def _jsonable(x):
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, Fraction):
        return str(x)
    return x


def _is_square(n: int) -> bool:
    if n < 0:
        return False
    s = math.isqrt(n)
    return s * s == n


def verify_curve(ainvs, label: str = "", known_gens=None, run_descent: bool = True, effort: int = 0,
                 saturation_bound: int = 1000, want_exact: bool = True, descent_timeout: float | None = 120.0,
                 analytic_timeout: float | None = 1800.0, cremona=None, ledger: bool = True,
                 padic_bound: int = 0, padic_prec: int = 5, padic_timeout: float | None = 600.0,
                 use_sage: bool = False, sage_timeout: float | None = 600.0,
                 l1_precomputed: dict | None = None) -> BSDRecord:
    pari = get_pari()
    t = {}
    errors: list[str] = []
    t0 = time.time()
    inv = compute_invariants(ainvs)
    t["invariants"] = time.time() - t0
    E = pari.ellinit(inv.ainvs)
    if not label:
        label = f"N{inv.conductor}:{inv.ainvs}"
    # analytic side
    t0 = time.time()
    try:
        with time_limit(analytic_timeout):
            an = analytic_data(E, inv.root_number, want_exact=want_exact, precomputed=l1_precomputed)
    except Exception as exc:
        errors.append(f"analytic side failed: {type(exc).__name__}: {str(exc)[:100]}")
        an = AnalyticData(-1, float("nan"), float("nan"), None, float("nan"), 0, False, "failed")
    t["analytic"] = time.time() - t0
    # algebraic side
    t0 = time.time()
    try:
        des = mordell_weil(E, known_gens=known_gens, effort=effort, saturation_bound=saturation_bound,
                           descent_timeout=descent_timeout, run_descent=run_descent)
    except Exception as exc:
        errors.append(f"algebraic side failed: {type(exc).__name__}: {str(exc)[:100]}")
        des = DescentData(0, 99, -1, [], [], "failed", 0, float("nan"), [], [str(exc)[:100]])
    t["descent"] = time.time() - t0
    # sharpen the analytic certificate with the independent points found (Gross–Zagier–Kolyvagin + parity)
    if an.analytic_rank >= 0 and des.rank_lower >= 1:
        l1_zero = {"zero": True, "nonzero": False}.get(an.L1_status)
        lower2, exact2, why2 = certify_order(an.analytic_rank, inv.root_number, an.L1_over_omega_plus_exact,
                                             l1_zero=l1_zero, l1_source=an.L1_source or "modular symbols",
                                             points_rank=des.rank_lower)
        if lower2 > an.rigorous_order_lower or (exact2 and not an.rigorous_order_exact):
            an.rigorous_order_lower, an.rigorous_order_exact, an.certification = lower2, exact2, why2
    # BSD quotient
    sha_an = sha_int = rel = None
    square = None
    r = an.analytic_rank
    ngens = len(des.gens)
    if r >= 0 and ngens == r and not math.isnan(des.regulator):
        sha_an = an.leading_term / (inv.bsd_constant * des.regulator)
        sha_int = int(round(sha_an))
        rel = abs(sha_an - sha_int) / max(sha_int, 1)
        square = _is_square(sha_int) and sha_int >= 1
        if rel > 1e-6:
            errors.append(f"Ш_an = {sha_an!r} is not close to an integer (rel. err {rel:.2e})")
        elif not square:
            errors.append(f"Ш_an = {sha_int} is not a perfect square — contradicts Cassels' theorem if Ш is finite")
    elif r >= 0 and ngens != r:
        errors.append(f"found {ngens} independent points but analytic rank is {r}: Ш_an not computed")
    if r >= 0 and not (des.rank_lower <= r <= des.rank_upper):
        errors.append(f"analytic rank {r} outside 2-descent bounds [{des.rank_lower},{des.rank_upper}]")
    rank_status, sha_status = rank_and_sha_status(r, an.rigorous_order_lower, an.rigorous_order_exact,
                                                  des.rank_lower, des.rank_upper, ngens,
                                                  l1_exact_nonzero=(an.L1_status == "nonzero"))
    pp = p_part_ledger(inv, r, sha_int) if (ledger and r >= 0) else []
    # p-adic BSD data and, for rank ≥ 2, Kato-based finiteness of Ш[p^∞]
    padic = []
    if padic_bound > 0 and r >= 0:
        gens = [[Fraction(x), Fraction(y)] for x, y in des.gens] if des.gens else None
        for p in [q for q in SMALL_PRIMES if 3 <= q <= padic_bound and q not in inv.bad_primes]:
            try:
                d = asdict(padic_bsd(inv.ainvs, p, n=padic_prec, gens=gens, timeout=padic_timeout))
            except Exception as exc:
                d = {"p": p, "reduction": "?", "r_p": -1, "L_p": "", "R_p": None, "sha_p": None,
                     "note": f"failed: {type(exc).__name__}: {str(exc)[:80]}"}
            d["kato"] = sha_p_finite_via_kato(inv, des.rank_lower, p, d["r_p"])
            padic.append(d)
    # optional Sage stage: provable saturation + prove_BSD (rank ≤ 1)
    sage_info = None
    if use_sage:
        sage_info = {}
        try:
            from .sage_bridge import prove_bsd, provable_saturation
            if des.gens:
                sat = provable_saturation(inv.ainvs, des.gens, timeout=sage_timeout)
                sage_info["saturation"] = sat
                if sat.get("index") not in (None, 1):
                    errors.append(f"Sage found the PARI-saturated basis to have index {sat['index']}: regulator corrected")
                    des.regulator = sat["regulator"]
                    des.gens = sat["gens"]
                    des.notes.append("basis replaced by Sage's provably saturated basis")
                    if r >= 0 and len(des.gens) == r:
                        sha_an = an.leading_term / (inv.bsd_constant * des.regulator)
                        sha_int = int(round(sha_an))
                        rel = abs(sha_an - sha_int) / max(sha_int, 1)
                        square = _is_square(sha_int) and sha_int >= 1
            if 0 <= r <= 1:
                sage_info["prove_bsd"] = prove_bsd(inv.ainvs, timeout=sage_timeout)
        except Exception as exc:  # noqa: BLE001
            sage_info["error"] = f"{type(exc).__name__}: {str(exc)[:160]}"
    # comparison with Cremona's tables
    crem = None
    if cremona is not None:
        crem = compare_with_cremona(inv, an, des, sha_an, cremona)
        if crem.get("mismatches"):
            errors.append("mismatch with Cremona's table: " + ", ".join(crem["mismatches"]))
    return BSDRecord(label=label, invariants=inv.to_dict(), analytic=_jsonable(asdict(an)), descent=asdict(des),
                     sha_an=sha_an, sha_an_int=sha_int, sha_is_square=square, sha_rel_err=rel,
                     rank_status=rank_status, sha_status=sha_status, p_part=pp, padic=padic, sage=sage_info, cremona=crem, timings=t,
                     errors=errors)


def compare_with_cremona(inv: CurveInvariants, an: AnalyticData, des: DescentData, sha_an, c) -> dict:
    """Compare our invariants with Cremona's allbsd values (independent software)."""
    def close(a, b, tol=1e-8):
        if a is None or b is None:
            return None
        return abs(a - b) <= tol * max(1.0, abs(b))
    out = {"label": c.label, "rank": c.rank, "sha_an": c.sha_an, "regulator": c.regulator,
           "omega": c.omega, "tamagawa_product": c.tamagawa_product, "L_value": c.L_value}
    mism = []
    if c.rank != an.analytic_rank:
        mism.append(f"rank {c.rank} vs analytic {an.analytic_rank}")
    if c.torsion_order is not None and c.torsion_order != inv.torsion_order:
        mism.append("torsion")
    if c.tamagawa_product is not None and c.tamagawa_product != inv.tamagawa_product:
        mism.append(f"tamagawa {c.tamagawa_product} vs {inv.tamagawa_product}")
    if close(inv.omega, c.omega, 1e-9) is False:
        mism.append(f"omega {inv.omega} vs {c.omega}")
    if close(an.leading_term, c.L_value, 1e-9) is False:
        mism.append(f"L-value {an.leading_term} vs {c.L_value}")
    if close(des.regulator, c.regulator, 1e-8) is False:
        mism.append(f"regulator {des.regulator} vs {c.regulator}")
    if sha_an is not None and close(sha_an, c.sha_an, 1e-6) is False:
        mism.append(f"sha {sha_an} vs {c.sha_an}")
    out["mismatches"] = mism
    return out
