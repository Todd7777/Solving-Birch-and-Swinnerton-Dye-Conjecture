"""The analytic side: L(E,s) near s = 1.

* analytic_rank      : numerical order of vanishing (PARI ellanalyticrank).
* leading_term       : L^{(r)}(E,1)/r!  (PARI ellL1 returns the raw derivative).
* exact_L1_ratio     : L(E,1)/Ω⁺ as an *exact rational number*, from modular
                       symbols (PARI msfromell).  This is how one *certifies*
                       L(E,1) = 0 — a numerical value of 1e-20 proves nothing.
* certify_order      : what is rigorously known about ord_{s=1} L(E,s) given
                       the exact ratio, the root number and the numerics.

Rigour caveat, stated once here and repeated in the docs: a non-zero
derivative L^{(k)}(E,1) ≠ 0 is "certified" only in the sense that PARI's
value is many orders of magnitude above its error; PARI does not use interval
arithmetic.  Vanishing of L^{(k)}(E,1) for k ≥ 1 can never be certified
numerically; for k = 0 it can (modular symbols), and for odd/even parity the
functional equation gives ord ≥ 1 for root number -1.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import factorial

from .pari_session import get_pari


@dataclass
class AnalyticData:
    analytic_rank: int
    leading_term: float               # L^{(r)}(E,1)/r!
    raw_derivative: float             # L^{(r)}(E,1)
    L1_over_omega_plus_exact: Fraction | None
    L1_over_omega_plus_numeric: float
    rigorous_order_lower: int         # proven lower bound for ord_{s=1} L(E,s)
    rigorous_order_exact: bool        # True if ord is rigorously pinned (modulo the numerical non-vanishing caveat)
    certification: str
    L1_status: str = ""               # "zero" / "nonzero" / "" (unknown exactly)
    L1_source: str = ""               # how the exact status of L(E,1) was obtained (or why it was not)


def analytic_rank(E, eps=None):
    pari = get_pari()
    res = pari.ellanalyticrank(E) if eps is None else pari.ellanalyticrank(E, eps)
    return int(res[0]), float(res[1])


def derivative(E, r: int) -> float:
    pari = get_pari()
    return float(pari.ellL1(E, r))


def exact_L1_ratio(E) -> Fraction:
    """|L(E,1)/Ω⁺| as an exact rational, via modular symbols."""
    pari = get_pari()
    Mx = pari.msfromell(E, 1)
    M, xp = Mx[0], Mx[1]
    val = pari.mseval(M, xp, pari("[0,oo]"))
    fr = Fraction(str(val))
    return abs(fr)


def certify_order(r_num: int, root_number: int, exact_ratio: Fraction | None = None,
                  l1_zero: bool | None = None, l1_source: str = "modular symbols",
                  points_rank: int = 0) -> tuple[int, bool, str]:
    """Return (rigorous lower bound for ord, whether ord is pinned, explanation).

    The exact knowledge about L(E,1) can be given either as `exact_ratio` (an exact rational
    L(E,1)/Ω⁺) or as `l1_zero` (True/False: L(E,1) = 0 known exactly), with `l1_source` saying how.

    `points_rank` is a proven lower bound for rank E(Q) (independent rational points).  By
    Gross–Zagier + Kolyvagin + Bump–Friedberg–Hoffstein/Murty–Murty, ord_{s=1} L(E,s) ≤ 1 implies
    rank E(Q) = ord; contrapositively rank ≥ k implies ord ≥ min(k, 2).  With the parity of ord fixed
    by the root number this gives: one point and root number +1 ⇒ ord ≥ 2; two points and root
    number −1 ⇒ ord ≥ 3.  (This is why analytic rank 3 is certifiable for curves such as 5077a1;
    Buhler–Gross–Zagier 1985 did it before Kolyvagin's theorem by showing a Heegner point is torsion.)
    Nothing known gives ord ≥ 4 for any curve.
    """
    if l1_zero is None and exact_ratio is not None:
        l1_zero = (exact_ratio == 0)
    parity = 0 if root_number == 1 else 1
    if l1_zero is False:
        return 0, True, f"L(E,1) ≠ 0 certified exactly ({l1_source}): ord = 0 rigorously."
    lower = 0
    notes = []
    if root_number == -1:
        lower = 1
        notes.append("root number −1 forces ord odd, so ord ≥ 1")
    elif l1_zero is True:
        lower = 2
        notes.append(f"L(E,1) = 0 exactly ({l1_source}); with root number +1, ord is even, so ord ≥ 2")
    if points_rank >= 1:
        gzk = min(points_rank, 2)
        if gzk % 2 != parity:
            gzk += 1
        if gzk > lower:
            lower = gzk
            notes.append(f"{points_rank} independent rational point(s) ⇒ rank ≥ {min(points_rank, 2)} ⇒ ord ≥ {min(points_rank, 2)} "
                         f"(Gross–Zagier–Kolyvagin: ord ≤ 1 forces rank = ord), and ord ≡ {parity} (mod 2) from the root number, "
                         f"so ord ≥ {gzk}")
    if r_num == lower:
        prefix = ("; ".join(notes) + "; ") if notes else ""
        return lower, True, prefix + f"L^({r_num})(E,1) ≠ 0 numerically, so ord = {r_num} (numerical non-vanishing)."
    # r_num > lower: vanishing of higher derivatives cannot be certified from numerics alone
    return lower, False, (("; ".join(notes) + "; ") if notes else "") + (
        f"numerically ord = {r_num}, but vanishing of the derivatives of order {lower}..{r_num-1} is not certified here — "
        f"rigorous statement: ord ≥ {lower} with ord ≡ {parity} (mod 2).")


MODSYM_MAX_CONDUCTOR = 300_000   # modular symbols of level N cost ~N memory/time


def analytic_data(E, root_number: int, want_exact: bool = True,
                  modsym_max_conductor: int = MODSYM_MAX_CONDUCTOR, precomputed: dict | None = None) -> AnalyticData:
    """Analytic invariants at s = 1.

    `precomputed` (optional) supplies exact knowledge of L(E,1) obtained elsewhere, e.g. once per
    isogeny class by the sweep: {"ratio": Fraction | None, "zero": bool | None, "source": str}.
    For root number −1 no modular-symbol computation is done: the functional equation forces L(E,1) = 0.
    """
    pari = get_pari()
    N = int(pari.ellglobalred(E)[0])
    r, raw = analytic_rank(E)
    lead = raw / factorial(r)
    om_plus = float(pari("(E)->E.omega[1]")(E))
    L1 = derivative(E, 0)
    exact = None
    zero = None
    source = ""
    if root_number == -1:
        zero, source = True, "forced by the functional equation (root number −1)"
    elif precomputed is not None:
        exact, zero, source = precomputed.get("ratio"), precomputed.get("zero"), precomputed.get("source", "")
    elif not want_exact:
        source = "not computed"
    elif N > modsym_max_conductor:
        source = f"not computed: conductor {N} exceeds the modular-symbol limit {modsym_max_conductor}"
    else:
        try:
            exact = exact_L1_ratio(E)
            zero, source = (exact == 0), "modular symbols (PARI msfromell)"
        except Exception as exc:  # e.g. PARI stack overflow at highly composite levels
            source = f"modular symbols failed ({type(exc).__name__}: {str(exc)[:60]})"
    if exact is None and zero is None and source and not source.startswith("not computed"):
        source += "; ord = 0 below rests on numerical non-vanishing only" if abs(L1) > 1e-10 else ""
    lower, is_exact, why = certify_order(r, root_number, exact, l1_zero=zero, l1_source=source or "modular symbols")
    if zero is None and source:
        why += f" [exact status of L(E,1): {source}]"
    status = "" if zero is None else ("zero" if zero else "nonzero")
    return AnalyticData(analytic_rank=r, leading_term=lead, raw_derivative=raw,
                        L1_over_omega_plus_exact=exact, L1_over_omega_plus_numeric=L1 / om_plus,
                        rigorous_order_lower=lower, rigorous_order_exact=is_exact, certification=why,
                        L1_status=status, L1_source=source)
