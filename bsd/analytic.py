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


def certify_order(r_num: int, root_number: int, exact_ratio: Fraction | None) -> tuple[int, bool, str]:
    """Return (rigorous lower bound for ord, whether ord is rigorously exact, explanation)."""
    parity = 0 if root_number == 1 else 1
    if exact_ratio is not None and exact_ratio != 0:
        return 0, True, "L(E,1)/Ω⁺ is a non-zero rational (modular symbols): ord = 0 rigorously."
    lower = 0
    notes = []
    if root_number == -1:
        lower = 1
        notes.append("root number -1 forces ord odd, so ord ≥ 1")
    if exact_ratio == 0:
        notes.append("L(E,1) = 0 exactly (modular symbols)")
        if root_number == 1:
            lower = 2
            notes.append("with root number +1, ord is even, so ord ≥ 2")
    if r_num == lower:
        return lower, True, "; ".join(notes) + f"; L^{{({r_num})}}(E,1) ≠ 0 numerically, so ord = {r_num} (numerical non-vanishing)."
    # r_num > lower: vanishing of higher derivatives cannot be certified numerically
    return lower, False, ("; ".join(notes) + f"; numerically ord = {r_num}, but vanishing of derivatives of order "
                          f"{lower}..{r_num-1} cannot be certified numerically — rigorous statement is only ord ≥ {lower} "
                          f"with ord ≡ {parity} (mod 2).")


MODSYM_MAX_CONDUCTOR = 300_000   # modular symbols of level N cost ~N memory/time


def analytic_data(E, root_number: int, want_exact: bool = True,
                  modsym_max_conductor: int = MODSYM_MAX_CONDUCTOR) -> AnalyticData:
    pari = get_pari()
    N = int(pari.ellglobalred(E)[0])
    if N > modsym_max_conductor:
        want_exact = False
    r, raw = analytic_rank(E)
    lead = raw / factorial(r)
    om_plus = float(pari("(E)->E.omega[1]")(E))
    L1 = derivative(E, 0)
    exact = None
    if want_exact:
        try:
            exact = exact_L1_ratio(E)
        except Exception:  # pragma: no cover
            exact = None
    lower, is_exact, why = certify_order(r, root_number, exact)
    if exact is None and want_exact is False and N > modsym_max_conductor:
        why += f" [exact L(E,1)/Ω⁺ not computed: conductor {N} exceeds the modular-symbol limit {modsym_max_conductor}]"
    return AnalyticData(analytic_rank=r, leading_term=lead, raw_derivative=raw,
                        L1_over_omega_plus_exact=exact, L1_over_omega_plus_numeric=L1 / om_plus,
                        rigorous_order_lower=lower, rigorous_order_exact=is_exact, certification=why)
