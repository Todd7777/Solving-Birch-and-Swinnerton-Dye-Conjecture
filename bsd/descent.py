"""The algebraic side: Mordell–Weil group via 2-descent, saturation, regulator.

PARI's ellrank performs a 2-descent (with the 2-part of the Cassels–Tate
pairing).  It returns [r1, r2, s, points] with r1 ≤ rank ≤ r2, where
r2 = dim Sel_2 − dim E(Q)[2] − s and s is the (even) 2-rank of Ш[2]/2Ш[4]
detected by the Cassels–Tate pairing.  r1 = r2 is only possible when Ш has
no 4-torsion.  The points returned are independent but NOT saturated: we
saturate at all primes below `saturation_bound` with ellsaturation.

Honesty note: proving that the saturation is complete at *all* primes needs
a provable lower bound for the canonical height on E(Q) (Cremona–Siksek),
which PARI does not expose.  We therefore record the bound used, and in the
systematic sweeps we compare the regulator with Cremona's independently
computed value.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

from .pari_session import get_pari, time_limit


@dataclass
class DescentData:
    rank_lower: int
    rank_upper: int
    sha2_mod_4_rank: int                       # PARI's s
    points_found: list[list[str]]
    gens: list[list[str]] = field(default_factory=list)   # saturated basis (as strings "a/b")
    gens_source: str = ""
    saturation_bound: int = 0
    regulator: float = 1.0
    heights: list[float] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def _pt_to_gp(P) -> str:
    return "[" + ",".join(str(Fraction(c)) for c in P) + "]"


def points_to_gp(points) -> str:
    return "[" + ",".join(_pt_to_gp(P) for P in points) + "]"


def _gp_points_to_py(V) -> list[list[str]]:
    out = []
    for P in V:
        out.append([str(P[0]), str(P[1])])
    return out


def two_descent(E, effort: int = 0, points=None, timeout: float | None = 120.0):
    """Run PARI's ellrank; `points` (affine rational points) are used if given."""
    pari = get_pari()
    with time_limit(timeout):
        if points:
            V = pari(points_to_gp(points))
            res = pari.ellrank(E, effort, V)
        else:
            res = pari.ellrank(E, effort)
    return int(res[0]), int(res[1]), int(res[2]), res[3]


def saturate(E, V, bound: int = 1000, timeout: float | None = 300.0):
    pari = get_pari()
    if len(V) == 0:
        return V
    with time_limit(timeout):
        return pari.ellsaturation(E, V, bound)


def regulator(E, V) -> tuple[float, list[float]]:
    pari = get_pari()
    if len(V) == 0:
        return 1.0, []
    H = pari.ellheightmatrix(E, V)
    n = len(V)
    heights = [float(H[i, i]) for i in range(n)]
    return float(pari.matdet(H)), heights


def is_torsion(E, P) -> bool:
    pari = get_pari()
    return int(pari.ellorder(E, P)) != 0


def mordell_weil(E, known_gens=None, effort: int = 0, saturation_bound: int = 1000,
                 descent_timeout: float | None = 120.0, run_descent: bool = True) -> DescentData:
    """Compute rank bounds and a saturated basis of E(Q)/torsion.

    If `known_gens` (e.g. from Cremona's tables) are supplied they are fed to
    the 2-descent (which then verifies the rank bounds independently) and
    re-saturated by PARI.
    """
    pari = get_pari()
    notes = []
    r1 = r2 = s = None
    V = None
    if run_descent:
        try:
            r1, r2, s, V = two_descent(E, effort=effort, points=known_gens, timeout=descent_timeout)
        except Exception as exc:  # timeout or PARI error
            notes.append(f"2-descent failed/timed out: {type(exc).__name__}: {str(exc)[:80]}")
    if V is None or (known_gens and len(V) < len(known_gens)):
        if known_gens:
            V = pari(points_to_gp(known_gens))
            notes.append("using supplied generators without a fresh 2-descent")
        elif V is None:
            V = pari("[]")
    # discard torsion points, if any slipped in
    Vlist = [P for P in V if not is_torsion(E, P)]
    V = pari("[" + ",".join(str(P) for P in Vlist) + "]") if Vlist else pari("[]")
    points_found = _gp_points_to_py(V)
    gens_source = "ellrank" if run_descent and not known_gens else ("supplied+ellrank" if run_descent else "supplied")
    W = V
    if len(V) > 0:
        try:
            W = saturate(E, V, saturation_bound)
        except Exception as exc:
            notes.append(f"saturation failed: {type(exc).__name__}")
            W = V
    reg, hts = regulator(E, W)
    if r1 is None:
        r1 = len(W)
        r2 = len(W)
        s = -1
        notes.append("rank bounds not independently established (no descent)")
    return DescentData(rank_lower=r1, rank_upper=r2, sha2_mod_4_rank=s, points_found=points_found,
                       gens=_gp_points_to_py(W), gens_source=gens_source, saturation_bound=saturation_bound,
                       regulator=reg, heights=hts, notes=notes)
