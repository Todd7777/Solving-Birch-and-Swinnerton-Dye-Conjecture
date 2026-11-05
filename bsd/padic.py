"""p-adic BSD (Mazur–Tate–Teitelbaum, Bernardi–Perrin-Riou) numerics via PARI.

ellpadicbsd(E, p, n) returns [r_p, L_p] where r_p is the order of vanishing
of the p-adic L-function at the trivial character (an upper bound, with the
r_p-th derivative certified non-zero mod p^N) and L_p is conjecturally
R_p · #Ш, with R_p the p-adic regulator (PARI's normalisation).  For split
multiplicative p the p-adic L-function has an extra "exceptional" zero
(Mazur–Tate–Teitelbaum; proved by Greenberg–Stevens 1993), so r_p = r + 1.

For supersingular p both L_p and R_p are pairs (the ± / (ω, Fω) components).
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

from .pari_session import get_pari, time_limit
from .descent import points_to_gp


@dataclass
class PadicBSD:
    p: int
    reduction: str
    r_p: int
    L_p: str
    R_p: str | None
    sha_p: str | None       # L_p / R_p, conjecturally #Ш (p-adically)
    note: str = ""


def padic_bsd(ainvs, p: int, n: int = 6, gens=None, timeout: float | None = 300.0) -> PadicBSD:
    pari = get_pari()
    E = pari.ellinit(list(ainvs))
    N = int(pari.ellglobalred(E)[0])
    a_p = int(pari.ellap(E, p))
    if N % p == 0:
        lr = pari.elllocalred(E, p)
        red = ("split multiplicative" if int(lr[0]) == 1 and a_p == 1 else
               "nonsplit multiplicative" if int(lr[0]) == 1 else "additive")
    else:
        red = "good ordinary" if a_p % p != 0 else "good supersingular"
    if red == "additive":
        return PadicBSD(p, red, -1, "", None, None, "additive reduction: not handled")
    with time_limit(timeout):
        res = pari.ellpadicbsd(E, p, n)
        r_p = int(res[0])
        Lp = res[1]
        Rp = None
        sha = None
        if gens:
            V = pari(points_to_gp(gens))
            Rp = pari.ellpadicregulator(E, p, n, V)
            try:
                if pari.type(Lp) == "t_VEC":
                    sha = pari("(a,b)->[a[1]/b[1], a[2]/b[2]]")(Lp, Rp)
                else:
                    sha = Lp / Rp
            except Exception as exc:  # pragma: no cover
                sha = None
        elif r_p == 0:
            sha = Lp
    return PadicBSD(p, red, r_p, str(Lp), None if Rp is None else str(Rp), None if sha is None else str(sha))


def padic_bsd_dict(*a, **k) -> dict:
    return asdict(padic_bsd(*a, **k))
