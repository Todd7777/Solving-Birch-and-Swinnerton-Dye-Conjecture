"""Arithmetic invariants of an elliptic curve over Q that enter the BSD formula.

Everything here is computed by PARI/GP from a (global minimal) Weierstrass
model.  Conventions (they matter — most "counterexamples" to BSD found by
beginners are normalisation errors):

* omega_plus  : least positive real period Ω⁺ = ∫ over the identity component
                of E(R) of the Néron differential ω = dx/(2y + a1 x + a3).
* omega       : Ω_E = ∫_{E(R)} |ω| = Ω⁺ · (number of real components)
                = Ω⁺ · (2 if Δ_min > 0 else 1).  This is the period in BSD.
* tamagawa_product : ∏_{p | N} c_p over *finite* primes only.
* heights     : Néron–Tate height normalised as in Cremona's book / PARI /
                LMFDB (twice Silverman's).  With this normalisation and Ω as
                above, the BSD formula reads
                    L^{(r)}(E,1)/r! = Ω_E · Reg(E/Q) · ∏ c_p · #Ш / (#E(Q)_tors)².
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from fractions import Fraction

from .pari_session import get_pari

CM_J_INVARIANTS = {
    0: -3, 1728: -4, -3375: -7, 8000: -8, -32768: -11, 54000: -12, 287496: -16,
    -884736: -19, -12288000: -27, 16581375: -28, -884736000: -43,
    -147197952000: -67, -262537412640768000: -163,
}


def kodaira_symbol(k: int) -> str:
    """PARI's integer encoding of Kodaira symbols -> string."""
    if k == 1:
        return "I0"
    if k == 2:
        return "II"
    if k == 3:
        return "III"
    if k == 4:
        return "IV"
    if k > 4:
        return f"I{k-4}"
    if k == -1:
        return "I0*"
    if k == -2:
        return "II*"
    if k == -3:
        return "III*"
    if k == -4:
        return "IV*"
    return f"I{-k-4}*"


@dataclass
class LocalData:
    p: int
    conductor_exponent: int
    kodaira: str
    ord_disc: int          # valuation of the minimal discriminant
    tamagawa: int          # c_p = [E(Q_p) : E_0(Q_p)]
    reduction: str         # good / split multiplicative / nonsplit multiplicative / additive
    ap: int


@dataclass
class CurveInvariants:
    ainvs: list[int]                       # global minimal model
    conductor: int
    discriminant: int                      # minimal discriminant
    j_invariant: str
    cm_discriminant: int | None
    bad_primes: list[int]
    local_data: list[LocalData]
    tamagawa_product: int                  # finite primes only
    real_components: int
    omega_plus: float
    omega: float
    torsion_order: int
    torsion_structure: list[int]
    root_number: int
    isogeny_degrees: list[int]             # minimal degrees of rational isogenies from E (incl. 1)
    reducible_primes: list[int]            # primes l with E[l] reducible as Galois module
    semistable: bool
    bsd_constant: float                    # Ω · ∏c_p / #tors²  (PARI ellbsd)

    def to_dict(self) -> dict:
        d = asdict(self)
        return d


def _primes_dividing(n: int) -> list[int]:
    pari = get_pari()
    if n in (0, 1, -1):
        return []
    F = pari.factor(abs(n))
    nrows = int(pari.matsize(F)[0])
    return [int(F[i, 0]) for i in range(nrows)]


def minimal_model(ainvs) -> list[int]:
    pari = get_pari()
    E = pari.ellinit(list(ainvs))
    Em = pari.ellminimalmodel(E)
    return [int(Em[i]) for i in range(5)]


def compute_invariants(ainvs) -> CurveInvariants:
    pari = get_pari()
    amin = minimal_model(ainvs)
    E = pari.ellinit(amin)
    gr = pari.ellglobalred(E)
    N = int(gr[0])
    tam_finite = int(gr[2])
    disc = int(pari("(E)->E.disc")(E))
    j = pari("(E)->E.j")(E)
    j_int = int(j) if pari.denominator(j) == 1 else None
    cm = CM_J_INVARIANTS.get(j_int) if j_int is not None else None
    bad = _primes_dividing(N)
    local = []
    for p in bad:
        lr = pari.elllocalred(E, p)
        f, kod, cp = int(lr[0]), int(lr[1]), int(lr[3])
        ordD = int(pari.valuation(disc, p))
        ap = int(pari.ellap(E, p))
        if f == 1:
            red = "split multiplicative" if ap == 1 else "nonsplit multiplicative"
        else:
            red = "additive"
        local.append(LocalData(p=p, conductor_exponent=f, kodaira=kodaira_symbol(kod),
                               ord_disc=ordD, tamagawa=cp, reduction=red, ap=ap))
    om = pari("(E)->E.omega")(E)
    omega_plus = float(om[0])
    comps = 2 if disc > 0 else 1
    tors = pari.elltors(E)
    tors_order = int(tors[0])
    tors_struct = [int(x) for x in tors[1]]
    w = int(pari.ellrootno(E))
    iso = pari.ellisomat(E, 0, 1)
    M = iso[1]
    ncols = int(pari.matsize(M)[1])
    degs = sorted({int(M[0, j]) for j in range(ncols)})
    red_primes = sorted({q for d in degs for q in _primes_dividing(d)})
    semistable = all(ld.conductor_exponent == 1 for ld in local)
    bsd_c = float(pari.ellbsd(E))
    return CurveInvariants(
        ainvs=amin, conductor=N, discriminant=disc, j_invariant=str(j), cm_discriminant=cm,
        bad_primes=bad, local_data=local, tamagawa_product=tam_finite, real_components=comps,
        omega_plus=omega_plus, omega=omega_plus * comps, torsion_order=tors_order,
        torsion_structure=tors_struct, root_number=w, isogeny_degrees=degs,
        reducible_primes=red_primes, semistable=semistable, bsd_constant=bsd_c,
    )


def ap(ainvs, p: int) -> int:
    pari = get_pari()
    return int(pari.ellap(pari.ellinit(list(ainvs)), p))
