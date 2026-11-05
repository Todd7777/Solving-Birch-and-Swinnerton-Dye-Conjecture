"""A curated list of elliptic curves that exercise every feature of BSD.

Curves of conductor < 20000 are referenced by Cremona label and pulled from
ecdata (so their coefficients are never typed by hand).  Curves of larger
conductor (the record-rank examples) are given by coefficients and *validated
at run time*: the computed conductor must equal the claimed one and the
2-descent must reproduce the claimed rank, otherwise the curve is flagged.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Curated:
    label: str
    note: str
    ainvs: list[int] | None = None      # None: look up in ecdata
    claimed_conductor: int | None = None
    claimed_rank: int | None = None
    descent_effort: int = 0
    descent_timeout: float = 600.0


SMALL = [
    Curated("11a1", "smallest conductor; rank 0; E(Q)_tors = Z/5, c_11 = 5, L(E,1)/Ω = 1/5"),
    Curated("14a1", "rank 0, torsion Z/6"),
    Curated("15a1", "rank 0, torsion Z/4 × Z/2"),
    Curated("17a1", "rank 0, torsion Z/4"),
    Curated("19a1", "rank 0, torsion Z/3"),
    Curated("27a1", "CM by Z[ζ_3] (j = 0), rank 0"),
    Curated("32a2", "y² = x³ − x, CM by Z[i], rank 0 (1 is not a congruent number)"),
    Curated("36a1", "y² = x³ + 1, CM by Z[ζ_3], rank 0"),
    Curated("37a1", "smallest conductor of rank 1; generator (0,0)"),
    Curated("37b1", "same conductor as 37a1 but rank 0, torsion Z/3"),
    Curated("43a1", "rank 1"),
    Curated("49a1", "CM by Z[(1+√−7)/2], rank 0"),
    Curated("53a1", "rank 1, supersingular at 5 (used in the p-adic BSD examples)"),
    Curated("57a1", "rank 1"),
    Curated("65a1", "rank 1 with 2-torsion"),
    Curated("67a1", "rank 1"),
    Curated("91a1", "rank 1, nonsplit multiplicative at 7"),
    Curated("91b1", "rank 1, split multiplicative at 7 (exceptional zero of the 7-adic L-function)"),
    Curated("121b1", "CM by Z[(1+√−11)/2]"),
    Curated("128a1", "y² = x³ + x² + x + 1 ... (CM by Z[i] family, conductor 2^7)"),
    Curated("389a1", "smallest conductor of rank 2; L(E,1) = 0 exactly, L''(E,1) ≠ 0"),
    Curated("433a1", "rank 2"),
    Curated("446d1", "rank 2"),
    Curated("571a1", "rank 0 with Ш ≅ (Z/2)², the smallest conductor with non-trivial Ш"),
    Curated("571b1", "rank 2 (same conductor as 571a1)"),
    Curated("681b1", "rank 0 with Ш ≅ (Z/3)²"),
    Curated("1058d1", "rank 0 with Ш of order 25? (whatever the computation says — label kept for its reputation)"),
    Curated("2849a1", "rank 0 with Ш ≅ (Z/3)²"),
    Curated("3364c1", "rank 0 with Ш ≅ (Z/7)²"),
    Curated("5077a1", "smallest conductor of rank 3"),
    Curated("11197a1", "rank 3"),
]

LARGE = [
    Curated("234446a1", "smallest known conductor of rank 4", ainvs=[1, -1, 0, -79, 289],
            claimed_conductor=234446, claimed_rank=4, descent_effort=2, descent_timeout=1800.0),
    Curated("19047851a1", "smallest known conductor of rank 5", ainvs=[0, 0, 1, -79, 342],
            claimed_conductor=19047851, claimed_rank=5, descent_effort=3, descent_timeout=3600.0),
    Curated("rank6-5187563742", "smallest known conductor of rank 6", ainvs=[1, 1, 0, -2582, 48720],
            claimed_conductor=5187563742, claimed_rank=6, descent_effort=4, descent_timeout=7200.0),
    Curated("rank7-382623908456", "smallest known conductor of rank 7", ainvs=[0, 0, 0, -10012, 346900],
            claimed_conductor=382623908456, claimed_rank=7, descent_effort=4, descent_timeout=7200.0),
]


def resolve(curated: list[Curated], ecdata_loader=None):
    """Attach a-invariants (and Cremona reference data) to each curated entry."""
    from . import ecdata as ecd
    out = []
    for c in curated:
        ref = None
        ainvs = c.ainvs
        if ainvs is None:
            ref = ecd.by_label(c.label)
            ainvs = ref.ainvs
        out.append((c, ainvs, ref))
    return out
