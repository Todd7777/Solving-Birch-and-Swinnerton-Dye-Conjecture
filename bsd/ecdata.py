"""Access to John Cremona's database of elliptic curves over Q (ecdata).

The files are fetched from https://github.com/JohnCremona/ecdata and cached
locally.  They are used (a) as a source of curves, (b) as *independent*
reference values (ranks, generators, regulators, analytic Sha) against which
our own PARI computations are compared.  Cremona's values were computed with
different software (eclib/mwrank), so agreement is a genuine cross-check.

File formats (one curve per line):
  allcurves: N iso num [a1,a2,a3,a4,a6] rank torsion_order
  allbsd:    N iso num ainvs rank torsion_order tamagawa_product omega L^{(r)}(1)/r! regulator sha_an
  allgens:   N iso num ainvs rank [torsion_structure] P1 P2 ... (projective [x:y:z])
"""
from __future__ import annotations

import os
import urllib.request
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Iterable

ECDATA_RAW = "https://raw.githubusercontent.com/JohnCremona/ecdata/master"
DEFAULT_CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "ecdata")


def _range_tag(lo: int) -> str:
    hi = lo + 9999
    return f"{lo:05d}-{hi:05d}"


def fetch(kind: str, lo: int, cache_dir: str = DEFAULT_CACHE) -> str:
    """Download (if needed) and return the local path of e.g. allbsd.00000-09999."""
    assert kind in ("allcurves", "allbsd", "allgens")
    os.makedirs(cache_dir, exist_ok=True)
    name = f"{kind}.{_range_tag(lo)}"
    path = os.path.join(cache_dir, name)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        url = f"{ECDATA_RAW}/{kind}/{name}"
        urllib.request.urlretrieve(url, path)
    return path


@dataclass
class CremonaCurve:
    label: str
    conductor: int
    iso: str
    number: int
    ainvs: list[int]
    rank: int
    torsion_order: int | None = None
    torsion_structure: list[int] = field(default_factory=list)
    tamagawa_product: int | None = None
    omega: float | None = None
    L_value: float | None = None          # L^{(r)}(E,1)/r!
    regulator: float | None = None
    sha_an: float | None = None
    gens: list[list[Fraction]] = field(default_factory=list)   # affine, infinite order only
    torsion_gens: list[list[Fraction]] = field(default_factory=list)


def _parse_ainvs(s: str) -> list[int]:
    return [int(t) for t in s.strip("[]").split(",")]


def _parse_point(s: str) -> list[Fraction]:
    x, y, z = (int(t) for t in s.strip("[]").split(":"))
    return [Fraction(x, z), Fraction(y, z)]


def load_range(lo: int, hi: int | None = None, cache_dir: str = DEFAULT_CACHE,
               with_bsd: bool = True, with_gens: bool = True) -> list[CremonaCurve]:
    """Load all curves with lo <= conductor <= hi (hi defaults to lo's 10000-block end)."""
    block = (lo // 10000) * 10000
    if hi is None:
        hi = block + 9999
    curves: dict[str, CremonaCurve] = {}
    order: list[str] = []
    while block <= hi:
        with open(fetch("allcurves", block, cache_dir)) as f:
            for line in f:
                parts = line.split()
                if len(parts) < 6:
                    continue
                N = int(parts[0])
                if N < lo or N > hi:
                    continue
                label = f"{parts[0]}{parts[1]}{parts[2]}"
                c = CremonaCurve(label=label, conductor=N, iso=parts[1], number=int(parts[2]),
                                 ainvs=_parse_ainvs(parts[3]), rank=int(parts[4]),
                                 torsion_order=int(parts[5]))
                curves[label] = c
                order.append(label)
        if with_bsd:
            with open(fetch("allbsd", block, cache_dir)) as f:
                for line in f:
                    parts = line.split()
                    if len(parts) < 11:
                        continue
                    label = f"{parts[0]}{parts[1]}{parts[2]}"
                    c = curves.get(label)
                    if c is None:
                        continue
                    c.tamagawa_product = int(parts[6])
                    c.omega = float(parts[7])
                    c.L_value = float(parts[8])
                    c.regulator = float(parts[9])
                    c.sha_an = float(parts[10])
        if with_gens:
            with open(fetch("allgens", block, cache_dir)) as f:
                for line in f:
                    parts = line.split()
                    if len(parts) < 6:
                        continue
                    label = f"{parts[0]}{parts[1]}{parts[2]}"
                    c = curves.get(label)
                    if c is None:
                        continue
                    r = int(parts[4])
                    tstruct = parts[5].strip("[]")
                    c.torsion_structure = [int(t) for t in tstruct.split(",")] if tstruct else []
                    pts = [_parse_point(p) for p in parts[6:]]
                    c.gens = pts[:r]
                    c.torsion_gens = pts[r:]
        block += 10000
    return [curves[l] for l in order]


def by_label(label: str, cache_dir: str = DEFAULT_CACHE) -> CremonaCurve:
    """Look a single curve up by Cremona label, e.g. '389a1'."""
    import re
    m = re.fullmatch(r"(\d+)([a-z]+)(\d+)", label)
    if not m:
        raise ValueError(f"not a Cremona label: {label}")
    N = int(m.group(1))
    for c in load_range(N, N, cache_dir):
        if c.label == label:
            return c
    raise KeyError(label)


def labels(curves: Iterable[CremonaCurve]) -> list[str]:
    return [c.label for c in curves]


def fetch_all(lo: int = 0, hi: int = 499999, cache_dir: str = DEFAULT_CACHE, kinds=("allcurves", "allbsd", "allgens")) -> list[str]:
    """Download every ecdata block covering conductors lo..hi (Cremona's database currently reaches 499999).

    Run this on a machine with internet access (e.g. a cluster login node) before running sweeps on
    compute nodes that cannot reach GitHub.
    """
    paths = []
    block = (lo // 10000) * 10000
    while block <= hi:
        for kind in kinds:
            paths.append(fetch(kind, block, cache_dir))
        block += 10000
    return paths
