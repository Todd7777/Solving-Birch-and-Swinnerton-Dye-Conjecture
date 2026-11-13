"""Tests for research/t2-389a1/level3: the Fricke symmetry used to halve the level-3 sum, and (where a compiler
and AVX-512 are available) the compiled pipeline against bsd.kato_bound.riemann_coefficients at small p."""
import math
import os
import random
import shutil
import subprocess
import sys

import pytest

from bsd import kato_bound as kb
from bsd.pari_session import get_pari

pari = get_pari()
HERE = os.path.dirname(os.path.abspath(__file__))
LEVEL3 = os.path.join(os.path.dirname(HERE), "research", "t2-389a1", "level3")
sys.path.insert(0, LEVEL3)
C389 = [0, 1, 1, -2, 0]


def test_fricke_symmetry_values():
    """m'(b/M) = -w_N m'(b'/M) when N b b' = -1 mod M (run_level3.py docstring); -w_N = root number = +1."""
    E = pari.ellinit(C389)
    sym = kb.PlusSymbol(E)
    eps = int(pari.ellrootno(E))
    rng = random.Random(1)
    for M in (125, 343, 1000, 1331, 12345, 101 ** 2, 16231 ** 2):
        for _ in range(40):
            b = rng.randrange(1, M)
            while math.gcd(b, M) != 1:
                b = rng.randrange(1, M)
            bp = (-pow(389 * b, -1, M)) % M
            assert sym(b, M) == eps * sym(bp, M)


@pytest.mark.parametrize("p", [5, 7])
def test_reflection_of_level3_sums(p):
    """A(c - k) = A(k) for the level-3 sums, c = -log_gamma<389> mod p^2 (as computed by run_level3.reflection)."""
    import run_level3
    E = pari.ellinit(C389)
    sym = kb.PlusSymbol(E)
    q, p2 = p ** 3, p ** 2
    V = sym.values_mod(q)
    teich = [pow(t, p2, q) for t in range(1, p)]
    A, g = [], 1
    for k in range(p2):
        A.append(sum(V[w * g % q] for w in teich))
        g = g * (1 + p) % q
    c, kstar, eps = run_level3.reflection(p, E)
    assert (2 * kstar - c) % p == 0
    assert all(A[(c - k) % p2] == eps * A[k] for k in range(p2))


def _have_avx512():
    try:
        return "avx512f" in open("/proc/cpuinfo").read()
    except OSError:
        return False


@pytest.mark.skipif(shutil.which("gcc") is None or not _have_avx512(), reason="needs gcc and AVX-512")
def test_compiled_pipeline_small_p(tmp_path):
    import run_level3
    binary, table = str(tmp_path / "level3"), str(tmp_path / "phi389.bin")
    run_level3.TABLE, run_level3.BINARY = table, binary
    run_level3.table()
    subprocess.run(["gcc", "-O2", "-march=native", "-o", binary, os.path.join(LEVEL3, "level3.c")], check=True)
    E = pari.ellinit(C389)
    sym = kb.PlusSymbol(E)
    for p in (7, 11):
        ap = int(pari.ellap(E, p))
        c, kstar, eps = run_level3.reflection(p, E)
        rows = run_level3.parse(run_level3.run_c(p, ap, c, kstar, (p + 1) // 2))
        coeffs, _, covered = run_level3.coefficients(p, ap, rows, eps)
        assert len(covered) == p and all(r["hecke"] and r["hecke_r"] for r in rows)
        assert coeffs == kb.riemann_coefficients(sym, p, ap, 3, 4, P=4)["coeffs"]


@pytest.mark.skipif(shutil.which("gcc") is None or not _have_avx512(), reason="needs gcc and AVX-512")
def test_compiled_level2_small_p(tmp_path):
    """research/t2-389a1/beyond/level2.c: theta_1 = kato_bound's GP loop, and the ordinary level-2 Riemann sum
    = kato_bound.riemann_coefficients(n = 2), both mod p^3; every row passes its Hecke check."""
    import run_level3
    sys.path.insert(0, os.path.join(os.path.dirname(LEVEL3), "beyond"))
    import run_beyond
    table, binary = str(tmp_path / "phi389.bin"), str(tmp_path / "level2")
    run_level3.TABLE = run_beyond.TABLE = table
    run_beyond.BINARY = binary
    run_level3.table()
    subprocess.run(["gcc", "-O2", "-march=native", "-o", binary,
                    os.path.join(os.path.dirname(LEVEL3), "beyond", "level2.c")], check=True)
    E = pari.ellinit(C389)
    sym = kb.PlusSymbol(E)
    for p in (101, 107):
        d = run_beyond.level2(p, E)
        assert d["hecke_fail"] == 0 and d["fricke_fail"] == 0
        assert [x % p ** 3 for x in kb.theta1_coefficients(sym, p, 4)] == d["C"]
    for p in (7, 11):
        d = run_beyond.level2(p, E)
        c, _ = run_beyond.ordinary_coeffs(d)
        assert c == kb.riemann_coefficients(sym, p, int(pari.ellap(E, p)), 2, 4, P=3)["coeffs"]
    # the lean kernel (default build) against run_blocks (-DREFERENCE_KERNEL build): whole output lines
    ref = str(tmp_path / "level2_ref")
    subprocess.run(["gcc", "-O2", "-march=native", "-DREFERENCE_KERNEL", "-o", ref,
                    os.path.join(os.path.dirname(LEVEL3), "beyond", "level2.c")], check=True)
    for p in (101, 1009):
        a, b = run_beyond.level2(p, E), run_beyond.level2(p, E, binary=ref)
        assert (a["kernel"], b["kernel"]) == ("lean", "run_blocks")
        assert all(a[k] == b[k] for k in a if k not in ("seconds", "kernel"))
    # and against both run_blocks and the scalar msym on random b/q, up to the lean kernel's bound q < 2^40
    for q in (16231 ** 2, 1048573 ** 2):
        r = subprocess.run([binary, table, "vtest", str(q), "200000", "7"], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
