#!/usr/bin/env python3
"""Level-3 Riemann sum of L_p(389a1) in compiled code (level3.c), for the exceptional prime p = 16231.

At p = 16231 the T^2 coefficient of L_p is divisible by p (ord_p Reg_gamma = 1, see p16231.py), so Sha(389a1)[p]
can only be decided from [T^2]L_p modulo p^2, which the level-3 sum gives (kato_bound.coefficient_precision:
e_2 = 3 - 1 = 2). If ord_p [T^2]L_p = 1, then ord_T L_p = 2, Sha[p^oo] is finite (Kato, SW Thm 6.1), and Lemma B
of bsd/sha_zero.py gives b_p = 1, hence Sha(389a1)[p] = 0 (Cassels); with PARI's regulator (valuation 3,
generators saturated) b_p = 0 and the main conjecture holds at p.

Fricke symmetry (Lemma; used to halve the work). Let f be the newform of level N, w_N its Atkin-Lehner eigenvalue,
m' its plus modular symbol, and gcd(M, N) = gcd(b, M) = 1 with N b b' = -1 (mod M). Then m'(b/M) = -w_N m'(b'/M).
Proof: W = (0 -1; N 0) satisfies int_{Wx}^{Wy} f(u) du = int_x^y (f|W)(z) dz = w_N int_x^y f(z) dz, and W commutes with
z -> -conj(z) on cusps, so the plus symbol phi satisfies phi(W D) = w_N phi(D). Put
gamma = ((1 + N b b')/M, b; N b', M), an element of Gamma_0(N) of determinant 1 with gamma(0) = b/M and
gamma(-M/(N b')) = oo. Then -m'(b/M) = phi({b/M -> oo}) = phi(gamma {0 -> -M/(N b')}) = phi(W {oo -> b'/M})
= w_N m'(b'/M). For 389a1, -w_N is the root number +1. With b = omega(i) gamma^k one has
b' = omega(-N^-1 i^-1) gamma^(c - k), c = -log_gamma<N>, so A(c - k) = -w_N A(k), and likewise S2(c - k) = -w_N S2(k).
Hence only the classes k0 = k* + t mod p (t = 0..(p-1)/2, k* = c/2 mod p) are computed; each also yields its mirror
class c - k0 mod p. The reflection is checked end to end at small p (validate (c)) and directly at p = 16231 on
sample classes (validate (d)).

Steps (each also a subcommand):
  table      Manin-symbol tables from PARI (phi389.bin, regenerated, not tracked)
  build      gcc -O3 -march=native level3.c
  validate   (a) scalar continued-fraction evaluation = PARI mseval at random b/q, q = 16231^3 included;
             (b) AVX-512 kernel = scalar evaluation at random b/16231^3;
             (c) the whole level-3 pipeline = kato_bound.riemann_coefficients(n = 3) exactly mod p^4 at small p
  run        shards over k0 in [0, p), resumable: data/level3_p16231_shard<i>.txt, one line per k0
  aggregate  completeness, Hecke checks, T^j coefficients mod p^4 -> data/level3_p16231.json
  certify    if ord_p [T^2]L_p = 1: the Sha certificate at p (bsd/sha_zero.py) -> data/sha_zero_p16231.json

Usage: python3 research/t2-389a1/level3/run_level3.py {table,build,validate,run,aggregate} [--p P] [--nshards K]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import random
import struct
import subprocess
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
T2 = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(T2))
sys.path.insert(0, ROOT)

from bsd import kato_bound as kb  # noqa: E402
from bsd.pari_session import get_pari  # noqa: E402

AINVS = [0, 1, 1, -2, 0]
TABLE = os.path.join(HERE, "phi389.bin")
BINARY = os.path.join(HERE, "level3")


def curve():
    pari = get_pari(stack=2 ** 30)
    E = pari.ellinit(AINVS)
    return pari, E, kb.PlusSymbol(E)


def table():
    """PHI[s][c][d] = m'-value of the path g.oo -> g.0, g in SL_2(Z) with bottom row (c, (-1)^s d) mod N."""
    pari, E, sym = curve()
    N = sym.N

    def value(c, d):
        a, b, c1, d1 = kb.lift_to_sl2z(c, d, N)
        w = pari.mseval(sym.M, sym.xp, pari.matrix(2, 2, [a, b, c1, d1]))
        v = Fraction(int(pari.numerator(w)), int(pari.denominator(w))) / sym.g
        if v.denominator != 1:
            raise ArithmeticError("non-integral Manin symbol value")
        return int(v)

    vals = [0] * (2 * N * N)
    for c in range(N):
        for d in range(N):
            if c or d:
                vals[c * N + d] = value(c, d)
                vals[N * N + c * N + d] = value(c, (-d) % N)
    with open(TABLE, "wb") as fh:
        fh.write(struct.pack(f"<{len(vals)}i", *vals))
    print(f"table: N = {N}, values in [{min(vals)}, {max(vals)}]")


def build():
    subprocess.run(["gcc", "-O3", "-march=native", "-Wall", "-o", BINARY, os.path.join(HERE, "level3.c")], check=True)
    print("built", BINARY)


def run_c(*args, inp=None):
    r = subprocess.run([BINARY, TABLE, *map(str, args)], input=inp, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"level3 {args}: exit {r.returncode}: {r.stderr}")
    return r.stdout


def reflection(p, E):
    """c = -log_gamma<N> mod p^2 (gamma = 1 + p), k* = c/2 mod p, and eps = -w_N (the root number)."""
    pari = get_pari()
    x = pari(f"389 + O({p}^6)")
    lg = pari.log(x / pari.teichmuller(x)) / pari.log(pari(f"1 + {p} + O({p}^6)"))
    l = int(pari.lift(pari.Mod(pari.truncate(lg), p * p)))
    c = (-l) % (p * p)
    return c, (c % p) * pow(2, -1, p) % p, int(pari.ellrootno(E))


def coefficients(p, ap, lines, eps, P=4):
    """T^j coefficients (j <= 4) of the level-3 Riemann sum mod p^P from the per-k0 lines (each line gives class
    k0 and, if different, its mirror class k0r; covered classes are returned too)."""
    pP = p ** P
    ainv = pow(kb.unit_root(ap, p, P), -1, pP)
    z3, z4 = pow(ainv, 3, pP), pow(ainv, 4, pP)
    c = [0] * 5
    covered = set()
    for f in lines:
        parts = [(f["k0"], f["S2"], f["C"], f["W"])]
        if f["k0r"] != f["k0"]:
            parts.append((f["k0r"], f["S2r"], [eps * x for x in f["R"]], f["Wr"]))
        for k0, S2, C, W in parts:
            if k0 in covered:
                continue
            covered.add(k0)
            for j in range(5):
                c[j] = (c[j] + z3 * C[j] - z4 * S2 * W[j]) % pP
    return c, kb.unit_root(ap, p, P), covered


def parse(out):
    rows = []
    for line in out.splitlines():
        f = line.split()
        if len(f) != 27:
            continue
        n = [int(x) for x in f[:26]]
        rows.append(dict(k0=n[0], S2=n[1], hecke=n[2], C=n[3:8], W=n[8:13], k0r=n[13], S2r=n[14], hecke_r=n[15],
                         R=n[16:21], Wr=n[21:26], seconds=float(f[26])))
    return rows


def validate(p_big=16231):
    pari, E, sym = curve()
    rng = random.Random(20261111)
    report = {}
    # (a) scalar evaluation vs PARI
    qs = [p_big ** 3, p_big ** 2, p_big] + [p ** 3 for p in (5, 7, 11, 101, 307)] + [rng.randrange(2, 2 ** 50) for _ in range(10)]
    pairs = []
    for q in qs:
        for _ in range(300):
            b = rng.randrange(1, q)
            while math.gcd(b, q) != 1:
                b = rng.randrange(1, q)
            pairs.append((b, q))
    out = run_c("test", inp="\n".join(f"{b} {q}" for b, q in pairs)).split()
    bad = sum(int(v) != sym(b, q) for (b, q), v in zip(pairs, out))
    report["scalar_vs_pari"] = dict(pairs=len(pairs), mismatches=bad)
    print(f"(a) scalar continued fractions vs PARI mseval: {len(pairs)} values, {bad} mismatches", flush=True)
    # (b) vector kernel vs scalar
    out = run_c("vtest", p_big ** 3, 4_000_000, 20261111)
    report["vector_vs_scalar"] = out.strip()
    print("(b)", out.strip(), flush=True)
    # (c) whole pipeline, half the classes plus mirrors, vs kato_bound.riemann_coefficients at level 3
    small = []
    for p in (5, 7, 11, 13, 101, 211, 307):
        ap = int(pari.ellap(E, p))
        c0, kstar, eps = reflection(p, E)
        rows = parse(run_c(p, ap, c0, kstar, (p + 1) // 2))
        c, _, covered = coefficients(p, ap, rows, eps)
        ref = kb.riemann_coefficients(sym, p, ap, 3, 4, P=4)["coeffs"]
        hk = all(r["hecke"] and r["hecke_r"] for r in rows)
        ok = c == ref and hk and len(covered) == p
        small.append(dict(p=p, a_p=ap, c=c0, agree=c == ref, hecke_all=hk, classes=len(covered)))
        print(f"(c) p = {p}: compiled (half + mirrors) {c} vs Python {ref}: {'agree' if ok else 'DISAGREE'}", flush=True)
    report["pipeline_vs_python"] = small
    # (d) at p_big: a few classes computed directly equal the mirrors of their partners
    ap = int(pari.ellap(E, p_big))
    c0, kstar, eps = reflection(p_big, E)
    direct = []
    for t in (1, 2):
        rep_row = parse(run_c(p_big, ap, c0, (kstar + t) % p_big, 1))[0]
        mir_row = parse(run_c(p_big, ap, c0, rep_row["k0r"], 1))[0]
        same = (mir_row["C"] == [eps * x % p_big ** 4 for x in rep_row["R"]] and mir_row["W"] == rep_row["Wr"]
                and mir_row["S2"] == rep_row["S2r"] and mir_row["hecke"] and rep_row["hecke"])
        direct.append(dict(k0=rep_row["k0"], mirror=rep_row["k0r"], agree=same,
                           seconds=round(rep_row["seconds"] + mir_row["seconds"], 1)))
        print(f"(d) p = {p_big}: class {rep_row['k0']} mirrored to {rep_row['k0r']} = direct computation: {same}", flush=True)
    report["mirror_direct_checks"] = direct
    ok = (bad == 0 and "differing: 0" in report["vector_vs_scalar"] and all(s["agree"] and s["hecke_all"] for s in small)
          and all(d["agree"] for d in direct))
    report["ok"] = ok
    with open(os.path.join(T2, "data", "level3_validation.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    print("validation", "PASSED" if ok else "FAILED")
    return ok


def shard_file(p, s):
    return os.path.join(T2, "data", f"level3_p{p}_shard{s}.txt")


def run(p, nshards):
    pari, E, sym = curve()
    ap = int(pari.ellap(E, p))
    c0, kstar, eps = reflection(p, E)
    n = (p + 1) // 2                                    # representatives k* + t, t < n
    procs = []
    for s in range(nshards):
        lo, hi = n * s // nshards, n * (s + 1) // nshards
        fn = shard_file(p, s)
        done = set()
        if os.path.exists(fn):
            done = {r["k0"] for r in parse(open(fn).read())}
        t = lo
        while t < hi and (kstar + t) % p in done:
            t += 1
        if t >= hi:
            continue
        fh = open(fn, "a")
        procs.append((subprocess.Popen([BINARY, TABLE, str(p), str(ap), str(c0), str((kstar + t) % p), str(hi - t)],
                                       stdout=fh), fh))
        print(f"shard {s}: t in [{t}, {hi}) (k0 = k* + t mod p, k* = {kstar}, c = {c0})", flush=True)
    for pr, fh in procs:
        pr.wait()
        fh.close()


def aggregate(p, nshards):
    pari, E, sym = curve()
    ap = int(pari.ellap(E, p))
    rows = {}
    import glob
    for fn in sorted(glob.glob(os.path.join(T2, "data", f"level3_p{p}_shard*.txt"))):   # extra shards allowed
        for r in parse(open(fn).read()):
            rows[r["k0"]] = r
    c0, kstar, eps = reflection(p, E)
    reps = [(kstar + t) % p for t in range((p + 1) // 2)]
    missing = [k for k in reps if k not in rows]
    hecke_fail = sorted(k for k, r in rows.items() if not (r["hecke"] and r["hecke_r"]))
    out = dict(p=p, a_p=ap, reflection_c=c0, k_star=kstar, root_number=eps, representatives=len(reps),
               representatives_done=len(reps) - len(missing), hecke_failures=hecke_fail[:20],
               cpu_seconds=round(sum(r["seconds"] for r in rows.values()), 1))
    if not missing:
        c, alpha, covered = coefficients(p, ap, [rows[k] for k in reps], eps)
        out["classes_covered"] = len(covered)
        p2 = p * p
        out["coeffs_mod_p4"] = c
        out["coeffs_mod_p2"] = [x % p2 for x in c]
        out["c0_is_exact_constant"] = c[0] == ((1 - pow(alpha, -1, p ** 4)) ** 2 * sym(0, 1)) % p ** 4
        out["c1_zero_mod_p2"] = c[1] % p2 == 0
        out["c2_mod_p2"] = c[2] % p2
        out["vp_c2_at_least"] = kb.vp(c[2] % p2, p) if c[2] % p2 else 2
        out["c2_over_p_mod_p"] = (c[2] % p2) // p if c[2] % p == 0 else None
        out["certified_ord_p_T2_is_1"] = c[2] % p == 0 and (c[2] % p2) != 0
    with open(os.path.join(T2, "data", f"level3_p{p}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
    return out


def certify(p):
    """Sha(389a1)[p] = 0 from the level-3 result, by the chain of bsd/sha_zero.py (Lemmas A and B)."""
    from bsd import sha_zero as sz
    pari, E, sym = curve()
    lv = json.load(open(os.path.join(T2, "data", f"level3_p{p}.json")))
    p2 = p * p
    c = lv.get("coeffs_mod_p4")
    out = dict(p=p, level3=dict((k, lv.get(k)) for k in ("c2_mod_p2", "c0_is_exact_constant", "c1_zero_mod_p2",
                                                            "classes_covered", "hecke_failures")))
    checks = dict(complete=lv.get("classes_covered") == p and not lv.get("hecke_failures"),
                  c0_exact=bool(lv.get("c0_is_exact_constant")), c1_zero_mod_p2=bool(lv.get("c1_zero_mod_p2")),
                  c2_zero_mod_p=c is not None and c[2] % p == 0)
    vL = None if c is None or c[2] % p2 == 0 else kb.vp(c[2] % p2, p)
    out["checks"] = checks
    out["vp_T2_coefficient"] = vL if vL is not None else ">= 2"
    if not all(checks.values()) or vL != 1:
        out.update(ok=False, reason="ord_p [T^2]L_p = 1 not established")
    else:
        gens = [[-1, 1], [0, 0]]
        out["ord_T_Lp"] = 2                       # [T^2]L_p != 0 (valuation 1); >= 2 by Kato and rank 2
        out["surjectivity"] = sz.rhobar_surjective_certificate(E, p)
        out["saturation"] = sz.saturation_certificate(E, gens, p)
        out["bound"] = sz.ordinary_bound(E, p, vL, 2, gens)
        b = out["bound"]
        out["ok"] = bool(out["surjectivity"]["ok"] and b["sha_p_trivial"] and b.get("consistent", True))
        out["main_conjecture"] = bool(out["ok"] and out["saturation"]["ok"] and b.get("b_p_from_pari_regulator") == 0)
    with open(os.path.join(T2, "data", f"sha_zero_p{p}.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["table", "build", "validate", "run", "aggregate", "certify"])
    ap.add_argument("--p", type=int, default=16231)
    ap.add_argument("--nshards", type=int, default=4)
    a = ap.parse_args()
    if a.cmd == "table":
        table()
    elif a.cmd == "build":
        build()
    elif a.cmd == "validate":
        sys.exit(0 if validate(a.p) else 1)
    elif a.cmd == "run":
        run(a.p, a.nshards)
    elif a.cmd == "aggregate":
        aggregate(a.p, a.nshards)
    else:
        certify(a.p)


if __name__ == "__main__":
    main()
