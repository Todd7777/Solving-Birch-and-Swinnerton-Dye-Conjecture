#!/usr/bin/env python3
"""Sha(389a1)[p] = 0 for primes beyond Stein-Wuthrich's bound 48859, prime by prime, with a compiled level-2 sum.

SW Thm 12.3 (cited) covers p = 2 and the good ordinary p < 48859 except 16231; ../sha_zero_report.md covers the rest
below 48859. From 48859 on (48859 is prime) every prime is certified here, by the same theorem chain as
bsd/sha_zero.py, with the level-2 data computed by level2.c (validated below against bsd/kato_bound.py):
  good ordinary p: c_2 = [T^2] of the level-2 Riemann sum (correct mod p) a unit gives ord_T L_p <= 2 (Kato, SW
    Thm 7.3/6.1: Sha[p^oo] finite); Lemma A gives rhobar_p surjective; Lemma B gives b_p = 0 with ord_p L_p^* = 0, no
    regulator needed (anomalous p, a_p = 1, included). So Sha(389a1)[p] = 0, and the main conjecture holds at p.
  supersingular p (a_p = 0): [T^2]L^- = -[T^2]theta_1 mod p a unit gives ord_T L_p <= 2 (kato_bound.SUPERSINGULAR_PROOF);
    then PARI's D_p-valued regulator, the saturation certificate and Lemma A, as in sha_zero.supersingular_bound.
A prime where c_2 = 0 mod p would need level 3 (see ../level3/); the certified range stops below the first such prime.

Usage: python3 research/t2-389a1/beyond/run_beyond.py build
       python3 research/t2-389a1/beyond/run_beyond.py validate
       python3 research/t2-389a1/beyond/run_beyond.py validate-lean   (lean kernel = run_blocks, whole output lines)
       python3 research/t2-389a1/beyond/run_beyond.py run --pmax 120000 [--shard i --nshards k]
       python3 research/t2-389a1/beyond/run_beyond.py aggregate --pmax 120000
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
T2 = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(T2))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(T2, "level3"))

from bsd import kato_bound as kb  # noqa: E402
from bsd import sha_zero as sz  # noqa: E402
from bsd.pari_session import get_pari  # noqa: E402
import run_level3  # noqa: E402

AINVS, GENS, PMIN = [0, 1, 1, -2, 0], [[-1, 1], [0, 0]], 48859
BINARY = os.path.join(HERE, "level2")
TABLE = run_level3.TABLE


def curve():
    pari = get_pari(stack=2 ** 30)
    # the D_p-valued regulator at supersingular p > 75000 needs more than 1 GB of PARI stack: allow it to grow
    pari.default("parisizemax", 12 * 10 ** 9)
    E = pari.ellinit(AINVS)
    return pari, E


def build():
    if not os.path.exists(TABLE):
        run_level3.table()
    subprocess.run(["gcc", "-O3", "-march=native", "-Wall", "-o", BINARY, os.path.join(HERE, "level2.c")], check=True)
    print("built", BINARY)


def level2(p, E, binary=None):
    pari = get_pari()
    ap = int(pari.ellap(E, p))
    c, _, eps = run_level3.reflection(p, E)
    r = subprocess.run([binary or BINARY, TABLE, str(p), str(ap), str(c % p), str(eps)], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"level2 p={p}: {r.stderr}")
    f = r.stdout.split()
    n = [int(x) for x in f[:19]]
    return dict(p=n[0], a_p=n[1], c=n[2], eps=n[3], S1=n[4], m0=n[5], rows=n[6], hecke_fail=n[7], fricke_fail=n[8],
                C=n[9:14], W=n[14:19], evals=int(f[19]), seconds=float(f[20]),
                kernel=f[21] if len(f) > 21 else "run_blocks")      # binaries before 2026-11-13 printed no tag


def ordinary_coeffs(d, P=3):
    p, pP = d["p"], d["p"] ** P
    ainv = pow(kb.unit_root(d["a_p"], p, P), -1, pP)
    z2, z3 = pow(ainv, 2, pP), pow(ainv, 3, pP)
    return [(z2 * d["C"][j] - z3 * d["S1"] * d["W"][j]) % pP for j in range(5)], kb.unit_root(d["a_p"], p, P)


def certify_prime(p, E):
    pari = get_pari()
    t0 = time.time()
    d = level2(p, E)
    out = dict(p=p, a_p=d["a_p"],
               level2={k: d[k] for k in ("rows", "hecke_fail", "fricke_fail", "evals", "seconds", "kernel")})
    checks_ok = d["hecke_fail"] == 0 and d["fricke_fail"] == 0
    if d["a_p"] % p == 0:                                     # supersingular, a_p = 0
        out["kind"] = "supersingular"
        th = d["C"]
        Lm = [(-x) % p for x in th]
        out["Lminus_mod_p"] = Lm
        out["theta1_c0_ok"] = th[0] % p ** 3 == (-(p - 1) * d["m0"]) % p ** 3
        out["T2_certified"] = bool(checks_ok and out["theta1_c0_ok"] and Lm[1] == 0 and Lm[2] != 0)
        if out["T2_certified"]:
            out["surjectivity"] = sz.rhobar_surjective_certificate(E, p)
            out["saturation"] = sz.saturation_certificate(E, GENS, p)
            out["bound"] = sz.supersingular_bound(E, p, Lm[2], 2, GENS)
            out["ok"] = bool(out["surjectivity"]["ok"] and out["saturation"]["ok"] and out["bound"]["sha_p_trivial"]
                             and out["bound"].get("consistent", True))
        else:
            out["ok"] = False
    else:
        out["kind"] = "ordinary"
        c, alpha = ordinary_coeffs(d)
        p3 = p ** 3
        out["coeffs_mod_p"] = [x % p for x in c[:5]]
        out["c0_exact"] = c[0] == ((1 - pow(alpha, -1, p3)) ** 2 * d["m0"]) % p3
        out["T2_certified"] = bool(checks_ok and out["c0_exact"] and c[1] % p == 0 and c[2] % p != 0)
        if out["T2_certified"]:
            out["surjectivity"] = sz.rhobar_surjective_certificate(E, p)
            out["bound"] = sz.ordinary_bound(E, p, 0, 2)          # Lemma B: no regulator, no basis needed
            out["ok"] = bool(out["surjectivity"]["ok"] and out["bound"]["sha_p_trivial"])
        else:
            out["ok"] = False
            out["reason"] = "c_2 = 0 mod p: needs level 3" if checks_ok and c[2] % p == 0 else "check failed"
    out["seconds"] = round(time.time() - t0, 1)
    return out


def validate():
    pari, E = curve()
    sym = kb.PlusSymbol(E)
    ok, rep = True, {}
    # (a) theta_1 against the GP loop of kato_bound (exact integers, compared mod p^3)
    rows = []
    for p in (101, 107, 211, 599, 1049, 2957):
        d = level2(p, E)
        ref = kb.theta1_coefficients(sym, p, 4)
        agree = [x % p ** 3 for x in ref] == d["C"] and d["hecke_fail"] == 0 and d["fricke_fail"] == 0
        rows.append(dict(p=p, agree=agree))
        ok &= agree
        print(f"(a) theta_1 at p = {p}: compiled = kato_bound: {agree}", flush=True)
    rep["theta1"] = rows
    # (b) ordinary level-2 Riemann sum against kato_bound.riemann_coefficients (n = 2), exactly mod p^3
    rows = []
    for p in (5, 7, 11, 13, 101, 211, 307, 1009):
        ap = int(pari.ellap(E, p))
        d = level2(p, E)
        c, _ = ordinary_coeffs(d)
        ref = kb.riemann_coefficients(sym, p, ap, 2, 4, P=3)["coeffs"]
        agree = c == ref and d["hecke_fail"] == 0 and d["fricke_fail"] == 0
        rows.append(dict(p=p, agree=agree))
        ok &= agree
        print(f"(b) level-2 sum at p = {p}: compiled {c} = Python {ref}: {agree}", flush=True)
    rep["riemann_level2"] = rows
    # (c) against the stored GP results: c_2 mod p (T2 run) and [T^2]L^- mod p (supersingular runs)
    stored = {}
    for fn in glob.glob(os.path.join(T2, "data", "t2_*shard*.jsonl")) + glob.glob(os.path.join(T2, "data", "t2_supersingular*.jsonl")):
        for line in open(fn):
            r = json.loads(line)
            stored[r["p"]] = r
    rows = []
    for p in (4931, 4993, 4999, 6661, 15101):
        d = level2(p, E)
        if d["a_p"] == 0:
            got, want = (-d["C"][2]) % p, stored[p]["Lminus_mod_p"][2]
        else:
            got, want = ordinary_coeffs(d)[0][2] % p, stored[p]["c2_mod_p"]
        rows.append(dict(p=p, got=got, stored=want, agree=got == want))
        ok &= got == want
        print(f"(c) p = {p}: compiled {got}, stored {want}: {got == want}", flush=True)
    rep["stored"] = rows
    rep["ok"] = bool(ok)
    with open(os.path.join(T2, "data", "beyond48859_validation.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    print("validation", "PASSED" if ok else "FAILED")
    return ok


def validate_lean(primes=(101, 1009, 4999, 15101, 20011, 30011, 48859, 58631, 75431, 78787, 99991, 100003)):
    """The lean kernel (default build) against run_blocks (-DREFERENCE_KERNEL build of the same source): the whole
    output line, every C_j and W_j mod p^3 included, must agree; then vtest at q up to 2^40."""
    import tempfile
    from concurrent.futures import ThreadPoolExecutor
    pari, E = curve()
    ref = os.path.join(tempfile.mkdtemp(), "level2_ref")
    subprocess.run(["gcc", "-O3", "-march=native", "-DREFERENCE_KERNEL", "-o", ref, os.path.join(HERE, "level2.c")],
                   check=True)
    keys = ("p", "a_p", "c", "eps", "S1", "m0", "rows", "hecke_fail", "fricke_fail", "C", "W", "evals")

    def both(p):
        a, b = level2(p, E), level2(p, E, binary=ref)
        return dict(p=p, a_p=a["a_p"], identical=all(a[k] == b[k] for k in keys) and a["kernel"] == "lean"
                    and b["kernel"] == "run_blocks", seconds_lean=a["seconds"], seconds_run_blocks=b["seconds"])

    with ThreadPoolExecutor(max_workers=4) as ex:
        rows = list(ex.map(both, sorted(primes, reverse=True)))
    rows.sort(key=lambda r: r["p"])
    for r in rows:
        print(f"p = {r['p']}: identical output {r['identical']} (lean {r['seconds_lean']} s, "
              f"run_blocks {r['seconds_run_blocks']} s)", flush=True)
    vt = []
    for q in (10 ** 6 + 3, 16231 ** 2, 524287 ** 2, 1048573 ** 2, 999999999989):
        r = subprocess.run([BINARY, TABLE, "vtest", str(q), "4000000", "20261113"], capture_output=True, text=True)
        vt.append(dict(q=q, ok=r.returncode == 0, line=r.stdout.strip()))
        print(r.stdout.strip(), flush=True)
    ok = all(r["identical"] for r in rows) and all(v["ok"] for v in vt)
    with open(os.path.join(T2, "data", "beyond_lean_validation.json"), "w") as fh:
        json.dump(dict(primes=rows, vtest=vt, ok=ok), fh, indent=1)
    print("lean validation", "PASSED" if ok else "FAILED")
    return ok


def shard_file(s):
    return os.path.join(T2, "data", f"beyond48859_shard{s}.jsonl")


def run(pmax, shard, nshards):
    pari, E = curve()
    primes = [int(q) for q in pari.primes([PMIN, pmax - 1])][shard::nshards]
    fn = shard_file(shard)
    done = {json.loads(l)["p"] for l in open(fn)} if os.path.exists(fn) else set()
    for p in primes:
        if p in done:
            continue
        rec = certify_prime(p, E)
        with open(fn, "a") as fh:
            fh.write(json.dumps(rec, default=str) + "\n")
        print(time.strftime("%H:%M:%S"), p, rec["kind"], rec["ok"], rec["seconds"], flush=True)


def aggregate(pmax):
    pari, E = curve()
    recs = {}
    for fn in glob.glob(os.path.join(T2, "data", "beyond48859_shard*.jsonl")):
        for line in open(fn):
            r = json.loads(line)
            recs[r["p"]] = r
    primes = [int(q) for q in pari.primes([PMIN, pmax - 1])]
    first_gap = next((p for p in primes if not recs.get(p, {}).get("ok")), None)
    reach = first_gap if first_gap is not None else pmax
    done = [p for p in primes if recs.get(p, {}).get("ok")]
    level3_needed = sorted(p for p, r in recs.items() if r.get("reason", "").startswith("c_2 = 0"))
    failures = sorted(p for p, r in recs.items() if not r.get("ok") and p not in level3_needed)
    ss = sorted(p for p, r in recs.items() if r.get("kind") == "supersingular")
    out = dict(pmin=PMIN, pmax=pmax, primes_in_range=len(primes), certified=len(done), every_prime_below=reach,
               first_gap=first_gap, level3_needed=level3_needed, other_failures=failures, supersingular=ss,
               anomalous=sorted(p for p, r in recs.items() if r.get("a_p") == 1),
               b_p_values=sorted({r.get("bound", {}).get("b_p") for r in recs.values() if r.get("ok")}),
               hecke_row_checks=sum(r["level2"]["rows"] for r in recs.values()),
               symbol_evaluations=sum(r["level2"]["evals"] for r in recs.values()),
               cpu_seconds=round(sum(r["seconds"] for r in recs.values()), 1))
    with open(os.path.join(T2, "data", "beyond48859_summary.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "validate", "validate-lean", "run", "aggregate"])
    ap.add_argument("--pmax", type=int, default=120000)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    a = ap.parse_args()
    if a.cmd == "build":
        build()
    elif a.cmd == "validate":
        sys.exit(0 if validate() else 1)
    elif a.cmd == "validate-lean":
        sys.exit(0 if validate_lean() else 1)
    elif a.cmd == "run":
        run(a.pmax, a.shard, a.nshards)
    else:
        aggregate(a.pmax)


if __name__ == "__main__":
    main()
