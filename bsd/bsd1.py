"""Certified BSD (I), curve by curve: rank E(Q) = ord_{s=1} L(E,s).

For one curve per isogeny class (both sides are isogeny invariants) this module proves
rank E(Q) = ord_{s=1} L(E,s) = r, where r is the rank recorded in Cremona's tables, whenever
the following chain closes. It proves nothing about curves it is not run on, and nothing
about formula (II) or the finiteness of Sha beyond what is stated.

Proof chain
-----------
Write w for the root number (PARI ellrootno, exact) and c_k = L^(k)(E,1)/k!.

* Parity. Lambda(1+t) = w Lambda(1-t), so ord_{s=1} L(E,s) = (1-w)/2 mod 2. [functional equation;
  modularity: Wiles, Taylor-Wiles, Breuil-Conrad-Diamond-Taylor]
* Upper bound for ord. A ball (bsd/rigorous_l.py, Arb interval arithmetic with proven tail
  bounds) for c_r that excludes 0 gives ord <= r.
* Lower bound for ord.
    r = 0: nothing to prove.
    r = 1: w = -1 gives L(E,1) = 0.
    r = 2: w = +1, and a certified point of infinite order gives L(E,1) = 0, because L(E,1) != 0
           implies E(Q) finite [Kolyvagin 1988 with Gross-Zagier 1986 and Bump-Friedberg-
           Hoffstein / Murty-Murty; independently Kato 2004]. Parity gives L'(E,1) = 0.
    r = 3: w = -1, and two certified independent points give ord != 1, because ord = 1 implies
           rank E(Q) = 1 [Gross-Zagier 1986, Kolyvagin 1988/1990]. Parity rules out ord = 0, 2.
           (This is the logic of Buhler-Gross-Zagier 1985 for 5077a1, with Kolyvagin's theorem in
           place of the Heegner-point height bound; see bsd/gz_certificate.py for the latter.)
* Rank. r <= 1: ord = r implies rank E(Q) = r [Gross-Zagier-Kolyvagin; Kato for r = 0].
  r >= 2: rank >= r from r points certified independent modulo torsion (exact reduction
  argument, bsd/kato_bound.independence_certificate), and rank <= r from
    (a) 2-descent: PARI ellrank on ellrankinit(E). Its 2-Selmer computation depends on GRH only
        through the class group and units of the cubic 2-division field (PARI 2.17 ellrank.c:
        makevbnf calls Buchall, bnfselmer uses its S-units); curves with a rational 2-torsion
        point use Hilbert symbols and quadratic forms only. Every such cubic bnf is certified
        here with bnfcertify, so the bound is unconditional. The bound r2 = C - T - s subtracts
        the 2-rank s of Sha[2]/2Sha[4] from Cassels-Tate pairing computations; the plain
        Selmer bound C - T = r2 + s is recorded as well;
    (b) optionally, for small conductors, Kato + Mazur control + exact p-adic L-function
        (bsd/kato_bound.py) at a good ordinary p >= 5 for E without CM:
        corank Sel_{p^oo}(E/Q) <= ord_T L_p <= r. Independent of (a).

Trusted software: PARI/GP 2.17 (ellrootno, ellan, ellrank, bnfcertify, msfromell/mseval,
reduction of points mod q) and FLINT/Arb (ball arithmetic). Trusted theorems: as cited above.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from fractions import Fraction
from typing import Sequence

from .pari_session import get_pari

ASSUMPTIONS = [
    "Modularity of E/Q (Wiles; Taylor-Wiles; Breuil-Conrad-Diamond-Taylor 2001): L(E,s) is entire with "
    "Lambda(s) = w Lambda(2-s), w = root number.",
    "Gross-Zagier 1986 + Kolyvagin 1988/1990 (+ Bump-Friedberg-Hoffstein 1990 / Murty-Murty 1991 for the "
    "auxiliary imaginary quadratic field): ord_{s=1} L(E,s) <= 1 implies rank E(Q) = ord. Kato 2004 "
    "independently for ord = 0.",
    "Kato 2004 Thm 17.4 and Mazur's control theorem, only for the optional p-adic bound (b).",
    "Correctness of PARI/GP 2.17 and FLINT/Arb (python-flint) as used.",
]

STACK = 2 ** 31


# ------------------------------------------------------------------------------------------
# pieces
# ------------------------------------------------------------------------------------------
def on_curve(ainvs: Sequence[int], P: Sequence[Fraction]) -> bool:
    a1, a2, a3, a4, a6 = ainvs
    x, y = Fraction(P[0]), Fraction(P[1])
    return y * y + a1 * x * y + a3 * y == x ** 3 + a2 * x * x + a4 * x + a6


def descent_bound(E) -> dict:
    """Unconditional 2-descent upper bound for the rank (see the module docstring, (a))."""
    pari = get_pari()
    t0 = time.time()
    R = pari.ellrankinit(E)
    vbnf = R[2]
    cubic = []
    for K in vbnf:
        # degree <= 2 factors are stored as nf (no class group used); degree 3 as bnf
        pol = pari("(K)->K.pol")(K)
        if int(pari.poldegree(pol)) == 3:
            cubic.append(dict(pol=str(pol), disc=int(pari("(K)->K.disc")(K)),
                              class_number=int(pari("(K)->K.no")(K)), certified=int(pari.bnfcertify(K)) == 1))
    res = pari.ellrank(R)
    r1, r2, s = int(res[0]), int(res[1]), int(res[2])
    ok = all(c["certified"] for c in cubic)
    return dict(rank_lower=r1, rank_upper=r2, sha2_rank_s=s, selmer_bound_C_minus_T=r2 + s,
                cubic_fields=cubic, unconditional=ok, seconds=round(time.time() - t0, 3),
                method="PARI ellrank(ellrankinit(E)) with bnfcertify on the cubic 2-division field")


def kato_upper_bound(E, N: int, target: int, primes: Sequence[int] = (5, 7, 11, 13, 17, 19, 23, 29)) -> dict:
    """Optional independent bound (b): first good ordinary p >= 5 giving ord_T L_p <= target.

    p = 3 is left out because Stein-Wuthrich's Thm 7.3/7.4 cover every case only for p >= 5. CM curves
    are skipped (outside their running hypothesis)."""
    from . import kato_bound as kb
    pari = get_pari()
    t0 = time.time()
    if int(pari.elliscm(E)) != 0:
        return dict(ok=False, error="CM curve: skipped", seconds=0.0)
    try:
        sym = kb.PlusSymbol(E, N)
    except Exception as exc:  # noqa: BLE001  (PARI stack overflow at large N)
        return dict(ok=False, error=f"{type(exc).__name__}: {str(exc)[:120]}", seconds=round(time.time() - t0, 2))
    tried = []
    for p in primes:
        if N % p == 0 or int(pari.ellap(E, p)) % p == 0:
            continue
        res = kb.ordinary_bound(E, p, target=target, sym=sym, max_evaluations=10 ** 6)
        tried.append(dict(p=p, ord_upper_bound=res["ord_upper_bound"], level=res.get("level")))
        r = res["ord_upper_bound"]
        if r is not None and r <= target:
            return dict(ok=True, p=p, ord_upper_bound=r, level=res["level"], tried=tried,
                        seconds=round(time.time() - t0, 2))
    return dict(ok=False, tried=tried, seconds=round(time.time() - t0, 2))


def analytic_order(ainvs, r: int, lower_reason: str | None, prec: int) -> dict:
    from .rigorous_l import certify_order_rigorous
    ext = (r, lower_reason) if lower_reason else None
    for pr in (prec, 2 * prec, 4 * prec):
        res = certify_order_rigorous(ainvs, kmax=r, prec=pr, use_modular_symbols=False,
                                     external_lower_bound=ext, digits=20)
        if res["order_certified"]:
            break
    out = {k: res[k] for k in ("order", "order_certified", "order_candidates", "order_lower", "order_upper",
                               "conclusion", "root_number", "M_terms", "prec_bits")}
    out["leading_ball"] = res["c_conditioned"][r]["ball"] if res["order_certified"] else None
    out["c_r_ball"] = res["c"][r]["ball"]
    return out


# ------------------------------------------------------------------------------------------
# one class
# ------------------------------------------------------------------------------------------
def certify_class(label: str, ainvs: Sequence[int], rank: int, gens: Sequence[Sequence[Fraction]],
                  torsion_order: int, *, prec: int = 64, kato_max_conductor: int = 0) -> dict:
    """Certify rank E(Q) = ord_{s=1} L(E,s) = rank (Cremona's value) for one curve."""
    pari = get_pari(stack=STACK)
    t0 = time.time()
    E = pari.ellinit([int(a) for a in ainvs])
    N = int(pari.ellglobalred(E)[0])
    w = int(pari.ellrootno(E))
    out = dict(label=label, conductor=N, ainvs=[int(a) for a in ainvs], cremona_rank=rank, root_number=w)
    problems = []
    if rank >= 4:
        problems.append("rank >= 4: no known method certifies ord >= 4")
    if (-1) ** rank != w:
        problems.append("root number contradicts the parity of Cremona's rank")
    pts = [tuple(P) for P in gens[:rank]]
    if rank >= 2:
        if len(pts) < rank:
            problems.append(f"only {len(pts)} generators available")
        if not all(on_curve(ainvs, P) for P in pts):
            problems.append("a generator is not on the curve")
    lower_reason = None
    indep = None
    if rank >= 2 and not problems:
        indep = pari_independence(E, pts, torsion_order)
        out["independence"] = indep
        if not indep["ok"]:
            problems.append("independence of the generators not certified")
        elif rank == 2:
            lower_reason = ("a certified point of infinite order: L(E,1) != 0 would force E(Q) finite "
                            "(Kolyvagin/Kato), so L(E,1) = 0; w = +1")
        elif rank == 3:
            lower_reason = ("two certified independent points: ord = 1 would force rank 1 "
                            "(Gross-Zagier-Kolyvagin), so ord >= 3 by parity (w = -1)")
    analytic = None
    if not problems:
        analytic = analytic_order(ainvs, rank, lower_reason, prec)
        out["analytic"] = analytic
        if not (analytic["order_certified"] and analytic["order"] == rank):
            problems.append("order of vanishing not certified")
    if rank >= 2 and not problems:
        db = descent_bound(E)
        out["descent"] = db
        kato = None
        if kato_max_conductor and N <= kato_max_conductor:
            kato = kato_upper_bound(E, N, rank)
            out["kato"] = kato
        ok_descent = db["unconditional"] and db["rank_upper"] <= rank
        ok_kato = bool(kato and kato.get("ok"))
        out["upper_bound_sources"] = ([f"2-descent (r2 = {db['rank_upper']})"] if ok_descent else []) + \
                                     ([f"Kato p = {kato['p']}"] if ok_kato else [])
        if not (ok_descent or ok_kato):
            problems.append("no unconditional rank upper bound <= rank")
        # Sha[p^oo] finiteness that comes for free
        sha = []
        if ok_descent and db["selmer_bound_C_minus_T"] == rank:
            sha.append(dict(p=2, statement="Sha[2] = 0, hence Sha[2^oo] = 0 (dim Sel_2 = rank + dim E(Q)[2])"))
        if ok_kato:
            sha.append(dict(p=kato["p"], statement="Sha[p^oo] finite (corank Sel_{p^oo} <= rank)"))
        out["sha_p_finite"] = sha
    out["bsd1_certified"] = not problems
    out["problems"] = problems
    if not problems:
        out["statement"] = f"rank E(Q) = ord_(s=1) L(E,s) = {rank}"
    out["seconds"] = round(time.time() - t0, 3)
    return out


def pari_independence(E, pts, torsion_order):
    from .kato_bound import independence_certificate
    return independence_certificate(E, pts, torsion_order)


# ------------------------------------------------------------------------------------------
# batch runner over Cremona's tables (resumable, parallel)
# ------------------------------------------------------------------------------------------
def load_classes(nmin: int, nmax: int, cache_dir: str | None = None) -> list[dict]:
    """Curve number 1 of every isogeny class with nmin <= N <= nmax, with Cremona's rank, torsion
    order and generators."""
    from . import ecdata
    kw = {} if cache_dir is None else {"cache_dir": cache_dir}
    rows: dict[str, dict] = {}
    for block in range((nmin // 10000) * 10000, nmax + 1, 10000):
        with open(ecdata.fetch("allcurves", block, **kw)) as f:
            for line in f:
                t = line.split()
                N, iso, num = int(t[0]), t[1], int(t[2])
                if num != 1 or not nmin <= N <= nmax:
                    continue
                rows[f"{N}{iso}1"] = dict(label=f"{N}{iso}1", N=N, ainvs=[int(a) for a in t[3].strip("[]").split(",")],
                                          rank=int(t[4]), torsion_order=int(t[5]), gens=[])
        with open(ecdata.fetch("allgens", block, **kw)) as f:
            for line in f:
                t = line.split()
                lab = f"{t[0]}{t[1]}{t[2]}"
                if lab in rows:
                    r = int(t[4])
                    gens = []
                    for tok in t[6:6 + r]:
                        X, Y, Z = (int(v) for v in tok.strip("[]").split(":"))
                        gens.append((Fraction(X, Z), Fraction(Y, Z)))
                    rows[lab]["gens"] = gens
    return sorted(rows.values(), key=lambda d: (d["N"], d["label"]))


def _worker(args):
    row, prec, kato_max = args
    try:
        res = certify_class(row["label"], row["ainvs"], row["rank"], row["gens"], row["torsion_order"],
                            prec=prec, kato_max_conductor=kato_max)
    except Exception as exc:  # noqa: BLE001
        res = dict(label=row["label"], error=f"{type(exc).__name__}: {str(exc)[:300]}", bsd1_certified=False)
    return json.dumps(res, default=str)


def run(nmin: int, nmax: int, out: str, workers: int = 3, prec: int = 64, kato_max_conductor: int = 3000,
        log=print) -> None:
    import multiprocessing as mp
    rows = load_classes(nmin, nmax)
    done = set()
    if os.path.exists(out):
        with open(out) as f:
            for line in f:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if "error" not in r:
                    done.add(r["label"])
    todo = [r for r in rows if r["label"] not in done]
    log(f"{len(rows)} classes with {nmin} <= N <= {nmax}; {len(done)} done; {len(todo)} to do")
    # big ranks first so that slow items do not straggle at the end
    todo.sort(key=lambda r: (-r["rank"], -r["N"]))
    t0 = time.time()
    with mp.get_context("fork").Pool(workers, maxtasksperchild=2000) as pool, open(out, "a") as fh:
        for i, line in enumerate(pool.imap_unordered(_worker, [(r, prec, kato_max_conductor) for r in todo],
                                                     chunksize=16), start=1):
            fh.write(line + "\n")
            if i % 2000 == 0 or i == len(todo):
                fh.flush()
                log(f"{i}/{len(todo)} classes, {time.time() - t0:.0f} s")


def summarize(path: str) -> dict:
    from collections import Counter, defaultdict
    recs = {}
    with open(path) as f:
        for line in f:
            r = json.loads(line)
            if r["label"] not in recs or "error" in recs[r["label"]]:
                recs[r["label"]] = r
    by_rank = defaultdict(Counter)
    failures = []
    for r in recs.values():
        k = r.get("cremona_rank", "?")
        by_rank[k]["classes"] += 1
        if r.get("bsd1_certified"):
            by_rank[k]["certified"] += 1
        else:
            failures.append(dict(label=r["label"], problems=r.get("problems"), error=r.get("error")))
        if r.get("kato", {}).get("ok"):
            by_rank[k]["kato_bound"] += 1
        d = r.get("descent")
        if d:
            by_rank[k]["descent_unconditional"] += bool(d["unconditional"] and d["rank_upper"] <= k)
            by_rank[k]["descent_needs_cassels"] += bool(d["selmer_bound_C_minus_T"] > k)
            by_rank[k]["cubic_fields_certified"] += len([c for c in d["cubic_fields"] if c["certified"]])
            by_rank[k]["sha2_zero"] += bool(d["selmer_bound_C_minus_T"] == k)
        by_rank[k]["seconds"] += r.get("seconds", 0)
        a = r.get("analytic")
        if a:
            by_rank[k][f"prec_{a['prec_bits']}"] += 1
    Ns = [r["conductor"] for r in recs.values() if "conductor" in r]
    return dict(classes=len(recs), conductor_range=[min(Ns), max(Ns)] if Ns else None,
                by_rank={str(k): dict(v) for k, v in sorted(by_rank.items(), key=lambda kv: str(kv[0]))},
                failures=failures)


def _cli(argv=None):  # pragma: no cover
    import argparse
    ap = argparse.ArgumentParser(description="Certify rank = analytic rank for every isogeny class in a range.")
    ap.add_argument("--min", type=int, default=1)
    ap.add_argument("--max", type=int, default=9999)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--prec", type=int, default=64)
    ap.add_argument("--kato-max-conductor", type=int, default=3000)
    ap.add_argument("--summary", action="store_true", help="only summarise an existing output file")
    a = ap.parse_args(argv)
    if not a.summary:
        run(a.min, a.max, a.out, a.workers, a.prec, a.kato_max_conductor,
            log=lambda *m: print(time.strftime("%H:%M:%S"), *m, flush=True))
    print(json.dumps(summarize(a.out), indent=1)[:20000])


if __name__ == "__main__":  # pragma: no cover
    _cli(sys.argv[1:])


# ------------------------------------------------------------------------------------------
# compact exports (the raw JSONL is large and is not tracked)
# ------------------------------------------------------------------------------------------
CSV_FIELDS = ["label", "conductor", "rank", "root_number", "bsd1_certified", "leading_coefficient_ball", "prec_bits",
              "M_terms", "independence_ell", "independence_primes", "descent_rank_upper", "selmer2_minus_T",
              "cubic_field_certified", "kato_p", "sha_p_finite"]


def export(raw: str, outdir: str) -> dict:
    """Write classes.csv.gz (one row per class), rank2plus.jsonl.gz (full records, rank >= 2) and
    summary.json from a raw bsd1 JSONL file."""
    import csv
    import gzip
    recs = {}
    with open(raw) as f:
        for line in f:
            r = json.loads(line)
            recs[r["label"]] = r
    rows = sorted(recs.values(), key=lambda r: (r.get("conductor", 0), r["label"]))
    os.makedirs(outdir, exist_ok=True)
    with gzip.open(os.path.join(outdir, "classes.csv.gz"), "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_FIELDS)
        for r in rows:
            a, ind, d, k = r.get("analytic") or {}, r.get("independence") or {}, r.get("descent") or {}, r.get("kato") or {}
            cub = d.get("cubic_fields")
            w.writerow([r["label"], r.get("conductor"), r.get("cremona_rank"), r.get("root_number"),
                        int(bool(r.get("bsd1_certified"))), a.get("leading_ball"), a.get("prec_bits"), a.get("M_terms"),
                        ind.get("ell", ""), " ".join(map(str, ind.get("primes", []))), d.get("rank_upper", ""),
                        d.get("selmer_bound_C_minus_T", ""),
                        "" if cub is None else ("none" if not cub else int(all(c["certified"] for c in cub))),
                        k.get("p", "") if k.get("ok") else "",
                        " ".join(str(s["p"]) for s in r.get("sha_p_finite", []))])
    with gzip.open(os.path.join(outdir, "rank2plus.jsonl.gz"), "wt") as fh:
        for r in rows:
            if (r.get("cremona_rank") or 0) >= 2:
                fh.write(json.dumps(r, default=str) + "\n")
    summ = summarize(raw)
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1)
    return summ
