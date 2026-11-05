"""Command-line interface.

  python -m bsd verify 11a1 37a1 389a1            # verify curves by Cremona label
  python -m bsd verify --ainvs 0,0,1,-7,6          # or by a-invariants
  python -m bsd curated --out results/curated      # the curated list (docs/05 is built from this)
  python -m bsd sweep --min 1 --max 1000 --workers 4 --out results/sweep_1_1000
  python -m bsd heuristic 11a1 37a1 389a1 5077a1 --X 200000
  python -m bsd padic 389a1 --primes 3,5,7
  python -m bsd report results/sweep_1_1000.jsonl --md results/sweep_1_1000.md
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import time


def _parse_ainvs(s: str) -> list[int]:
    return [int(t) for t in s.strip("[] ").split(",")]


def _curve_from_arg(arg: str):
    from . import ecdata
    if arg.startswith("[") or "," in arg:
        return arg, _parse_ainvs(arg), None
    c = ecdata.by_label(arg)
    return c.label, c.ainvs, c


def cmd_verify(a):
    from .verify import verify_curve
    from .report import curve_section, records_table, write_jsonl
    recs = []
    for arg in a.curves:
        label, ainvs, ref = _curve_from_arg(arg)
        gens = ref.gens if (ref is not None and a.use_known_gens) else None
        rec = verify_curve(ainvs, label=label, known_gens=gens, run_descent=not a.no_descent, effort=a.effort,
                           saturation_bound=a.sat_bound, descent_timeout=a.timeout, cremona=ref,
                           padic_bound=a.padic, padic_prec=a.padic_prec).to_dict()
        recs.append(rec)
        print(curve_section(rec))
    print(records_table(recs))
    if a.out:
        write_jsonl(a.out + ".jsonl", recs)
        with open(a.out + ".md", "w") as f:
            f.write(records_table(recs) + "\n\n" + "\n".join(curve_section(r) for r in recs))


def rec_wants_padic(c) -> bool:
    """p-adic L-functions need modular symbols of level N: only feasible for moderate conductors."""
    return (c.claimed_conductor or 0) <= 300_000


def cmd_curated(a):
    from .curated import SMALL, LARGE, resolve
    from .verify import verify_curve
    from .report import curve_section, records_table, write_jsonl, summary
    items = resolve(SMALL + (LARGE if not a.small_only else []))
    recs = []
    for c, ainvs, ref in items:
        t0 = time.time()
        gens = ref.gens if (ref is not None and a.use_known_gens) else None
        rec = verify_curve(ainvs, label=c.label, known_gens=gens, effort=c.descent_effort, saturation_bound=a.sat_bound,
                           descent_timeout=c.descent_timeout, cremona=ref,
                           padic_bound=(a.padic if rec_wants_padic(c) else 0), padic_prec=a.padic_prec).to_dict()
        rec["note"] = c.note
        if c.claimed_conductor is not None and rec["invariants"]["conductor"] != c.claimed_conductor:
            rec["errors"].append(f"claimed conductor {c.claimed_conductor} ≠ computed {rec['invariants']['conductor']}")
        if c.claimed_rank is not None and rec["analytic"]["analytic_rank"] != c.claimed_rank:
            rec["errors"].append(f"claimed rank {c.claimed_rank} ≠ analytic rank {rec['analytic']['analytic_rank']}")
        recs.append(rec)
        print(f"[{time.time()-t0:7.1f}s] {c.label:>18}  r_an={rec['analytic']['analytic_rank']}  Ш_an={rec['sha_an']}  {rec['errors'] or ''}", flush=True)
        if a.out:
            write_jsonl(a.out + ".jsonl", recs)
    if a.out:
        with open(a.out + ".md", "w") as f:
            f.write("## Summary\n\n" + summary(recs) + "\n\n## Table\n\n" + records_table(recs) + "\n\n## Curve by curve\n\n"
                    + "\n".join(curve_section(r) for r in recs))
    print(summary(recs))


# ---- sweep over Cremona's database, in parallel ----
_SWEEP_OPTS = {}


def _sweep_init(opts):
    global _SWEEP_OPTS
    _SWEEP_OPTS = opts
    from .pari_session import get_pari
    get_pari()


def _sweep_one(item):
    from .verify import verify_curve
    label, ainvs, ref = item
    o = _SWEEP_OPTS
    t0 = time.time()
    try:
        rec = verify_curve(ainvs, label=label, known_gens=(ref.gens if o["use_known_gens"] else None),
                           run_descent=not o["no_descent"], effort=o["effort"], saturation_bound=o["sat_bound"],
                           descent_timeout=o["timeout"], cremona=ref, ledger=o["ledger"]).to_dict()
    except BaseException as exc:  # noqa: BLE001 - keep the sweep alive
        rec = {"label": label, "errors": [f"fatal: {type(exc).__name__}: {str(exc)[:200]}"], "invariants": {"conductor": ref.conductor},
               "analytic": {"analytic_rank": -1}, "descent": {}, "sha_an": None, "sha_an_int": None, "sha_is_square": None,
               "sha_rel_err": None, "rank_status": "failed", "sha_status": "failed", "p_part": [], "cremona": None, "timings": {}}
    rec["wall"] = time.time() - t0
    return rec


def cmd_sweep(a):
    from . import ecdata
    curves = ecdata.load_range(a.min, a.max)
    if a.rank is not None:
        curves = [c for c in curves if c.rank == a.rank]
    if a.limit:
        curves = curves[: a.limit]
    done = set()
    outpath = a.out + ".jsonl"
    if a.resume and os.path.exists(outpath):
        with open(outpath) as f:
            for line in f:
                try:
                    done.add(json.loads(line)["label"])
                except Exception:
                    pass
    items = [(c.label, c.ainvs, c) for c in curves if c.label not in done]
    print(f"sweep: {len(curves)} curves in [{a.min},{a.max}], {len(items)} to do, {a.workers} workers", flush=True)
    opts = dict(use_known_gens=a.use_known_gens, no_descent=a.no_descent, effort=a.effort, sat_bound=a.sat_bound,
                timeout=a.timeout, ledger=not a.no_ledger)
    t0 = time.time()
    n = 0
    with open(outpath, "a" if a.resume else "w") as f:
        ctx = mp.get_context("fork")
        with ctx.Pool(a.workers, initializer=_sweep_init, initargs=(opts,)) as pool:
            for rec in pool.imap_unordered(_sweep_one, items, chunksize=4):
                f.write(json.dumps(rec) + "\n")
                f.flush()
                n += 1
                if n % 200 == 0 or rec.get("errors"):
                    msg = f"[{time.time()-t0:8.0f}s] {n}/{len(items)} {rec['label']}"
                    if rec.get("errors"):
                        msg += "  " + " | ".join(rec["errors"])[:160]
                    print(msg, flush=True)
    print(f"done in {time.time()-t0:.0f}s", flush=True)


def cmd_report(a):
    from .report import read_jsonl, records_table, summary, curve_section
    recs = read_jsonl(a.jsonl)
    recs.sort(key=lambda r: (r["invariants"].get("conductor", 0), r["label"]))
    md = "## Summary\n\n" + summary(recs) + "\n"
    if not a.summary_only:
        md += "\n## Table\n\n" + records_table(recs) + "\n"
    if a.sections:
        md += "\n## Curve by curve\n\n" + "\n".join(curve_section(r) for r in recs)
    if a.md:
        with open(a.md, "w") as f:
            f.write(md)
    print(summary(recs))


def cmd_heuristic(a):
    from .heuristic import partial_products, fitted_exponent
    rows = []
    for arg in a.curves:
        label, ainvs, ref = _curve_from_arg(arg)
        pts = partial_products(ainvs, a.X, n_samples=a.samples)
        slope = fitted_exponent(pts)
        rows.append((label, ref.rank if ref else None, slope, pts))
        print(f"{label}: rank {ref.rank if ref else '?'}; fitted exponent of ∏_(p≤X) N_p/p against log X: {slope:.3f}")
    if a.out:
        with open(a.out, "w") as f:
            json.dump([{"label": l, "rank": r, "fitted_exponent": s, "points": p} for l, r, s, p in rows], f)
    if a.plot:
        import math
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.5, 5))
        for label, r, slope, pts in rows:
            ax.plot([math.log(x) for x, _ in pts], [math.exp(v) for _, v in pts], label=f"{label} (rank {r}, fit {slope:.2f})")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("log X")
        ax.set_ylabel("∏_{p ≤ X, p ∤ N} N_p / p")
        ax.set_title("The Birch–Swinnerton-Dyer heuristic: slope ≈ rank on log–log axes")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=130)
        print("plot written to", a.plot)


def cmd_padic(a):
    from .padic import padic_bsd_dict
    label, ainvs, ref = _curve_from_arg(a.curve)
    gens = ref.gens if ref is not None else None
    out = []
    for p in [int(t) for t in a.primes.split(",")]:
        d = padic_bsd_dict(ainvs, p, n=a.prec, gens=gens)
        out.append(d)
        print(f"{label} p={p} ({d['reduction']}): r_p = {d['r_p']}, L_p = {d['L_p']}, R_p = {d['R_p']}, L_p/R_p = {d['sha_p']} {d['note']}")
    if a.out:
        with open(a.out, "w") as f:
            json.dump({"label": label, "results": out}, f, indent=1)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="bsd", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--no-descent", action="store_true", help="do not run the 2-descent (use known generators only)")
        p.add_argument("--use-known-gens", action="store_true", default=True)
        p.add_argument("--fresh", dest="use_known_gens", action="store_false", help="ignore Cremona's generators")
        p.add_argument("--effort", type=int, default=0)
        p.add_argument("--sat-bound", type=int, default=1000)
        p.add_argument("--timeout", type=float, default=120.0)
        p.add_argument("--out", default=None)
        p.add_argument("--padic", type=int, default=0, help="also run p-adic BSD for good primes 3 ≤ p ≤ this bound")
        p.add_argument("--padic-prec", type=int, default=5)

    p = sub.add_parser("verify"); p.add_argument("curves", nargs="+"); common(p); p.set_defaults(fn=cmd_verify)
    p = sub.add_parser("curated"); common(p); p.add_argument("--small-only", action="store_true"); p.set_defaults(fn=cmd_curated)
    p = sub.add_parser("sweep"); common(p)
    p.add_argument("--min", type=int, default=1); p.add_argument("--max", type=int, default=1000)
    p.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 1)); p.add_argument("--rank", type=int, default=None)
    p.add_argument("--limit", type=int, default=0); p.add_argument("--resume", action="store_true"); p.add_argument("--no-ledger", action="store_true")
    p.set_defaults(fn=cmd_sweep)
    p = sub.add_parser("report"); p.add_argument("jsonl"); p.add_argument("--md", default=None)
    p.add_argument("--summary-only", action="store_true"); p.add_argument("--sections", action="store_true"); p.set_defaults(fn=cmd_report)
    p = sub.add_parser("heuristic"); p.add_argument("curves", nargs="+"); p.add_argument("--X", type=int, default=100000)
    p.add_argument("--samples", type=int, default=40); p.add_argument("--out", default=None); p.add_argument("--plot", default=None); p.set_defaults(fn=cmd_heuristic)
    p = sub.add_parser("padic"); p.add_argument("curve"); p.add_argument("--primes", default="3,5,7"); p.add_argument("--prec", type=int, default=6)
    p.add_argument("--out", default=None); p.set_defaults(fn=cmd_padic)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
