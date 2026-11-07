#!/usr/bin/env python3
"""Re-run sweep records whose exact L(E,1) status is missing (e.g. PARI ran out of stack at a highly
composite level in an older version of the sweep) using the current class-based method, and rewrite
the JSONL file atomically with the repaired records in place.

usage: python scripts/repair_records.py results/sweep_1001_10000.jsonl [--workers 2]
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from bsd import ecdata  # noqa: E402
from bsd import cli  # noqa: E402


def needs_repair(r: dict) -> bool:
    inv, an = r.get("invariants", {}), r.get("analytic", {})
    if r.get("errors"):
        return True
    if inv.get("root_number") != 1:
        return False
    return an.get("L1_over_omega_plus_exact") is None and not an.get("L1_status")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jsonl")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    recs = [json.loads(l) for l in open(a.jsonl) if l.strip()]
    bad = [r for r in recs if needs_repair(r)]
    print(f"{len(recs)} records, {len(bad)} need repair")
    if not bad:
        return
    by_label = {}
    for N in sorted({r["invariants"]["conductor"] for r in bad}):
        for c in ecdata.load_range(N, N):
            by_label[c.label] = c
    classes: dict = {}
    for r in bad:
        c = by_label[r["label"]]
        classes.setdefault((c.conductor, c.iso), []).append((c.label, c.ainvs, c))
    opts = dict(use_known_gens=True, no_descent=False, effort=0, sat_bound=1000, timeout=120.0, ledger=True,
                use_sage=False, sage_timeout=600.0, padic_bound=0, padic_prec=5)
    fixed = {}
    with mp.get_context("fork").Pool(a.workers, initializer=cli._sweep_init, initargs=(opts,)) as pool:
        for out in pool.imap_unordered(cli._sweep_class, list(classes.values()), chunksize=1):
            for rec in out:
                fixed[rec["label"]] = rec
    tmp = a.jsonl + ".repair_tmp"
    with open(tmp, "w") as f:
        for r in recs:
            f.write(json.dumps(fixed.get(r["label"], r)) + "\n")
    os.replace(tmp, a.jsonl)
    still = sum(1 for r in fixed.values() if needs_repair(r))
    print(f"repaired {len(fixed)} records; {still} still lack an exact L(E,1) status")


if __name__ == "__main__":
    main()
