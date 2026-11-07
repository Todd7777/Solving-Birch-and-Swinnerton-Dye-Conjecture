#!/usr/bin/env python3
"""Compact the (large) sweep JSONL files into one gzipped CSV row per curve, suitable for the repository.

Columns: label, conductor, analytic_rank, rigorous_ord_lower, rigorous_ord_exact, rank_lower, rank_upper,
sha2_mod4_rank, torsion, tamagawa, omega, regulator, leading_term, L1_over_omega_plus_exact, sha_an, sha_an_int,
sha_is_square, cremona_match, p_part_proven_primes, p_part_not_covered, sage_prove_bsd_unproven, sage_sat_index, errors
"""
from __future__ import annotations

import csv
import glob
import gzip
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from bsd.report import COMPACT_FIELDS as COLS, compact_row  # noqa: E402


def row(r: dict) -> list:
    d = r if "rigorous_ord_lower" in r else compact_row(r)   # accept full or already-compact records
    return [d.get(c) for c in COLS]


def main(paths):
    for path in paths:
        out = path[:-6] + ".csv.gz"
        n = 0
        with open(path) as f, gzip.open(out, "wt", newline="") as g:
            w = csv.writer(g)
            w.writerow(COLS)
            recs = []
            for line in f:
                line = line.strip()
                if line:
                    recs.append(json.loads(line))
            recs.sort(key=lambda r: (r.get("invariants", {}).get("conductor", 0), r.get("label", "")))
            for r in recs:
                w.writerow(row(r))
                n += 1
        print(f"{out}: {n} rows, {os.path.getsize(out)} bytes")


if __name__ == "__main__":
    main(sys.argv[1:] or sorted(glob.glob(os.path.join(ROOT, "results", "sweep_*.jsonl"))))
