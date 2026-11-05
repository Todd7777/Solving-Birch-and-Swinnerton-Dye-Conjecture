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
COLS = ["label", "conductor", "analytic_rank", "rigorous_ord_lower", "rigorous_ord_exact", "rank_lower", "rank_upper",
        "sha2_mod4_rank", "torsion", "tamagawa", "omega", "regulator", "leading_term", "L1_over_omega_plus_exact",
        "sha_an", "sha_an_int", "sha_is_square", "cremona_match", "p_part_proven_primes", "p_part_not_covered",
        "sage_prove_bsd_unproven", "sage_sat_index", "errors"]


def row(r: dict) -> list:
    inv, an, des = r.get("invariants", {}), r.get("analytic", {}), r.get("descent", {})
    pp = r.get("p_part") or []
    sg = r.get("sage") or {}
    pb = sg.get("prove_bsd") or {}
    st = sg.get("saturation") or {}
    crem = r.get("cremona")
    return [r.get("label"), inv.get("conductor"), an.get("analytic_rank"), an.get("rigorous_order_lower"),
            an.get("rigorous_order_exact"), des.get("rank_lower"), des.get("rank_upper"), des.get("sha2_mod_4_rank"),
            inv.get("torsion_order"), inv.get("tamagawa_product"), inv.get("omega"), des.get("regulator"),
            an.get("leading_term"), an.get("L1_over_omega_plus_exact"), r.get("sha_an"), r.get("sha_an_int"),
            r.get("sha_is_square"), (None if crem is None else not crem.get("mismatches")),
            " ".join(str(p["p"]) for p in pp if p["status"] == "proven"),
            " ".join(str(p["p"]) for p in pp if p["status"] == "not covered"),
            ("" if pb.get("primes_unproven") is None else " ".join(map(str, pb["primes_unproven"])) or "none"),
            st.get("index"), " | ".join(r.get("errors") or [])]


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
