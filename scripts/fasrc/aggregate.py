#!/usr/bin/env python3
"""Aggregate the compact sweep outputs in results/fasrc/ into a small summary that can be committed.

Writes results/fasrc/summary.json, results/fasrc/summary.md, results/fasrc/anomalies.jsonl (every curve
with an error, a mismatch with Cremona's tables, a non-square or non-integral Sha_an, or an analytic rank
outside the 2-descent bounds) and results/fasrc/rank_ge2.csv.gz (all curves of analytic rank >= 2).
The full per-curve files are large; keep them on the cluster (or archive them elsewhere), not in git.
"""
from __future__ import annotations

import csv
import glob
import gzip
import json
import os
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "results", "fasrc")


def anomalous(r: dict) -> list[str]:
    why = []
    if r.get("errors"):
        why.append("error")
    if r.get("cremona_match") is False:
        why.append("cremona-mismatch")
    if r.get("sha_is_square") is False:
        why.append("sha-not-square")
    if (r.get("sha_rel_err") or 0) > 1e-6:
        why.append("sha-not-integral")
    lo, hi, ra = r.get("rank_lower"), r.get("rank_upper"), r.get("analytic_rank")
    if None not in (lo, hi, ra) and not (lo <= ra <= hi):
        why.append("rank-outside-descent-bounds")
    return why


def main(kind: str = "sweep"):
    files = sorted(glob.glob(os.path.join(OUT, f"{kind}_*.jsonl")))
    if not files:
        sys.exit(f"no {kind}_*.jsonl files in {OUT}")
    blocks, total = [], Counter()
    ranks, shas, status, unproven = Counter(), Counter(), Counter(), Counter()
    seen = set()
    anomalies, rank2 = [], []
    maxrel, maxsha = 0.0, []
    for path in files:
        n = 0
        bc = Counter()
        for line in open(path):
            try:
                r = json.loads(line)
            except Exception:
                bc["unparseable"] += 1
                continue
            if r["label"] in seen:
                continue
            seen.add(r["label"])
            n += 1
            ranks[r.get("analytic_rank")] += 1
            if r.get("sha_an_int") is not None:
                shas[r["sha_an_int"]] += 1
                maxsha.append((r["sha_an_int"], r["label"]))
            status[r.get("rank_status", "")] += 1
            maxrel = max(maxrel, r.get("sha_rel_err") or 0.0)
            for q in (r.get("sage_prove_bsd_unproven") or "").split():
                if q != "none":
                    unproven[int(q)] += 1
            if r.get("sage_prove_bsd_unproven") == "none":
                bc["full_bsd_proven"] += 1
            why = anomalous(r)
            if why:
                bc["anomalies"] += 1
                anomalies.append({"why": why, **r})
            if (r.get("analytic_rank") or 0) >= 2:
                rank2.append(r)
        blocks.append({"file": os.path.basename(path), "curves": n, **bc})
        total["curves"] += n
        total.update({k: v for k, v in bc.items()})
    maxsha.sort(reverse=True)
    summary = {"kind": kind, "files": len(files), "curves": total["curves"], "anomalies": len(anomalies),
               "analytic_ranks": dict(sorted((str(k), v) for k, v in ranks.items())),
               "sha_an": {str(k): v for k, v in sorted(shas.items())}, "rank_status": dict(status),
               "max_rel_deviation_of_sha_an_from_integer": maxrel, "largest_sha_an": maxsha[:20],
               "full_bsd_proven": total.get("full_bsd_proven", 0), "prove_bsd_unproven_primes": dict(unproven),
               "blocks": blocks}
    with open(os.path.join(OUT, f"{kind}_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    with open(os.path.join(OUT, f"{kind}_anomalies.jsonl"), "w") as f:
        for a in anomalies:
            f.write(json.dumps(a) + "\n")
    if rank2:
        from bsd.report import COMPACT_FIELDS  # noqa: E402
        with gzip.open(os.path.join(OUT, f"{kind}_rank_ge2.csv.gz"), "wt", newline="") as g:
            w = csv.writer(g)
            w.writerow(COMPACT_FIELDS)
            for r in sorted(rank2, key=lambda r: (r["conductor"], r["label"])):
                w.writerow([r.get(c) for c in COMPACT_FIELDS])
    md = [f"# FASRC {kind} summary", "",
          f"- curves: {summary['curves']} in {len(files)} files",
          f"- analytic ranks: {summary['analytic_ranks']}",
          f"- Sha_an values: {summary['sha_an']}",
          f"- max relative deviation of Sha_an from an integer: {maxrel:.2e}",
          f"- anomalies (errors, mismatches with Cremona, non-square/non-integral Sha_an, rank outside descent bounds): {len(anomalies)}",
          f"- status counts: {dict(status)}"]
    if kind == "prove":
        md.append(f"- full BSD formula proven (Sage prove_BSD) for {summary['full_bsd_proven']} curves; unproven primes: {dict(unproven)}")
    md.append(f"- largest Sha_an: {maxsha[:10]}")
    with open(os.path.join(OUT, f"{kind}_summary.md"), "w") as f:
        f.write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    sys.path.insert(0, ROOT)
    main(sys.argv[1] if len(sys.argv) > 1 else "sweep")
