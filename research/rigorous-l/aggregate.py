#!/usr/bin/env python3
"""Aggregate validation.jsonl, timing.jsonl and certify.jsonl into results.json (and print report tables).

Usage: python3 research/rigorous-l/aggregate.py
Restart-safe: it only reads the .jsonl files.  The last good record per key wins.
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from bsd.rigorous_l import ASSUMPTIONS  # noqa: E402


def load(name, key):
    path = os.path.join(HERE, name)
    out = {}
    errors = {}
    if not os.path.exists(path):
        return out, errors
    with open(path) as f:
        for line in f:
            try:
                r = json.loads(line)
            except Exception:
                continue
            if "error" in r:
                errors[r[key]] = r
            else:
                out[r[key]] = r
                errors.pop(r[key], None)
    return out, errors


def compact_cert(r):
    rr = (r.get("cremona") or {}).get("rank")
    lead = None
    if r.get("order_certified"):
        lead = r["c_conditioned"][r["order"]]["ball"]
    ms = r.get("modular_symbols") or None
    return {
        "label": r["label"], "conductor": r["conductor"], "root_number": r["root_number"], "sets": r["sets"],
        "cremona_rank": rr, "order": r["order"], "order_certified": r["order_certified"],
        "order_candidates": r["order_candidates"], "order_lower": r["order_lower"], "order_upper": r["order_upper"],
        "c_zero_exact": r["c_zero_exact"], "c_nonzero_certified": r["c_nonzero_certified"],
        "c": [x["ball"] for x in r["c"]],
        "leading_coefficient_ball": lead,
        "modular_symbols": ({"ratio_abs": ms.get("ratio_abs"), "zero": ms.get("zero"), "seconds": ms.get("seconds")}
                            if ms and ms.get("computed") else (ms if ms else None)),
        "undetermined": r.get("undetermined") or {},
        "L_prime_abs_upper": r.get("L_prime_abs_upper"),
        "L_prime_abs_upper_prec512": r.get("L_prime_abs_upper_prec512"),
        "rel_diff_leading_vs_cremona": (r.get("cremona") or {}).get("rel_diff_c_r_vs_cremona"),
        "M_terms": r["M_terms"], "seconds": r["seconds"],
        "conclusion": r["conclusion"],
    }


def main():
    val, val_err = load("validation.jsonl", "key")
    tim, _ = load("timing.jsonl", "key")
    cert, cert_err = load("certify.jsonl", "label")
    r3, _ = load("rank3.jsonl", "label")
    xc, _ = load("crosscheck.jsonl", "key")
    bey, _ = load("beyond.jsonl", "label")

    # ---------------- validation
    vrows = []
    n_checked = n_ok = 0
    for key, r in val.items():
        for c in r["coefficients"]:
            n_checked += 2
            n_ok += int(c["contains_pari"]) + int(c["lambda_contains_pari"])
        vrows.append({"key": key, "label": r["label"], "conductor": r["conductor"], "root_number": r["root_number"],
                      "cremona_rank": r["cremona_rank"], "prec_bits": r["prec_bits"], "M_terms": r["M_terms"],
                      "seconds": r["seconds"], "all_contained": r["all_contained"], "max_rad": r["max_rad"],
                      "coefficients": [{"k": c["k"], "ball": c["c"]["ball"], "pari": c["pari_lfun_over_kfact"][:72],
                                        "contains_pari": c["contains_pari"],
                                        "abs_mid_minus_pari": c["abs_mid_minus_pari"],
                                        "exact_zero_parity": c["c_exact_zero_parity"]} for c in r["coefficients"]]})
    validation = {
        "reference": "PARI/GP 2.17 lfun(E,1,k)/k! and lfunlambda(E,1,k)/(2 k!) at 450 bits (~135 digits); "
                     "floating point, not rigorous",
        "n_curves": len({r['label'] for r in val.values()}), "n_runs": len(val),
        "n_values_checked": n_checked, "n_values_inside_ball": n_ok,
        "all_contained": all(r["all_contained"] for r in val.values()) if val else None,
        "max_rad_prec128": max((r["max_rad"] for r in val.values() if r["prec_bits"] == 128), default=None),
        "max_rad_prec256": max((r["max_rad"] for r in val.values() if r["prec_bits"] == 256), default=None),
        "runs": sorted(vrows, key=lambda x: (x["prec_bits"], x["conductor"])),
        "errors": list(val_err),
    }

    # ---------------- timing
    timing = sorted(({"label": r["label"], "conductor": r["conductor"], "prec_bits": r["prec_bits"],
                      "M_terms": r["M_terms"], "seconds": r["seconds"], "max_rad_c": r["max_rad_c"]}
                     for r in tim.values()), key=lambda x: (x["prec_bits"], x["conductor"]))

    # ---------------- certification
    by_rank = defaultdict(Counter)
    by_set = defaultdict(Counter)
    rank3 = []
    not_cert = []
    disagreements = []
    modsym_times = []
    rel_diffs = []
    for lab, r in cert.items():
        rr = (r.get("cremona") or {}).get("rank")
        key = f"cremona_rank_{rr}" if rr is not None else "no_cremona_entry"
        status = ("certified" if r["order_certified"]
                  else ("candidates_1_3" if r["order_candidates"] == [1, 3] else "not_certified"))
        by_rank[key]["total"] += 1
        by_rank[key][status] += 1
        for s in r["sets"]:
            by_set[s]["total"] += 1
            by_set[s][status] += 1
        if r["order_certified"] and rr is not None and r["order"] != rr:
            disagreements.append(lab)
        if r["order_candidates"] == [1, 3]:
            rank3.append({"label": lab, "conductor": r["conductor"], "c3_ball": r["c"][3]["ball"],
                          "L_prime_abs_upper_prec128": r.get("L_prime_abs_upper"),
                          "L_prime_abs_upper_prec512": r.get("L_prime_abs_upper_prec512"),
                          "c3_ball_prec512": r.get("c3_prec512")})
        elif not r["order_certified"]:
            not_cert.append({"label": lab, "candidates": r["order_candidates"], "conclusion": r["conclusion"],
                             "modular_symbols": r.get("modular_symbols")})
        ms = r.get("modular_symbols")
        if ms and ms.get("computed"):
            modsym_times.append((ms["seconds"], r["conductor"], lab))
        d = (r.get("cremona") or {}).get("rel_diff_c_r_vs_cremona")
        if d is not None:
            rel_diffs.append(d)
    modsym_times.sort(reverse=True)
    enc_times = sorted(r["timings"]["total_s"] for r in cert.values())
    tot_times = sorted(r["timings"].get("certify_total_s", r["seconds"]) for r in cert.values())

    def _q(v, f):
        return v[min(len(v) - 1, int(f * len(v)))] if v else None
    certification = {
        "prec_bits": 128, "kmax": 4,
        "sets": {
            "curated": "every curve of results/curated.jsonl with conductor < 20000",
            "rank2_class_rep / rank3_class_rep": "curve number 1 of every isogeny class with Cremona rank 2 or 3 "
                                                 "and conductor <= 10000 (allcurves.00000-09999, column 5)",
            "random_rank01": "random sample (seed 20261107) of 200 curves of Cremona rank 0 or 1, conductor <= 10000",
        },
        "n_curves": len(cert), "n_errors": len(cert_err),
        "summary_by_cremona_rank": {k: dict(v) for k, v in sorted(by_rank.items())},
        "summary_by_set": {k: dict(v) for k, v in sorted(by_set.items())},
        "certified_order_disagrees_with_cremona_rank": disagreements,
        "rank3_candidates_1_3": rank3,
        "not_certified_other": not_cert,
        "modular_symbols": {"n_computed": len(modsym_times),
                            "total_seconds": round(sum(t for t, _, _ in modsym_times), 1),
                            "slowest": [{"seconds": t, "conductor": N, "label": l} for t, N, l in modsym_times[:5]]},
        "max_rel_diff_leading_vs_cremona_table": max(rel_diffs) if rel_diffs else None,
        "enclosure_seconds": {"median": _q(enc_times, 0.5), "p90": _q(enc_times, 0.9),
                              "max": enc_times[-1] if enc_times else None},
        "certificate_seconds_incl_modsym": {"median": _q(tot_times, 0.5), "p90": _q(tot_times, 0.9),
                                            "max": tot_times[-1] if tot_times else None},
        "errors": [{"label": l, "error": e.get("error")} for l, e in cert_err.items()],
        "curves": [compact_cert(r) for r in sorted(cert.values(), key=lambda r: (r["conductor"], r["label"]))],
    }
    out = {
        "meta": {
            "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
            "module": "bsd/rigorous_l.py",
            "what": "Rigorous Arb ball enclosures of c_k = L^(k)(E,1)/k!, k = 0..4, and certified orders of vanishing "
                    "of L(E,s) at s = 1.  Nothing here proves BSD: these are statements about the analytic side only.",
            "ball_format": "Arb decimal strings '[mid +/- rad]'. Each printed ball is guaranteed to contain the true "
                           "value (arb_get_str); '0' is an exact zero.",
            "assumptions": ASSUMPTIONS,
        },
        "validation": validation,
        "timing": timing,
        "certification": certification,
        "crosscheck_two_rigorous_methods": {
            "note": "lambda_k = Lambda^(k)(1)/k! by two independent rigorous algorithms: (a) the K_k series of "
                    "taylor_enclosures, (b) Arb's certified integrator acb.integral on int_1^U g_M(u)(log u)^k du "
                    "plus proven u- and n-tails (lambda_by_integration). Balls must overlap.",
            "n_curves": len(xc), "all_overlap": all(r["all_overlap"] for r in xc.values()) if xc else None,
            "runs": [{"label": r["label"], "conductor": r["conductor"], "root_number": r["root_number"],
                      "all_overlap": r["all_overlap"], "seconds_series": r["seconds_series"],
                      "seconds_integral": r["seconds_integral"],
                      "max_rad_series": max(c["rad_series"] for c in r["coefficients"]),
                      "max_rad_integral": max(c["rad_integral"] for c in r["coefficients"])}
                     for r in sorted(xc.values(), key=lambda r: r["conductor"])],
        },
        "beyond_N_20000_record_ranks": {
            "note": "curated record-rank curves with N >= 20000. Only UPPER bounds are certifiable here: c_r != 0 "
                    "gives ord <= r. Modular symbols were not run at these levels. Lower bounds ord >= 4 are not known "
                    "to be provable for any curve.",
            "curves": [{"label": r["label"], "conductor": r["conductor"], "root_number": r["root_number"],
                        "numerical_analytic_rank_pari": r.get("numerical_analytic_rank_pari"),
                        "kmax": r["kmax"], "M_terms": r["M_terms"], "seconds": r["seconds"],
                        "order_upper": r["order_upper"], "order_candidates": r["order_candidates"],
                        "c_nonzero_certified": r["c_nonzero_certified"],
                        "undetermined": {k: v["c_abs_upper"] for k, v in r["undetermined"].items()},
                        "c_top_ball": r["c"][r["kmax"]]["ball"]}
                       for r in sorted(bey.values(), key=lambda r: r["conductor"])],
        },
        "rank3_all_N_below_20000": {
            "note": "all Cremona rank-3 curves with N < 20000. Each: c_0 = 0 exactly (parity), c_3 != 0 certified, "
                    "c_1 = L'(E,1) only bounded (a ball cannot prove 0). The bound is what a Gross-Zagier "
                    "certificate must turn into L'(E,1) = 0; with it, ord = 3 follows.",
            "curves": [r3[k] for k in sorted(r3, key=lambda l: r3[l]["conductor"])],
        },
    }
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    # short printout for the report
    print(json.dumps({"validation": {k: v for k, v in validation.items() if k != "runs"},
                      "summary_by_cremona_rank": certification["summary_by_cremona_rank"],
                      "summary_by_set": certification["summary_by_set"],
                      "disagreements": disagreements, "n_errors": len(cert_err),
                      "rank3": rank3, "not_certified_other": not_cert[:10],
                      "rank3_all": {l: (v["by_prec"]["128"]["L_prime_abs_upper"], v["by_prec"]["512"]["L_prime_abs_upper"],
                                        v["by_prec"]["512"]["seconds"]) for l, v in r3.items()},
                      "modsym": certification["modular_symbols"],
                      "max_rel_diff_vs_cremona": certification["max_rel_diff_leading_vs_cremona_table"]},
                     indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
