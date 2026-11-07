#!/usr/bin/env python3
"""Resumable batch driver for research/rigorous-l (validation, timing, certification).

Usage:  python3 research/rigorous-l/run_batch.py [validation|timing|certify|rank3|crosscheck|beyond|all]

Every phase appends one JSON line per item to its own .jsonl file in this
directory and skips items already present.  A restart or kill loses at most
the item in progress.  One process; PARI stack 2 GiB (the cap for this project).
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from flint import arb, ctx  # noqa: E402

from bsd.pari_session import get_pari  # noqa: E402
from bsd.rigorous_l import (ball_to_json, certify_order_rigorous, taylor_enclosures,  # noqa: E402
                            lambda_by_integration, _dec_upper)

ECDIR = os.path.join(ROOT, "data", "ecdata")
SEED = 20261107
PARI_BITS = 450          # PARI reference values at 450 bits (~135 digits)
STACK = 2 ** 31          # 2 GiB PARI stack (project cap)


def log(msg: str) -> None:
    with open(os.path.join(HERE, "batch.log"), "a") as f:
        f.write(time.strftime("%H:%M:%S ") + msg + "\n")


def done_keys(path: str, key: str = "key", skip_errors: bool = False) -> set:
    out = set()
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                except Exception:
                    continue      # a partially written last line (killed mid-write) is ignored
                if skip_errors and "error" in rec:
                    continue      # failed items are retried on the next run
                out.add(rec[key])
    return out


def ensure_newline(path: str) -> None:
    """After a kill mid-write the file may end in a partial line; terminate it so new records stay parseable."""
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path, "rb") as f:
            f.seek(-1, os.SEEK_END)
            last = f.read(1)
        if last != b"\n":
            with open(path, "a") as f:
                f.write("\n")


def append(path: str, rec: dict) -> None:
    with open(path, "a") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        f.flush()
        os.fsync(f.fileno())


def parse_allcurves(block: str):
    rows = []
    with open(os.path.join(ECDIR, f"allcurves.{block}")) as f:
        for line in f:
            p = line.split()
            if len(p) < 6:
                continue
            rows.append({"label": p[0] + p[1] + p[2], "N": int(p[0]), "iso": p[1], "num": int(p[2]),
                         "ainvs": [int(t) for t in p[3].strip("[]").split(",")], "rank": int(p[4]),
                         "torsion": int(p[5])})
    return rows


def cremona_bsd(labels: set) -> dict:
    """Cremona's allbsd data (rank, L^{(r)}(1)/r!) for the given labels."""
    out = {}
    for block in ("00000-09999", "10000-19999"):
        with open(os.path.join(ECDIR, f"allbsd.{block}")) as f:
            for line in f:
                p = line.split()
                if len(p) < 11:
                    continue
                lab = p[0] + p[1] + p[2]
                if lab in labels:
                    out[lab] = {"rank": int(p[4]), "L_value": p[8], "omega": p[7], "sha_an": p[10]}
    return out


def pari_reference(pari, E, k: int, lam: bool = False) -> tuple[arb, str]:
    v = pari.lfunlambda(E, 1, k, precision=PARI_BITS) if lam else pari.lfun(E, 1, k, precision=PARI_BITS)
    s = str(v).replace(" E", "e")
    with ctx.workprec(PARI_BITS):
        x = arb(s) / math.factorial(k)
        if lam:
            x = x / 2          # PARI's Lambda has gamma factor Gamma_R(s)Gamma_R(s+1) = 2 (2pi)^-s Gamma(s)
    return x, s


# ----------------------------------------------------------------------------
def phase_validation(pari):
    path = os.path.join(HERE, "validation.jsonl")
    ensure_newline(path)
    done = done_keys(path)
    rows = parse_allcurves("00000-09999") + parse_allcurves("10000-19999")
    bylab = {r["label"]: r for r in rows}
    named = ["11a1", "37a1", "389a1", "5077a1", "571a1"]
    rnd = random.Random(SEED)
    sample = [r["label"] for r in rnd.sample(rows, 20)]
    items = [(lab, 128) for lab in named] + [(lab, 256) for lab in named] + [(lab, 128) for lab in sample]
    for lab, prec in items:
        key = f"{lab}@{prec}"
        if key in done:
            continue
        r = bylab[lab]
        E = pari.ellinit(r["ainvs"])
        t = time.perf_counter()
        enc = taylor_enclosures(r["ainvs"], kmax=4, prec=prec, E=E)
        dt = time.perf_counter() - t
        ks = []
        for k in range(5):
            pv, ps = pari_reference(pari, E, k)
            lv, ls = pari_reference(pari, E, k, lam=True)
            ck, lk = enc.c[k], enc.lam[k]
            with ctx.workprec(PARI_BITS):
                if ck.is_zero():
                    c_ok = abs(float(pv.mid())) < 1e-60
                    c_diff = abs(pv.mid())
                else:
                    c_ok = bool(ck.contains(pv))
                    c_diff = abs(ck.mid() - pv.mid())
                if lk.is_zero():
                    l_ok = abs(float(lv.mid())) < 1e-60
                else:
                    l_ok = bool(lk.contains(lv))
            ks.append({"k": k, "c": ball_to_json(ck, 45), "pari_lfun_over_kfact": pv.str(70, radius=False),
                       "pari_string": ps[:80], "contains_pari": c_ok, "abs_mid_minus_pari": float(c_diff),
                       "c_exact_zero_parity": bool(ck.is_zero()),
                       "lambda": lk.str(45, radius=True),
                       "pari_lfunlambda_over_2kfact": lv.str(50, radius=False), "lambda_contains_pari": l_ok})
        rec = {"key": key, "label": lab, "ainvs": r["ainvs"], "conductor": enc.conductor,
               "root_number": enc.root_number, "cremona_rank": r["rank"], "prec_bits": prec, "M_terms": enc.M,
               "seconds": round(dt, 4), "coefficients": ks,
               "all_contained": all(x["contains_pari"] and x["lambda_contains_pari"] for x in ks),
               "max_rad": max(x["c"]["rad_float"] for x in ks)}
        append(path, rec)
        log(f"validation {key}: all_contained={rec['all_contained']} max_rad={rec['max_rad']:.2e} t={dt:.3f}s")


# ----------------------------------------------------------------------------
def phase_timing(pari):
    path = os.path.join(HERE, "timing.jsonl")
    ensure_newline(path)
    done = done_keys(path)
    rows = parse_allcurves("00000-09999") + parse_allcurves("10000-19999")
    targets = [11, 37, 100, 389, 1000, 2000, 5077, 7500, 10000, 12500, 15000, 17500, 19990]
    picks = []
    for T in targets:
        r = next((r for r in rows if r["N"] >= T and r["num"] == 1), None)
        if r is not None and r["label"] not in [p["label"] for p in picks]:
            picks.append(r)
    for r in picks:
        for prec in (128, 256, 512):
            key = f"{r['label']}@{prec}"
            if key in done:
                continue
            E = pari.ellinit(r["ainvs"])
            t = time.perf_counter()
            enc = taylor_enclosures(r["ainvs"], kmax=4, prec=prec, E=E)
            dt = time.perf_counter() - t
            rec = {"key": key, "label": r["label"], "conductor": enc.conductor, "root_number": enc.root_number,
                   "prec_bits": prec, "M_terms": enc.M, "seconds": round(dt, 4), "timings": enc.timings,
                   "max_rad_c": max(float(x.rad()) for x in enc.c)}
            append(path, rec)
            log(f"timing {key}: N={enc.conductor} M={enc.M} t={dt:.3f}s")


# ----------------------------------------------------------------------------
def certify_items():
    """(label, ainvs, sets) for the certification run, in processing order."""
    items: dict[str, dict] = {}
    order: list[str] = []

    def add(lab, ainvs, tag, N):
        if lab not in items:
            items[lab] = {"label": lab, "ainvs": ainvs, "sets": [], "N": N}
            order.append(lab)
        if tag not in items[lab]["sets"]:
            items[lab]["sets"].append(tag)

    with open(os.path.join(ROOT, "results", "curated.jsonl")) as f:
        for line in f:
            rec = json.loads(line)
            inv = rec["invariants"]
            if inv["conductor"] < 20000:
                add(rec["label"], inv["ainvs"], "curated", inv["conductor"])
    rows = parse_allcurves("00000-09999")
    pool = [r for r in rows if r["rank"] in (0, 1)]
    rnd = random.Random(SEED)
    for r in rnd.sample(pool, 200):
        add(r["label"], r["ainvs"], "random_rank01", r["N"])
    reps = sorted([r for r in rows if r["num"] == 1 and r["rank"] in (2, 3)], key=lambda r: (r["N"], r["iso"]))
    for r in reps:
        add(r["label"], r["ainvs"], f"rank{r['rank']}_class_rep", r["N"])
    return [items[l] for l in order]


def phase_certify(pari, prec: int = 128, modsym_timeout: float = 900.0):
    path = os.path.join(HERE, "certify.jsonl")
    ensure_newline(path)
    done = done_keys(path, "label", skip_errors=True)
    items = certify_items()
    crem = cremona_bsd({it["label"] for it in items})
    todo = [it for it in items if it["label"] not in done]
    log(f"certify: {len(items)} curves, {len(todo)} to do")
    t_start = time.time()
    for i, it in enumerate(todo):
        lab = it["label"]
        t = time.perf_counter()
        try:
            E = pari.ellinit(it["ainvs"])
            res = certify_order_rigorous(it["ainvs"], kmax=4, prec=prec, E=E, label=lab,
                                         modsym_timeout=modsym_timeout, digits=40)
            res.pop("assumptions", None)          # stated once in results.json
            res.pop("c_conditioned_note", None)
            cr = crem.get(lab)
            if cr is not None:
                res["cremona"] = dict(cr)
                rr = cr["rank"]
                if rr <= 4:
                    with ctx.workprec(prec + 64):
                        cmid = res["c_conditioned"][rr]["mid"]
                        cm = arb(cmid) if cmid != "0" else arb(0)
                        Lc = arb(cr["L_value"])
                        rel = float(abs(cm - Lc) / abs(Lc)) if not Lc.is_zero() else None
                    res["cremona"]["rel_diff_c_r_vs_cremona"] = rel
                res["cremona_rank_agrees"] = (res["order"] == rr) if res["order_certified"] else None
            if res["root_number"] == -1 and res.get("order_candidates") == [1, 3]:
                hi = certify_order_rigorous(it["ainvs"], kmax=4, prec=512, E=E, label=lab, digits=60)
                res["L_prime_abs_upper_prec512"] = hi.get("L_prime_abs_upper")
                res["c3_prec512"] = hi["c"][3]["ball"]
            res["sets"] = it["sets"]
            res["seconds"] = round(time.perf_counter() - t, 4)
            append(path, res)
        except Exception as exc:
            append(path, {"label": lab, "ainvs": it["ainvs"], "sets": it["sets"], "error": f"{type(exc).__name__}: {exc}",
                          "traceback": traceback.format_exc()[-1500:], "seconds": round(time.perf_counter() - t, 4)})
            log(f"certify {lab}: ERROR {type(exc).__name__}: {str(exc)[:200]}")
        if i % 25 == 0 or time.perf_counter() - t > 60:
            log(f"certify [{i + 1}/{len(todo)}] {lab} N={it['N']} t={time.perf_counter() - t:.2f}s "
                f"elapsed={time.time() - t_start:.0f}s")
    log("certify: phase complete")


# ----------------------------------------------------------------------------
def phase_rank3(pari):
    """All curves of Cremona rank 3 with N < 20000: rigorous bounds |L'(E,1)| <= eps at 128/256/512 bits."""
    path = os.path.join(HERE, "rank3.jsonl")
    ensure_newline(path)
    done = done_keys(path, "label", skip_errors=True)
    rows = parse_allcurves("00000-09999") + parse_allcurves("10000-19999")
    for r in [r for r in rows if r["rank"] == 3]:
        if r["label"] in done:
            continue
        E = pari.ellinit(r["ainvs"])
        res = certify_order_rigorous(r["ainvs"], kmax=4, prec=128, E=E, label=r["label"], digits=40)
        rec = {"label": r["label"], "ainvs": r["ainvs"], "conductor": res["conductor"],
               "root_number": res["root_number"], "order_candidates": res["order_candidates"],
               "order_certified": res["order_certified"], "conclusion": res["conclusion"], "by_prec": {}}
        for prec in (128, 256, 512):
            t = time.perf_counter()
            enc = taylor_enclosures(r["ainvs"], kmax=4, prec=prec, E=E)
            rec["by_prec"][str(prec)] = {
                "seconds": round(time.perf_counter() - t, 3), "M_terms": enc.M,
                "L_prime_abs_upper": _dec_upper(enc.c[1]),
                "c2_abs_upper": _dec_upper(enc.c[2]),
                "c3_ball": enc.c[3].str(prec // 4, radius=True),
                "c3_nonzero": not enc.c[3].contains(0),
                "lambda1_abs_upper": _dec_upper(enc.lam[1]),
                "lambda3_ball": enc.lam[3].str(prec // 4, radius=True)}
        append(path, rec)
        log(f"rank3 {r['label']}: candidates={rec['order_candidates']} "
            f"|L'(E,1)| <= {rec['by_prec']['512']['L_prime_abs_upper']} (512 bits)")


# ----------------------------------------------------------------------------
def phase_crosscheck(pari):
    """Two independent rigorous methods for lambda_k: K_k series (taylor_enclosures) vs acb.integral."""
    path = os.path.join(HERE, "crosscheck.jsonl")
    ensure_newline(path)
    done = done_keys(path)
    rows = parse_allcurves("00000-09999") + parse_allcurves("10000-19999")
    bylab = {r["label"]: r for r in rows}
    named = ["11a1", "37a1", "389a1", "5077a1", "571a1"]
    rnd = random.Random(SEED)
    sample = [r["label"] for r in rnd.sample(rows, 20)]          # the same 20 curves as phase_validation
    for lab in named + sample:
        key = f"{lab}@128"
        if key in done:
            continue
        r = bylab[lab]
        E = pari.ellinit(r["ainvs"])
        t = time.perf_counter()
        enc = taylor_enclosures(r["ainvs"], kmax=4, prec=128, E=E)
        t1 = time.perf_counter()
        li = lambda_by_integration(r["ainvs"], kmax=4, prec=128, E=E)
        t2 = time.perf_counter()
        ks = [{"k": k, "series": enc.lam[k].str(35, radius=True), "integral": li[k].str(35, radius=True),
               "overlap": bool(li[k].overlaps(enc.lam[k])), "rad_series": float(enc.lam[k].rad()),
               "rad_integral": float(li[k].rad())} for k in range(5)]
        rec = {"key": key, "label": lab, "conductor": enc.conductor, "root_number": enc.root_number,
               "seconds_series": round(t1 - t, 4), "seconds_integral": round(t2 - t1, 4), "coefficients": ks,
               "all_overlap": all(x["overlap"] for x in ks)}
        append(path, rec)
        log(f"crosscheck {key}: all_overlap={rec['all_overlap']} series {t1 - t:.2f}s integral {t2 - t1:.2f}s")


# ----------------------------------------------------------------------------
def phase_beyond(pari):
    """Curated record-rank curves with N >= 20000 (ranks 4, 5, 6): rigorous UPPER bounds ord <= r.

    Modular symbols are not attempted at these levels.  Rank 7 (N ~ 3.8e11, ~1.1e7 terms) is skipped:
    it is feasible but needs ~0.8 GB and ~40 min with the current non-streaming a_n handling.
    """
    path = os.path.join(HERE, "beyond.jsonl")
    ensure_newline(path)
    done = done_keys(path, "label", skip_errors=True)
    with open(os.path.join(ROOT, "results", "curated.jsonl")) as f:
        recs = [json.loads(line) for line in f]
    for rec in recs:
        inv = rec["invariants"]
        r_num = rec["analytic"]["analytic_rank"]
        if inv["conductor"] < 20000 or r_num > 6 or rec["label"] in done:
            continue
        t = time.perf_counter()
        res = certify_order_rigorous(inv["ainvs"], kmax=max(4, r_num), prec=128, label=rec["label"],
                                     use_modular_symbols=False, digits=40)
        res.pop("assumptions", None)
        res.pop("c_conditioned_note", None)
        res["numerical_analytic_rank_pari"] = r_num
        res["seconds"] = round(time.perf_counter() - t, 2)
        append(path, res)
        log(f"beyond {rec['label']}: N={res['conductor']} candidates={res['order_candidates']} t={res['seconds']}s")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    pari = get_pari()
    try:
        pari.allocatemem(STACK, silent=True)
    except TypeError:
        pari.allocatemem(STACK)
    pari.set_real_precision(140)
    log(f"=== run_batch {which} pid={os.getpid()} ===")
    if which in ("validation", "all"):
        phase_validation(pari)
    if which in ("timing", "all"):
        phase_timing(pari)
    if which in ("certify", "all"):
        phase_certify(pari)
    if which in ("rank3", "all"):
        phase_rank3(pari)
    if which in ("crosscheck", "all"):
        phase_crosscheck(pari)
    if which in ("beyond", "all"):
        phase_beyond(pari)
    log(f"=== run_batch {which} finished ===")


if __name__ == "__main__":
    main()
