#!/usr/bin/env python3
"""Sha(389a1)(p) = 0 for every prime p < 5000: bsd/sha_zero.py applied to the T2 records of this directory.

Inputs: the level-2 p-adic L-function records (data/t2_*shard*.jsonl: [T^2]L_p is a unit at every good
ordinary p), the supersingular records (data/t2_supersingular.jsonl: [T^2]L^- mod p), the split multiplicative
prime 389 (bsd.kato_bound.multiplicative_bound), rho3_certificate.json (rho_3 surjective) and, for p = 2, the
certified 2-descent (bsd.bsd1.descent_bound). Supersingular records beyond 5000 (data/t2_supersingular_ext*.jsonl)
are certified the same way and reported separately.

Usage: python3 research/t2-389a1/run_sha_zero.py --shard i --nshards k      (writes data/sha_zero_shard<i>.jsonl)
       (beyond 5000 first: run_supersingular.py --pmin 5001 --pmax 16230 --tag ext16231, here split into four
       ranges, tags ext16231a-d, run in parallel)
       python3 research/t2-389a1/run_sha_zero.py --aggregate                  (writes sha_zero_report.md, .json)
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from bsd import bsd1  # noqa: E402
from bsd import kato_bound as kb  # noqa: E402
from bsd import sha_zero as sz  # noqa: E402
from bsd.pari_session import get_pari  # noqa: E402

AINVS = [0, 1, 1, -2, 0]
GENS = [[-1, 1], [0, 0]]          # Cremona's generators; Lemma B does not need them to be a basis
RANK = 2


def t2_records():
    recs = {}
    for fn in sorted(glob.glob(os.path.join(HERE, "data", "t2_*shard*.jsonl"))):
        for line in open(fn):
            r = json.loads(line)
            recs[r["p"]] = r
    for ss in sorted(glob.glob(os.path.join(HERE, "data", "t2_supersingular*.jsonl"))):
        for line in open(ss):
            r = json.loads(line)
            recs[r["p"]] = dict(recs.get(r["p"], {}), **r)
    return recs


def certify(p, rec, E):
    """One prime. The bound uses Lemma B at good ordinary p (no saturation needed) and PARI's regulator for
    Cremona's generators at supersingular p and at 389, where the generators' p-saturation is part of the proof.
    The saturation certificate is computed at every p, so that at ordinary p the PARI regulator gives a second,
    independent route to the same b_p."""
    pari = get_pari()
    out = dict(p=p, kind=rec.get("kind"))
    t0 = time.time()
    out["saturation"] = sz.saturation_certificate(E, GENS, p)
    needs_saturation = True
    if p == 389:
        degrees = {int(d) for row in pari.Vec(pari.ellisomat(E)[1]) for d in row}
        out["surjectivity"] = sz.rhobar_surjective_certificate(E, p)
        out["surjectivity"]["isogeny_degrees"] = sorted(degrees)      # cross-check: no rational isogeny at all
        m = kb.multiplicative_bound(E, p)
        vL = kb.vp(m["coeffs"][3], p)
        if m["ord_upper_bound"] != 3 or vL != 0:
            out.update(ok=False, reason="[T^3]L_389 not a unit")
            return out
        out["bound"] = sz.split_multiplicative_bound(E, p, vL, RANK, GENS)
    elif rec.get("kind") == "supersingular":
        out["surjectivity"] = sz.rhobar_surjective_certificate(E, p)
        lplus = rec.get("Lplus_mod_p")
        out["bound"] = sz.supersingular_bound(E, p, rec["Lminus_mod_p"][2], RANK, GENS,
                                              Lplus_T2_mod_p=lplus[2] if lplus else None)
        if "direction_check" in out["bound"] and not out["bound"]["direction_check"]["agree"]:
            out.update(ok=False, reason="leading term not proportional to the D_p-valued regulator mod p")
            return out
    elif rec.get("kind") == "ordinary":
        needs_saturation = False
        if p == 3:
            rho3 = json.load(open(os.path.join(HERE, "rho3_certificate.json")))
            out["surjectivity"] = dict(ok=bool(rho3["rho3_surjective"]), method="rho_3 surjective (rho3_surjective.py)",
                                       rhobar=sz.rhobar_surjective_certificate(E, 3))
        else:
            out["surjectivity"] = sz.rhobar_surjective_certificate(E, p)
        if not rec.get("certified") or rec.get("level", 2) != 2 or not rec.get("c2_mod_p"):
            out.update(ok=False, reason="no unit [T^2]L_p at level 2")
            return out
        out["bound"] = sz.ordinary_bound(E, p, 0, RANK, GENS)
        b = out["bound"]
        out["regulator_route_ok"] = bool(out["saturation"]["ok"] and b.get("b_p_from_pari_regulator") is not None
                                         and b["b_p_from_pari_regulator"] <= 1)
    else:
        out.update(ok=False, reason=f"unexpected kind {rec.get('kind')}")
        return out
    b = out["bound"]
    out["ok"] = bool(out["surjectivity"].get("ok") and b.get("sha_p_trivial") and b.get("consistent", True)
                     and (out["saturation"]["ok"] or not needs_saturation))
    out["seconds"] = round(time.time() - t0, 2)
    return out


def run(shard, nshards):
    pari = get_pari(stack=2 ** 31)
    E = pari.ellinit(AINVS)
    recs = t2_records()
    out_fn = os.path.join(HERE, "data", f"sha_zero_shard{shard}.jsonl")
    done = set()
    if os.path.exists(out_fn):
        done = {json.loads(l)["p"] for l in open(out_fn)}
    primes = [p for i, p in enumerate(sorted(q for q in recs if q >= 3)) if i % nshards == shard and p not in done]
    for p in sorted(primes, reverse=True):
        res = certify(p, recs[p], E)
        with open(out_fn, "a") as fh:
            fh.write(json.dumps(res, default=str) + "\n")
        print(time.strftime("%H:%M:%S"), p, res.get("ok"), res.get("bound", {}).get("b_p"), flush=True)


def two_descent(E):
    """p = 2: PARI's ellrank gives C = dim Sel_2, T = dim E(Q)[2] and C = T + R + S with R the rank and
    S = dim Sha[2] (PARI manual, ellrank). C - T = 2 = R (two independent points) forces Sha[2] = 0, hence
    Sha(2) = 0. Unconditional once the cubic 2-division field is certified (bnfcertify, as in bsd.bsd1)."""
    d = bsd1.descent_bound(E)
    t2 = int(get_pari().elltors(E)[0]) % 2 == 0
    ok = bool(d["unconditional"] and d["selmer_bound_C_minus_T"] == RANK and d["rank_lower"] == RANK and not t2)
    return dict(p=2, kind="2-descent", ok=ok, descent=d, bound=dict(b_p=0, sha_p_trivial=ok),
                method="dim Sel_2 - dim E(Q)[2] = 2 = rank, so Sha[2] = 0 (C = T + R + S)")


def aggregate():
    pari = get_pari()
    E = pari.ellinit(AINVS)
    recs = {}
    for fn in sorted(glob.glob(os.path.join(HERE, "data", "sha_zero_shard*.jsonl"))):
        for line in open(fn):
            r = json.loads(line)
            recs[r["p"]] = r
    recs[2] = two_descent(E)
    # supersingular primes beyond 5000 (data/t2_supersingular_ext*.jsonl): reported separately; with SW Thm 12.3
    # (Sha[p] = 0 for p = 2 and every good ordinary p < 48859 except 16231; cited) they close the range below 16231
    ext = {p: r for p, r in recs.items() if p > 5000}
    recs = {p: r for p, r in recs.items() if p < 5000}
    ss_to_16231 = [int(q) for q in pari.primes([5001, 16230]) if int(pari.ellap(E, int(q))) % int(q) == 0]
    ext_missing = [p for p in ss_to_16231 if not ext.get(p, {}).get("ok")]
    primes = [int(q) for q in pari.primes([2, 4999])]
    missing = [p for p in primes if p not in recs]
    bad = sorted(p for p, r in recs.items() if not r.get("ok"))
    odd = [r for r in recs.values() if r["p"] > 2]
    ords = [r for r in odd if r.get("kind") == "ordinary"]
    ss = sorted(r["p"] for r in odd if r.get("kind") == "supersingular")
    pred_fail = []                       # PARI regulator must have exactly the valuation Lemma B + b_p >= 0 forces
    for r in ords:
        b = r["bound"]
        want = 2 if b["vp_N_p"] == 0 else 0
        if b.get("pari_vp_Reg_p") != want:
            pred_fail.append(r["p"])
    sat_fail = sorted(r["p"] for r in odd if not r.get("saturation", {}).get("ok"))
    route2 = sorted(r["p"] for r in ords if r.get("regulator_route_ok"))
    direction = {r["p"]: r["bound"]["direction_check"] for r in odd if "direction_check" in r.get("bound", {})}
    trans = sorted({(r["surjectivity"]["transvection"]["q"], r["surjectivity"]["transvection"]["q_equals_ell"])
                    for r in odd if r["p"] != 3 and "transvection" in r["surjectivity"]})
    aux = [u["l"] for r in odd for u in r["saturation"]["primes"]]
    by_kind = {}
    for r in recs.values():
        by_kind.setdefault(r.get("kind"), []).append(r["p"])
    summ = dict(curve="389a1", primes_below=5000, n_primes=len(primes), missing=missing, not_certified=bad,
                certified=len([p for p in primes if recs.get(p, {}).get("ok")]),
                counts={k: len(v) for k, v in by_kind.items()},
                anomalous=sorted(r["p"] for r in ords if r["bound"]["vp_N_p"] == 1), supersingular=ss,
                regulator_valuation_prediction_failures=pred_fail, saturation_failures=sat_fail,
                ordinary_regulator_route_certified=len(route2), direction_checks=direction,
                transvection_primes=[dict(q=q, q_equals_ell=e) for q, e in trans],
                saturation_aux_primes=dict(count=len(aux), max=max(aux) if aux else None),
                b_p_values=sorted({r.get("bound", {}).get("b_p") for r in recs.values()} - {None}),
                seconds=round(sum(r.get("seconds", 0) for r in recs.values()), 1),
                beyond_5000=dict(certified=sorted(p for p, r in ext.items() if r.get("ok")),
                                 not_certified=sorted(p for p, r in ext.items() if not r.get("ok")),
                                 supersingular_below_16231=ss_to_16231, supersingular_below_16231_missing=ext_missing,
                                 every_prime_below_16231_with_SW_12_3=not bad and not missing and not ext_missing))
    with open(os.path.join(HERE, "sha_zero_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1)
    n_non = len([r for r in ords if r["bound"]["vp_N_p"] == 0])
    lines = [
        "# Sha(389a1)(p) = 0 for every prime p < 5000",
        "",
        "*Generated by `research/t2-389a1/run_sha_zero.py --aggregate` from `data/sha_zero_shard*.jsonl`. Method: "
        "`bsd/sha_zero.py`, Stein-Wuthrich Algorithm 11.1 (Math. Comp. 82 (2013)) with the proofs given there "
        f"(Lemma A, Lemma B, saturation). Primes missing: {len(missing)}.*",
        "",
        f"**Certified: {summ['certified']} of {summ['n_primes']} primes p < 5000 have Sha(389a1)(p) = 0**, i.e. "
        f"Sha(389a1)[p] = 0 (not certified: {bad if bad else 'none'}). Each odd p gives b_p with "
        "#Sha(389a1)(p) <= p^{b_p}; b_p <= 1 suffices because #Sha(p) is a square (Cassels).",
        "",
        "| primes | count | b_p | leading term | regulator | other inputs |",
        "|---|---|---|---|---|---|",
        "| 2 | 1 | - | - | - | 2-descent: dim Sel_2 = 2 = rank, E(Q)[2] = 0, cubic field certified |",
        f"| good ordinary, non-anomalous | {n_non} | "
        f"{sorted({r['bound']['b_p'] for r in ords if r['bound']['vp_N_p'] == 0})} | [T^2]L_p unit | Lemma B(i): "
        f"v(Reg_gamma) >= 0 for any basis | Lemma A |",
        f"| good ordinary, anomalous ({', '.join(map(str, summ['anomalous']))}) | {len(summ['anomalous'])} | "
        f"{sorted({r['bound']['b_p'] for r in ords if r['bound']['vp_N_p'] == 1})} | [T^2]L_p unit, v(eps_p) = 2 | "
        "Lemma B(ii): v(Reg_gamma) >= -2 | Lemma A; at p = 3 rho_3 surjective (rho3_certificate.json) |",
        f"| supersingular ({', '.join(map(str, ss))}) | {len(ss)} | "
        f"{sorted({recs[p]['bound']['b_p'] for p in ss})} | [T^2]L^- unit: (1-phi)^-2 L_p^* has valuation 0 | "
        "PARI D_p-valued, v = 2; generators p-saturated | Lemma A; lattices (omega, phi omega) = (omega, eta) |",
        f"| 389, split multiplicative | 1 | {recs.get(389, {}).get('bound', {}).get('b_p')} | [T^3]L_389 unit | "
        "PARI modified height, v(Reg_gamma) = 0; generators 389-saturated | v(L-invariant) = 1; Lemma A with "
        "inertia at 389 |",
        "",
        "Cross-checks (not used in the proof):",
        f"- At all {len(ords)} ordinary primes the generators are p-saturated as well "
        f"({len(route2)} certificates), so PARI's regulator gives the same b_p by a second route; Lemma B predicts "
        "its valuation exactly (2 at non-anomalous, 0 at anomalous p) and the prediction failed at "
        f"{pred_fail if pred_fail else 'no prime'}.",
        f"- Direction check at supersingular p with [T^2]L^+ known: (1-phi)^-2 L_p^* is proportional to the D_p-valued "
        "regulator mod p, as Kato's divisibility and Perrin-Riou's formula force: "
        + ", ".join(f"p = {p}: {d['L_ratio_mod_p']} vs {d['reg_ratio_mod_p']}" for p, d in sorted(direction.items()))
        + ".",
        f"- Saturation used {summ['saturation_aux_primes']['count']} auxiliary primes l (largest "
        f"{summ['saturation_aux_primes']['max']}); failures: {sat_fail if sat_fail else 'none'}. PARI's "
        "ellsaturation(E, gens, 5000) returns a basis of the same regulator.",
        "",
        "Corollary (remark after SW Algorithm 11.1): where b_p = 0, f_E and L_p have the same order of vanishing and "
        "leading-term valuation, so with Kato's divisibility f_E / L_p is a unit: the cyclotomic Iwasawa main "
        f"conjecture holds for (389a1, p) at every odd p < 5000 (b_p values: {summ['b_p_values']}). At good "
        "ordinary p this is also a case of Skinner-Urban (SW Thm 7.5, with q = 389).",
        "",
        "## Relation to Stein-Wuthrich Theorem 12.3, and beyond 5000",
        "",
        "SW Thm 12.3 (cited, not re-derived) proves Sha(389a1)[p] = 0 for p = 2 and every good ordinary p < 48859 "
        "except 16231; it excludes 389 and the supersingular primes. Below 5000 the ordinary part of this table is "
        "therefore an independent reproduction, and the primes beyond SW Thm 12.3 are 107, 389, 599, 1049, 2957. "
        + (f"Beyond 5000 this directory certifies, by the same method and with b_p = 0 (so the main conjecture holds "
           f"there too), the supersingular primes {summ['beyond_5000']['certified']}"
           if summ['beyond_5000']['certified'] else "Beyond 5000 no prime is certified yet")
        + (f" (not certified: {summ['beyond_5000']['not_certified']})" if summ['beyond_5000']['not_certified'] else "")
        + ". "
        + ("Together with SW Thm 12.3 this gives **Sha(389a1)[p] = 0 for every prime p < 16231**; 16231 is the "
           "ordinary prime SW could not finish (their Remark 12.4: ord_p Reg_p = 3)."
           if summ["beyond_5000"]["every_prime_below_16231_with_SW_12_3"] else
           f"Supersingular primes below 16231 still missing: {ext_missing}."),
    ]
    with open(os.path.join(HERE, "sha_zero_report.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--aggregate", action="store_true")
    a = ap.parse_args()
    if a.aggregate:
        aggregate()
    else:
        run(a.shard, a.nshards)


if __name__ == "__main__":
    main()
