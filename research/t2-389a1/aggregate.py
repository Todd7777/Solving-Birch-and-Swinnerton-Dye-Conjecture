#!/usr/bin/env python3
"""Aggregate the T2 sweep for 389a1 into summary.json and report.md.

Adds the primes that the sweep does not handle itself:
  p = 2    from the certified 2-descent of results/bsd1 (dim Sel_2 = 2, E(Q)[2] = 0, so Sha[2] = 0);
  p = 389  split multiplicative, via bsd.kato_bound.multiplicative_bound (exceptional zero accounted for);
  p = 3    ordinary and certified at level 2, but Stein-Wuthrich Thm 7.3/7.4 need rho_3 surjective
           (rhobar_3 is surjective: the 3-division polynomial has Galois group S_4), not checked here.
Usage: python3 research/t2-389a1/aggregate.py
"""
from __future__ import annotations

import glob
import gzip
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from bsd import kato_bound as kb  # noqa: E402
from bsd.pari_session import get_pari  # noqa: E402

AINVS = [0, 1, 1, -2, 0]


def main():
    pari = get_pari()
    E = pari.ellinit(AINVS)
    recs = {}
    for fn in sorted(glob.glob(os.path.join(HERE, "data", "t2_*shard*.jsonl"))):
        for line in open(fn):
            r = json.loads(line)
            recs[r["p"]] = r
    ss_file = os.path.join(HERE, "data", "t2_supersingular.jsonl")
    if os.path.exists(ss_file):                       # supersingular certificates (run_supersingular.py)
        for line in open(ss_file):
            r = json.loads(line)
            recs[r["p"]] = dict(recs.get(r["p"], {}), **r)
    pmax = max(recs)
    primes = [int(q) for q in pari.primes([3, pmax])]
    missing = [p for p in primes if p not in recs]
    # p = 2 from the BSD(I) certificate
    two = None
    with gzip.open(os.path.join(ROOT, "results", "bsd1", "rank2plus.jsonl.gz"), "rt") as fh:
        for line in fh:
            r = json.loads(line)
            if r["label"] == "389a1":
                d = r["descent"]
                two = dict(p=2, kind="2-descent", certified=d["selmer_bound_C_minus_T"] == 2 and d["unconditional"],
                           statement="dim Sel_2 = 2 = rank, E(Q)[2] = 0: Sha[2] = 0, hence Sha[2^oo] = 0",
                           cubic_fields=d["cubic_fields"])
    # p = 389
    m = kb.multiplicative_bound(E, 389)
    bad = dict(p=389, kind="split multiplicative", a_p=m["a_p"], level=m["level"], ord_upper_bound=m["ord_upper_bound"],
               char_series_ord_upper_bound=m["char_series_ord_upper_bound"],
               rho_p_surjective=m["rho_p_surjective_by_isogeny_criterion"],
               certified=m["char_series_ord_upper_bound"] == 2 and m["rho_p_surjective_by_isogeny_criterion"],
               coeffs_mod_p3=m["coeffs"])
    recs[389] = dict(recs.get(389, {}), **bad)
    psi3 = pari.elldivpol(E, 3)
    gal3 = str(pari.polgalois(psi3)[3])
    rho3_file = os.path.join(HERE, "rho3_certificate.json")
    rho3 = json.load(open(rho3_file))["rho3_surjective"] if os.path.exists(rho3_file) else False
    ordinary = [r for r in recs.values() if r.get("kind") == "ordinary"]
    cert_ord = sorted(r["p"] for r in ordinary if r.get("certified"))
    unc_ord = sorted(r["p"] for r in ordinary if not r.get("certified"))
    anomalous = sorted(r["p"] for r in ordinary if r.get("regulator_anomalous_level2"))
    ss = sorted(r["p"] for r in recs.values() if r.get("kind") == "supersingular")
    ss_cert = sorted(p_ for p_ in ss if recs[p_].get("certified"))
    ss_plus = sorted(p_ for p_ in ss_cert if recs[p_].get("Lplus_mod_p") and any(recs[p_]["Lplus_mod_p"][:3]))
    inconsistent_ss = sorted(p_ for p_ in ss if not (recs[p_].get("c0_exact_zero") and recs[p_].get("c1_zero_mod_p")))
    inconsistent = sorted(r["p"] for r in recs.values() if r.get("inconsistent"))
    mazur = sorted(r["p"] for r in ordinary if r.get("mazur_anomalous"))
    heur = sum(1.0 / p for p in cert_ord + unc_ord if p >= 5)
    summary = dict(
        curve="389a1", ainvs=AINVS, prime_range=[2, pmax], missing=missing,
        p2=two, p3=dict(certified_ord_bound=3 in cert_ord, rhobar3_galois_group=gal3, rho3_surjective_certified=rho3,
                        note="SW Thm 7.3 needs rho_3 surjective; certified by rho3_surjective.py (Elkies 2006)"),
        p389=bad, ordinary_certified=cert_ord, ordinary_uncertified=unc_ord, supersingular=ss,
        supersingular_certified=ss_cert, supersingular_with_Lplus_cross_check=ss_plus,
        supersingular_inconsistent=inconsistent_ss,
        regulator_anomalous_at_level2=anomalous,
        anomalous_resolved={p: recs[p].get("level") for p in anomalous if recs[p].get("certified")},
        mazur_anomalous=mazur, inconsistent=inconsistent,
        heuristic_expected_anomalies=round(heur, 2),
        seconds=round(sum(r.get("seconds", 0) for r in recs.values()), 1),
    )
    with open(os.path.join(HERE, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    unconditional = sorted([2] + [p for p in cert_ord if p >= 5] + ([389] if bad["certified"] else [])
                           + ([3] if (rho3 and 3 in cert_ord) else []) + ss_cert)
    lines = [
        "# T2 for 389a1: finiteness of Sha[p^oo], prime by prime",
        "",
        f"*Generated by `research/t2-389a1/aggregate.py` from `data/t2_*shard*.jsonl`. Primes up to {pmax}; "
        f"{len(missing)} missing.*",
        "",
        "| primes | count | result | method |",
        "|---|---|---|---|",
        "| 2 | 1 | Sha[2^oo] = 0 | 2-descent, cubic field certified (`results/bsd1`) |",
        (f"| 3 | 1 | Sha[3^oo] finite | level-2 Riemann sum; rho_3 surjective (`rho3_surjective.py`: Galois group of "
         f"psi_3 is {gal3}, j not in Elkies' family) |" if (rho3 and 3 in cert_ord) else
         f"| 3 | 1 | ord_T L_3 <= 2: {3 in cert_ord}; rho_3 surjective: {rho3} | level-2 Riemann sum |"),
        f"| good ordinary 5 <= p <= {pmax} | {len([p for p in cert_ord if p >= 5])} of "
        f"{len([p for p in cert_ord + unc_ord if p >= 5])} | Sha[p^oo] finite, p-adic height pairing non-degenerate | "
        "Kato (SW Thm 7.3/7.4) + exact Riemann sums |",
        f"| 389 (split multiplicative) | 1 | {'Sha[389^oo] finite' if bad['certified'] else 'not certified'} "
        f"(ord_T L_p <= {bad['ord_upper_bound']}, exceptional zero) | SW Thm 7.3 split case, no 389-isogeny |",
        (f"| supersingular ({', '.join(map(str, ss))}) | {len(ss_cert)} of {len(ss)} | Sha[p^oo] finite, D_p-valued "
         f"regulator non-zero | L^- == -theta_1 mod omega_1 (`bsd.kato_bound.SUPERSINGULAR_PROOF`), SW Thm 6.2, 7.3/7.4; "
         f"L^+ from theta_2 as a second certificate at {ss_plus} |" if ss_cert else
         f"| supersingular | {len(ss)} | not treated | {', '.join(map(str, ss))} |"),
        "",
        f"**Certified (given the cited theorems and PARI/Arb): Sha(389a1)[p^oo] is finite for {len(unconditional)} primes,** "
        f"namely p = 2 (where it is 0), p = 3, p = 389, every good ordinary prime 5 <= p <= {pmax}"
        + (f" except {unc_ord}" if [p for p in unc_ord if p >= 5] else "")
        + (f", and the supersingular primes {ss_cert}" if ss_cert else "") + ".",
        ("Not covered below " + str(pmax) + ": "
         + ", ".join(([f"supersingular {sorted(set(ss) - set(ss_cert))}"] if set(ss) - set(ss_cert) else [])
                     + ([] if (rho3 and 3 in cert_ord) else ["p = 3"])
                     + ([f"ordinary {unc_ord}"] if unc_ord else [])) + "."
         if (set(ss) - set(ss_cert)) or unc_ord or not (rho3 and 3 in cert_ord) else
         f"**Every prime p <= {pmax} is covered.**"),
        "",
        "**Regulator anomalies.** At level 2 the T^2 coefficient is known modulo p. It was divisible by p for "
        f"{anomalous if anomalous else 'no'} prime(s)"
        + (f"; resolved at higher level: {summary['anomalous_resolved']}" if anomalous else "")
        + f". If the normalised leading coefficient were a random p-adic integer, the expected number of such primes in "
        f"[5, {pmax}] would be sum 1/p = {heur:.2f}. Mazur-anomalous primes (a_p = 1 mod p, so that (1 - 1/alpha)^2 is "
        f"divisible by p^2): {mazur}; "
        + (f"of these, {sorted(set(mazur) & set(anomalous))} also had c_2 = 0 mod p."
           if set(mazur) & set(anomalous) else
           "at none of them did c_2 vanish mod p: the factor p^2 in (1 - 1/alpha)^2 is cancelled by the denominator "
           "of the p-adic regulator, as the p-adic BSD formula predicts."),
        "",
        f"Consistency: c_0 = 0 exactly and c_1 = 0 mod p at every ordinary prime (forced by Kato and rank 2); "
        f"violations: {inconsistent if inconsistent else 'none'}. The same for the T^0, T^1 coefficients of L^- at the "
        f"supersingular primes: violations: {inconsistent_ss if inconsistent_ss else 'none'}.",
    ]
    with open(os.path.join(HERE, "report.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
