#!/usr/bin/env python3
"""T2 for 389a1 at the supersingular primes p < 3000 (107, 599, 1049, 2957; all have a_p = 0). There are no
others below 5000; the ordinary run classifies every prime and would list them.

For each p, bsd.kato_bound.supersingular_bound certifies ord_T L_p <= 2 for the D_p-valued p-adic L-function from
[T^j] L^- == -[T^j] theta_1 (mod p) (proof: bsd.kato_bound.SUPERSINGULAR_PROOF; cost p^2 symbol values). With
--plus-max-p, primes up to that bound also get the independent certificate from L^+ == theta_2 / Phi_p(1+T) (mod p),
which costs p^3 symbol values.

Usage: python3 research/t2-389a1/run_supersingular.py [--plus-max-p 600]
       python3 research/t2-389a1/run_supersingular.py --pmin 5000 --pmax 16231 --tag ext16231
Output: research/t2-389a1/data/t2_supersingular.jsonl, or t2_supersingular_<tag>.jsonl with --tag (one record per
prime; reruns replace the file).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from bsd import kato_bound as kb  # noqa: E402
from bsd.pari_session import get_pari  # noqa: E402

AINVS = [0, 1, 1, -2, 0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmin", type=int, default=5)
    ap.add_argument("--pmax", type=int, default=3000)
    ap.add_argument("--tag", default="")
    ap.add_argument("--plus-max-p", type=int, default=600)
    args = ap.parse_args()
    pari = get_pari(stack=2 ** 31)
    E = pari.ellinit(AINVS)
    sym = kb.PlusSymbol(E)
    N = sym.N
    primes = [int(q) for q in pari.primes([max(5, args.pmin), args.pmax]) if N % int(q) and int(pari.ellap(E, int(q))) == 0]
    out = os.path.join(HERE, "data", f"t2_supersingular{'_' + args.tag if args.tag else ''}.jsonl")
    recs = []
    for p in primes:
        t0 = time.time()
        res = kb.supersingular_bound(E, p, sym=sym, target=2, cross_check_plus=p <= args.plus_max_p,
                                     max_evaluations=max(3 * 10 ** 8, 2 * p * p))
        rec = dict(p=p, kind="supersingular", a_p=res["a_p"], ord_upper_bound=res["ord_upper_bound"],
                   certified=res["ord_upper_bound"] is not None and res["ord_upper_bound"] <= 2,
                   theta1_coeffs=res["theta1_coeffs"], Lminus_mod_p=res["Lminus_mod_p"],
                   Lplus_mod_p=res.get("Lplus_mod_p"),
                   c0_exact_zero=res["theta1_coeffs"][0] == 0,
                   c1_zero_mod_p=res["Lminus_mod_p"][1] == 0,
                   attempts=[{k: a[k] for k in ("sign", "k", "level", "precision", "coeffs_mod")}
                             for a in res["attempts"]],
                   seconds=round(time.time() - t0, 1))
        recs.append(rec)
        print(time.strftime("%H:%M:%S"), json.dumps({k: rec[k] for k in ("p", "ord_upper_bound", "Lminus_mod_p",
                                                                          "Lplus_mod_p", "seconds")}), flush=True)
    with open(out, "w") as fh:
        for rec in recs:
            fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
