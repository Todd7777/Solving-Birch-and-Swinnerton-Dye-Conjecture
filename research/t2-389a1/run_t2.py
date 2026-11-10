#!/usr/bin/env python3
"""Target T2 for 389a1: certify corank Sel_{p^oo}(E/Q) <= 2, hence Sha(E/Q)[p^oo] finite, prime by prime.

For each odd prime p of good ordinary reduction this computes the level-2 Riemann sum of the exact
p-adic L-function of bsd/kato_bound.py inside GP (same mathematics, one GP loop instead of Python),
giving c_0, c_1, c_2 modulo p^3 with c_j(L) = c_j(L^(2)) mod p for j = 1, 2 and c_0 exact.
c_2 != 0 mod p certifies ord_T L_p <= 2. When c_2 = 0 mod p ("regulator-anomalous" p), the
Python implementation at level 3 (or 4) is tried within a budget.

Usage:  python3 research/t2-389a1/run_t2.py --pmin 3 --pmax 3000 --shard 0 --nshards 1 [--tag NAME]
Output: research/t2-389a1/data/t2_shard<i>.jsonl, or t2_<NAME>_shard<i>.jsonl with --tag (resumable).
The p < 3000 run used no tag; the extension 3000 < p < 5000 used --tag ext5000 with 4 shards.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from bsd import kato_bound as kb  # noqa: E402
from bsd.pari_session import get_pari  # noqa: E402

AINVS = [0, 1, 1, -2, 0]           # 389a1
GP_LEVEL2 = """(M,xp,p,gnum,gden,ai2,ai3)->
{
  my(q=p^2, P=p^3, V1, gk, b, s, c0=0, c1=0, c2=0, tch);
  V1 = vector(p-1, a, mseval(M,xp,[oo,a/p])*gden/gnum);
  tch = vector(p-1, t, lift(Mod(t,q)^p));
  gk = Mod(1,q);
  for(k=0, p-1,
    s = 0;
    for(i=1, p-1,
      b = lift(tch[i]*gk);
      s += ai2*(mseval(M,xp,[oo,b/q])*gden/gnum) - ai3*V1[b%p]);
    s = s % P;
    c0 += s; c1 += s*k; c2 += s*(k*(k-1)/2);
    gk *= (1+p));
  [c0 % P, c1 % P, c2 % P]
}"""


def level2_gp(sym, p, ap):
    pari = get_pari()
    P = p ** 3
    alpha = kb.unit_root(ap, p, 3)
    ainv = pow(alpha, -1, P)
    f = pari(GP_LEVEL2)
    c = f(sym.M, sym.xp, p, sym.g.numerator, sym.g.denominator, pow(ainv, 2, P), pow(ainv, 3, P))
    return [int(x) for x in c]


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--pmin", type=int, default=3)
    ap_.add_argument("--pmax", type=int, default=3000)
    ap_.add_argument("--shard", type=int, default=0)
    ap_.add_argument("--nshards", type=int, default=1)
    ap_.add_argument("--level3-max-p", type=int, default=250)
    ap_.add_argument("--tag", default="", help="name for a separate set of output files")
    args = ap_.parse_args()
    pari = get_pari(stack=2 ** 31)
    E = pari.ellinit(AINVS)
    N = int(pari.ellglobalred(E)[0])
    sym = kb.PlusSymbol(E, N)
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    out = os.path.join(HERE, "data", f"t2_{args.tag + '_' if args.tag else ''}shard{args.shard}.jsonl")
    done = set()
    if os.path.exists(out):
        for line in open(out):
            done.add(json.loads(line)["p"])
    primes = [int(q) for q in pari.primes([args.pmin, args.pmax])]
    primes = [p for i, p in enumerate(primes) if i % args.nshards == args.shard and p not in done]
    # largest first balances the shards' running times
    for p in sorted(primes, reverse=True):
        t0 = time.time()
        rec = dict(p=p)
        if N % p == 0:
            rec.update(kind="bad", a_p=int(pari.ellap(E, p)), note="multiplicative; not treated here")
        else:
            ap = int(pari.ellap(E, p))
            rec["a_p"] = ap
            if ap % p == 0:
                rec.update(kind="supersingular", note="needs Kobayashi/Pollack +- theory; not treated here")
            else:
                rec["kind"] = "ordinary"
                rec["mazur_anomalous"] = (ap - 1) % p == 0
                c = level2_gp(sym, p, ap)
                rec["level2_coeffs_mod_p3"] = c
                rec["c0_exact_zero"] = c[0] == 0
                rec["c1_zero_mod_p"] = c[1] % p == 0
                rec["c2_mod_p"] = c[2] % p
                if c[0] != 0 or c[1] % p != 0:
                    rec["inconsistent"] = True          # Kato + rank 2 force c_0 = 0 and c_1 = 0
                if c[2] % p:
                    rec.update(ord_upper_bound=2, level=2, certified=True)
                else:
                    rec["regulator_anomalous_level2"] = True
                    rec["certified"] = False
                    if p <= args.level3_max_p:
                        for n in (3, 4):
                            if p ** n > 3 * 10 ** 7:
                                break
                            rc = kb.riemann_coefficients(sym, p, ap, n, 2)
                            v2 = kb.vp(rc["coeffs"][2], p)
                            e2 = rc["precision"][2]
                            rec[f"level{n}"] = dict(c2_mod=rc["coeffs"][2], modulus=f"{p}^{rc['P']}", v2=v2, e2=e2)
                            if v2 is not None and v2 < e2:
                                rec.update(ord_upper_bound=2, level=n, certified=True, c2_valuation=v2)
                                break
        rec["seconds"] = round(time.time() - t0, 2)
        with open(out, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        print(time.strftime("%H:%M:%S"), json.dumps(rec)[:200], flush=True)


if __name__ == "__main__":
    main()
