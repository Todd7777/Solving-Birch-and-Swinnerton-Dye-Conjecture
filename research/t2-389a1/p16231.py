#!/usr/bin/env python3
"""p = 16231 for 389a1: the one prime below 48859 where Sha(389a1)[p] = 0 is not proved (Stein-Wuthrich Remark 12.4).

What this script records (output: data/p16231.json):
1. The algebraic side. PARI's p-adic regulator of Cremona's generators has valuation 3, and the generators are
   p-saturated (bsd.sha_zero.saturation_certificate). p is good ordinary and not anomalous, so ord_p Reg_gamma = 1
   for E(Q).
2. The analytic side at level 2: c_0..c_4 of the level-2 Riemann sum, correct mod p for 1 <= j < p. This is the GP
   loop of run_t2.py extended to T^4 and halved by evenness, with the full sum recomputed at k = 0 as a guard.
   Kato's divisibility and SW Thm 6.1 force c_2 == 0 mod p, since ord_p Reg_gamma = 1. Numerically,
   c_3 == -(1 + log_gamma<389>/2) c_2 mod p holds wherever it was tested, consistent with the functional equation
   of L_p. So c_3 == 0 too, and level-2 data can show at best ord_T L_p <= 4. That does not even give finiteness,
   because the p-parity theorem allows Selmer corank 2 or 4. Deciding Sha(p) needs [T^2]L_p modulo p^2.
3. With --oms: the memory of PARI's overconvergent-symbol route (mspadicinit / mspadicmoments) at small p for this
   level, to show why it cannot reach p = 16231.

Usage: python3 research/t2-389a1/p16231.py [--oms]      (about 8 CPU-minutes; --oms adds about a minute)
       python3 research/t2-389a1/p16231.py --oms-only       (only the overconvergent memory table)
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from bsd import kato_bound as kb  # noqa: E402
from bsd import sha_zero as sz  # noqa: E402
from bsd.pari_session import get_pari  # noqa: E402

AINVS, GENS, P0 = [0, 1, 1, -2, 0], [[-1, 1], [0, 0]], 16231
GP_LEVEL2_J4 = """(M,xp,p,gnum,gden,ai2,ai3)->
{
  my(q=p^2, P=p^3, V1, gk, b, s, c=vector(5), tch, h=(p-1)/2, full=0);
  V1 = vector(p-1, a, mseval(M,xp,[oo,a/p])*gden/gnum);
  tch = vector(p-1, t, lift(Mod(t,q)^p));
  for(i=1, p-1, b = tch[i]; full += ai2*(mseval(M,xp,[oo,b/q])*gden/gnum) - ai3*V1[b%p]);
  gk = Mod(1,q);
  for(k=0, p-1,
    s = 0;
    for(i=1, h, b = lift(tch[i]*gk); s += ai2*(mseval(M,xp,[oo,b/q])*gden/gnum) - ai3*V1[b%p]);
    s = 2*s;
    if(k==0 && (s-full) % P, error("evenness guard failed"));
    s = s % P;
    for(j=0, 4, c[j+1] += s*binomial(k,j));
    gk *= (1+p));
  vector(5, j, c[j] % P)
}"""

OMS_PROBE = """
import sys, cypari2
pari = cypari2.Pari(); pari.allocatemem(10**8, silent=True); pari.default("parisizemax", 12 * 10**9)
E = pari.ellinit([0, 1, 1, -2, 0]); M, phi = pari.msfromell(E, 1)
Mp = pari.mspadicinit(M, int(sys.argv[1]), 3, 0); mu = pari.mspadicmoments(Mp, phi, 1); pari.mspadicseries(mu)
print([l for l in open("/proc/self/status") if l.startswith("VmHWM")][0].split()[1])   # peak RSS of this process, kB
"""


def level2_coeffs(sym, p, ap):
    pari = get_pari()
    P = p ** 3
    ainv = pow(kb.unit_root(ap, p, 3), -1, P)
    f = pari(GP_LEVEL2_J4)
    return [int(x) for x in f(sym.M, sym.xp, p, sym.g.numerator, sym.g.denominator, pow(ainv, 2, P), pow(ainv, 3, P))]


def fe_constant(p):
    """-(1 + log_gamma<389>/2) mod p, gamma = 1 + p."""
    pari = get_pari()
    x = pari(f"389 + O({p}^6)")
    ratio = pari.log(x / pari.teichmuller(x)) / pari.log(pari(f"1 + {p} + O({p}^6)"))
    return int(pari.lift(pari.Mod(pari.truncate(-(1 + ratio / 2)), p)))


def oms_rows():
    """Peak memory of mspadicinit + mspadicmoments + mspadicseries (3 digits) for 389a1, one fresh process per p."""
    rows = []
    for p in (5, 7, 11, 13, 17, 23):
        t = time.time()
        kb_ = subprocess.run([sys.executable, "-c", OMS_PROBE, str(p)], capture_output=True, text=True).stdout
        rows.append(dict(p=p, n=3, peak_rss_GB=round(int(kb_.split()[-1]) / 1e6, 2), seconds=round(time.time() - t, 1)))
    return rows


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--oms", action="store_true", help="also measure PARI's overconvergent route at small p")
    ap_.add_argument("--oms-only", action="store_true", help="only (re)measure the overconvergent route; keep the rest")
    ap_.add_argument("--check-primes", default="101,307,2969", help="primes where the c_3/c_2 relation is checked")
    a = ap_.parse_args()
    fn = os.path.join(HERE, "data", "p16231.json")
    if a.oms_only:
        out = json.load(open(fn))
        out["pari_overconvergent_memory"] = oms_rows()
        json.dump(out, open(fn, "w"), indent=1)
        print(json.dumps(out["pari_overconvergent_memory"], indent=1))
        return
    pari = get_pari(stack=2 ** 30)
    E = pari.ellinit(AINVS)
    sym = kb.PlusSymbol(E)
    out = dict(p=P0)
    ap = int(pari.ellap(E, P0))
    R = pari.ellpadicregulator(E, P0, 6, pari(GENS))
    out["algebraic"] = dict(a_p=ap, N_p=P0 + 1 - ap, anomalous=(P0 + 1 - ap) % P0 == 0, pari_regulator=str(R),
                            vp_Reg_p=int(pari.valuation(R, P0)),
                            saturated=sz.saturation_certificate(E, GENS, P0)["ok"],
                            rhobar_surjective=sz.rhobar_surjective_certificate(E, P0)["ok"])
    checks = []
    for p in [int(x) for x in a.check_primes.split(",") if x]:
        c = level2_coeffs(sym, p, int(pari.ellap(E, p)))
        checks.append(dict(p=p, c_mod_p=[x % p for x in c], c3_over_c2=c[3] * pow(c[2], -1, p) % p,
                           fe_constant=fe_constant(p)))
    out["c3_c2_relation_checks"] = checks
    t = time.time()
    c = level2_coeffs(sym, P0, ap)
    out["analytic_level2"] = dict(c_mod_p=[x % P0 for x in c], fe_constant=fe_constant(P0),
                                  seconds=round(time.time() - t, 1))
    if a.oms:
        out["pari_overconvergent_memory"] = oms_rows()
    with open(fn, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
