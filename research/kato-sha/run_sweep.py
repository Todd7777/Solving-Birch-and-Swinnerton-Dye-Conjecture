"""Resumable sweep: for every rank-2 isogeny class with N <= 9999 (curve #1 of the class) and 5077a1,
certify ord_{T=0} L_p <= r at odd primes p < PMAX, recording everything as JSON lines.

usage: python3 run_sweep.py --shard i --nshards k [--pmax 100] [--out data/sweep] [--labels a,b]

Each finished class is appended as one JSON line to  <out>_shard<i>.jsonl ; classes already present
are skipped, so the job can be killed and restarted at any time.
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cysignals.alarm import alarm, cancel_alarm  # noqa: E402
from cysignals.signals import AlarmInterrupt  # noqa: E402
from sage.all import EllipticCurve, prime_range  # noqa: E402

import kato_lib as K  # noqa: E402

SMALL = (3, 5, 7, 11, 13)


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def cost_ordinary(p, n):
    # new modular-symbol evaluations at level n (level n-1 values are cached)
    return (p - 1) * p ** (n - 1)


def cost_theta(p, n):
    return (p - 1) * p ** n


def do_ordinary(ms, p, ap, r, cost_cap, nstart=2, nmax=40):
    rec = dict(p=p, type='good ordinary', a_p=ap, method='MSD Riemann sums (exact), Kato',
               anomalous=(ap % p == 1))
    spent = 0
    tries = []
    n = nstart
    while n <= nmax:
        c = cost_ordinary(p, n)
        if r >= p ** (n - 1):        # need j = r < p^{n-1}
            n += 1
            continue
        if spent + c > cost_cap and tries:
            break
        t0 = time.time()
        M = n + 3
        jmax = r + 1
        co = K.ordinary_riemann(ms, p, ap, n, jmax, M)
        spent += c
        precs = [K.ordinary_precision(p, n, j) for j in range(jmax + 1)]
        vals = [K.vp_int(x, p, cap=M) for x in co]
        lower_ok = all(vals[j] >= min(precs[j], M) for j in range(r))
        cert = vals[r] < precs[r]
        tries.append(dict(n=n, coeffs_mod_pM=[int(x) for x in co], M=M, val=vals, prec=precs,
                          lower_coeffs_vanish=lower_ok, certified=cert, secs=round(time.time() - t0, 3)))
        if not lower_ok:
            rec.update(certified=False, reason='INCONSISTENT: a lower coefficient is nonzero mod its precision',
                       tries=tries, evals=spent)
            return rec
        if cert:
            rec.update(certified=True, n=n, val_ar=vals[r], prec_ar=precs[r], tries=tries, evals=spent,
                       reason=f'coeff of T^{r} has valuation {vals[r]} < guaranteed precision {precs[r]} (level p^{n})')
            return rec
        n += 1
    rec.update(certified=False, reason='precision: T^r coefficient vanished to the guaranteed precision at every level tried',
               tries=tries, evals=spent)
    return rec


def ss_window(p, n):
    """(precision exponent, max j) for the congruence of theta_n-derived G_n with L^-(n odd)/L^+(n even)."""
    if n % 2 == 1:
        return (n + 1) // 2, p - 1
    return n // 2, p * p - p


def do_supersingular(ms, p, r, cost_cap, nmax=12):
    rec = dict(p=p, type='good supersingular', a_p=0,
               method='Mazur-Tate elements -> Pollack L^+-, Kobayashi')
    spent = 0
    tries = []
    for n in range(1, nmax + 1):
        c = cost_theta(p, n)
        if spent + c > cost_cap and tries:
            break
        prec, jlim = ss_window(p, n)
        if r > jlim:
            continue
        t0 = time.time()
        cX = K.theta_X_coeffs(ms, p, n)
        spent += c
        # divide by omega~_n^{-+}: product of Phi_{p^m}, m < n, m of parity opposite to n
        rem_zero = True
        Q = cX
        for m in range(n - 1, 0, -2):
            Q, rem = K.divide_by_cyclotomic(Q, p, m)
            if any(rem):
                rem_zero = False
        sign = (-1) ** ((n + 1) // 2) if n % 2 == 1 else (-1) ** ((n + 2) // 2)
        G = [sign * x for x in K.X_to_T(Q, r + 1)]
        mod = p ** prec
        vals = [K.vp_int(x, p, cap=prec) if x % mod else prec for x in G]
        lower_ok = all(v >= prec for v in vals[:r])
        cert = vals[r] < prec
        which = 'L^-' if n % 2 == 1 else 'L^+'
        tries.append(dict(n=n, function=which, coeffs_mod=[int(x % mod) for x in G], modulus=f'{p}^{prec}',
                          val=vals, prec=prec, division_exact=rem_zero, lower_coeffs_vanish=lower_ok,
                          certified=cert, secs=round(time.time() - t0, 3)))
        if not rem_zero:
            rec.update(certified=False, reason='INCONSISTENT: cyclotomic division not exact', tries=tries, evals=spent)
            return rec
        if not lower_ok:
            rec.update(certified=False, reason='INCONSISTENT: a lower coefficient is nonzero mod p^prec',
                       tries=tries, evals=spent)
            return rec
        if cert:
            rec.update(certified=True, n=n, function=which, val_ar=vals[r], prec_ar=prec, tries=tries, evals=spent,
                       reason=f'T^{r} coefficient of {which} is nonzero mod {p}^{prec} (from theta_{n})')
            return rec
    rec.update(certified=False, reason='precision: T^r coefficients of L^+ and L^- vanished to the precision tried',
               tries=tries, evals=spent)
    return rec


def process_class(c, pmax, cap_small, cap_large, class_timeout):
    t_start = time.time()
    E = EllipticCurve(c['ainvs'])
    N = int(E.conductor())
    assert N == c['N']
    r = c['rank']
    out = dict(label=c['label'], iso=c['iso'], N=N, ainvs=c['ainvs'], rank_cremona=r, torsion=c['torsion'])
    cm = E.has_cm()
    out['cm_disc'] = int(E.cm_discriminant()) if cm else None
    # generators and exact independence
    pts = [K.point_from_proj(E, g) for g in c['gens']]
    out['gens'] = [[str(P.xy()[0]), str(P.xy()[1])] for P in pts]
    out['independence'] = K.independence_certificate(E, pts, c['torsion'])
    r_cert = len(pts) if out['independence']['ok'] else 0
    out['r_points'] = r_cert
    # Galois image info (only used to compare with Stein-Wuthrich's coverage)
    rho = E.galois_representation()
    try:
        out['non_surjective_list'] = [int(x) for x in rho.non_surjective()]
    except Exception as exc:  # pragma: no cover
        out['non_surjective_list'] = f'error {exc!r}'
    out['reducible_primes'] = [int(p) for p in prime_range(3, pmax) if rho.is_reducible(p)]
    # modular symbol
    ms = K.PrimitiveModSym(E)
    out['modsym'] = dict(content_g=str(ms.g), n_manin=ms.n_manin, t_build=round(ms.t_build, 3),
                         t_content=round(ms.t_content, 3), sanity=ms.sanity_checks())
    primes = []
    for p in prime_range(3, pmax):
        p = int(p)
        t0 = time.time()
        if N % p == 0:
            ap = int(E.ap(p))
            kind = ('split multiplicative' if ap == 1 else 'nonsplit multiplicative' if ap == -1 else 'additive')
            primes.append(dict(p=p, type='bad: ' + kind, a_p=ap, certified=False, reason='bad reduction (not treated)'))
            continue
        ap = int(E.ap(p))
        cap = cap_small if p in SMALL else cap_large
        try:
            alarm(class_timeout)
            if ap % p != 0:
                rec = do_ordinary(ms, p, ap, r_cert or r, cap)
            elif ap == 0:
                rec = do_supersingular(ms, p, r_cert or r, cap)
            else:
                rec = dict(p=p, type='good supersingular', a_p=ap, certified=False,
                           reason='supersingular with a_p != 0 (p=3, a_3=+-3): needs Sprung sharp/flat theory, not treated')
            cancel_alarm()
        except AlarmInterrupt:
            rec = dict(p=p, type='timeout', a_p=ap, certified=False, reason=f'timeout {class_timeout}s')
        rec['secs'] = round(time.time() - t0, 3)
        if not out['independence']['ok']:
            rec['certified_finiteness'] = False
        else:
            rec['certified_finiteness'] = bool(rec.get('certified'))
        primes.append(rec)
    out['primes'] = primes
    out['secs'] = round(time.time() - t_start, 2)
    out['rss_mb'] = round(rss_mb(), 1)
    del ms
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shard', type=int, default=0)
    ap.add_argument('--nshards', type=int, default=1)
    ap.add_argument('--pmax', type=int, default=100)
    ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'sweep'))
    ap.add_argument('--labels', default='')
    ap.add_argument('--cap-small', type=int, default=400000)
    ap.add_argument('--cap-large', type=int, default=150000)
    ap.add_argument('--timeout', type=int, default=600)
    ap.add_argument('--max-classes', type=int, default=10**9)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fn = f'{args.out}_shard{args.shard}.jsonl'
    done = set()
    if os.path.exists(fn):
        for line in open(fn):
            try:
                done.add(json.loads(line)['label'])
            except Exception:
                pass
    classes = K.load_classes()
    if args.labels:
        want = set(args.labels.split(','))
        classes = [c for c in classes if c['label'] in want]
    mine = [c for i, c in enumerate(classes) if i % args.nshards == args.shard and c['label'] not in done]
    print(f'[{time.strftime("%H:%M:%S")}] shard {args.shard}/{args.nshards}: {len(mine)} classes to do, '
          f'{len(done)} already done', flush=True)
    count = 0
    with open(fn, 'a') as fh:
        for c in mine:
            if count >= args.max_classes:
                break
            try:
                rec = process_class(c, args.pmax, args.cap_small, args.cap_large, args.timeout)
            except Exception as exc:
                rec = dict(label=c['label'], iso=c['iso'], N=c['N'], error=repr(exc), tb=traceback.format_exc())
            fh.write(json.dumps(rec) + '\n')
            fh.flush()
            count += 1
            ncert = sum(1 for q in rec.get('primes', []) if q.get('certified_finiteness'))
            ngood = sum(1 for q in rec.get('primes', []) if not q.get('type', '').startswith('bad'))
            print(f'[{time.strftime("%H:%M:%S")}] {rec["label"]} N={rec["N"]} cert {ncert}/{ngood} good odd p<{args.pmax}'
                  f' {rec.get("secs")}s rss {rec.get("rss_mb")}MB {"ERROR " + rec["error"] if "error" in rec else ""}',
                  flush=True)


if __name__ == '__main__':
    main()
