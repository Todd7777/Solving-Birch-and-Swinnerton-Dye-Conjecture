"""Aggregate data/sweep*_shard*.jsonl (+ data/pass2*.jsonl overrides) into results.json, results.csv
and class_table.md.  Needs Sage only to decide surjectivity at primes in Sage's non_surjective() list.

usage: python3 aggregate.py
"""
from __future__ import annotations

import csv
import glob
import json
import os
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    recs = {}
    for fn in sorted(glob.glob(os.path.join(HERE, 'data', 'sweep_shard*.jsonl'))):
        for line in open(fn):
            d = json.loads(line)
            recs[d['label']] = d
    # pass-2 overrides: one JSON line per (label, p) with a new prime record
    for fn in sorted(glob.glob(os.path.join(HERE, 'data', 'pass2*.jsonl'))):
        for line in open(fn):
            d = json.loads(line)
            lab, p = d['label'], d['p']
            if lab not in recs:
                continue
            prs = recs[lab]['primes']
            for i, q in enumerate(prs):
                if q['p'] == p:
                    if d.get('certified_finiteness') or not q.get('certified_finiteness'):
                        d['pass'] = 2
                        prs[i] = d
            else:
                if all(q['p'] != p for q in prs):
                    d['pass'] = 2
                    prs.append(d)
            prs.sort(key=lambda q: q['p'])
    # pass-3 (larger p) additions
    for fn in sorted(glob.glob(os.path.join(HERE, 'data', 'pass3*.jsonl'))):
        for line in open(fn):
            d = json.loads(line)
            lab = d['label']
            if lab not in recs:
                continue
            prs = recs[lab]['primes']
            if all(q['p'] != d['p'] for q in prs):
                d['pass'] = 3
                prs.append(d)
            prs.sort(key=lambda q: q['p'])
    return recs


def surjectivity(recs):
    """returns dict label -> set of odd primes p<100 (good) with rho-bar NOT surjective (or undecided)."""
    from sage.all import EllipticCurve
    out = {}
    for lab, d in recs.items():
        ns = d.get('non_surjective_list')
        if d.get('cm_disc') is not None:
            out[lab] = 'CM'
            continue
        bad = set()
        if isinstance(ns, list) and ns:
            E = EllipticCurve(d['ainvs'])
            rho = E.galois_representation()
            for p in ns:
                if p == 2:
                    continue
                s = rho.is_surjective(p)
                if s is not True:
                    bad.add(int(p))
        out[lab] = bad
    return out


def classify(d, q, nonsurj):
    """category of a (class, prime) record."""
    p = q['p']
    t = q.get('type', '')
    if t.startswith('bad'):
        return 'bad reduction'
    cm = d.get('cm_disc') is not None
    if not q.get('certified_finiteness'):
        if 'a_p != 0' in q.get('reason', ''):
            return 'not certified: supersingular a_3=+-3'
        if 'precision' in q.get('reason', ''):
            return 'not certified: precision'
        if 'INCONSISTENT' in q.get('reason', ''):
            return 'INCONSISTENT'
        return 'not certified: other'
    if cm:
        return 'certified (CM: needs Rubin / Pollack-Rubin)'
    if t == 'good ordinary':
        if p == 3:
            return 'certified ordinary p=3 (outside SW13 range)'
        if p in nonsurj:
            return 'certified ordinary, rho-bar not surjective (outside SW13)'
        return 'certified ordinary, surjective p>=5 (reproduces SW13, who prove Sha[p]=0)'
    if t == 'good supersingular':
        return 'certified supersingular a_p=0 (Kobayashi)'
    return 'certified other'


def main():
    recs = load()
    ns = surjectivity(recs)
    rows = []
    cat_count = Counter()
    per_p = defaultdict(Counter)
    classes = []
    for lab in sorted(recs, key=lambda l: (recs[l]['N'], recs[l]['iso'])):
        d = recs[lab]
        if 'error' in d:
            classes.append(dict(label=lab, N=d['N'], error=d['error']))
            continue
        nsj = ns.get(lab, set())
        cert, notc = [], []
        for q in d['primes']:
            cat = classify(d, q, nsj if nsj != 'CM' else set())
            cat_count[cat] += 1
            per_p[q['p']][cat] += 1
            tries = q.get('tries', [])
            rows.append(dict(label=lab, N=d['N'], rank=d['rank_cremona'], cm=d.get('cm_disc'), p=q['p'], type=q.get('type'),
                             a_p=q.get('a_p'), certified=bool(q.get('certified_finiteness')), category=cat,
                             method=('ordinary Riemann sum level p^%s' % q.get('n') if q.get('type') == 'good ordinary' and q.get('certified')
                                     else ('%s from theta_%s' % (q.get('function'), q.get('n')) if q.get('type') == 'good supersingular' and q.get('certified') else '')),
                             val_ar=q.get('val_ar'), prec_ar=q.get('prec_ar'),
                             levels_tried=';'.join(str(t.get('n')) for t in tries),
                             evals=q.get('evals'), secs=q.get('secs'), reason=q.get('reason', ''),
                             rho_surjective=(None if nsj == 'CM' else (q['p'] not in nsj)),
                             pass_=q.get('pass', 1)))
            if q.get('certified_finiteness'):
                cert.append(q['p'])
            else:
                notc.append((q['p'], q.get('reason', '')))
        classes.append(dict(label=lab, N=d['N'], rank=d['rank_cremona'], cm_disc=d.get('cm_disc'),
                            independence=d['independence'], content_g=d['modsym']['content_g'],
                            sanity=d['modsym']['sanity'], nonsurjective_odd_primes=(sorted(nsj) if nsj != 'CM' else 'CM'),
                            certified_primes=cert, not_certified=notc, secs=d['secs']))
    summary = dict(n_classes=len(classes), categories=dict(cat_count),
                   per_prime={p: dict(c) for p, c in sorted(per_p.items())},
                   n_class_prime_certified=sum(1 for r in rows if r['certified']),
                   n_errors=sum(1 for c in classes if 'error' in c))
    json.dump(dict(summary=summary, classes=classes), open(os.path.join(HERE, 'results.json'), 'w'), indent=1)
    with open(os.path.join(HERE, 'results.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    # markdown table per class
    with open(os.path.join(HERE, 'class_table.md'), 'w') as fh:
        fh.write('# Per-class table (generated by aggregate.py)\n\n')
        fh.write('Odd primes p < 100 (and the pass-3 primes, if any). "cert" = Sha(E/Q)[p^oo] proven finite for every curve in the class '
                 '(see report.md for the exact chain of theorems; supersingular and CM entries are conditional as explained there).\n'
                 'Reasons: bad = bad reduction (not treated); ss3 = supersingular with a_3 = +-3 (not treated); '
                 'prec = T^r coefficient vanished to the guaranteed precision at all levels tried.\n\n')
        fh.write('| class | N | r | CM | primes certified | not certified (reason) | secs |\n|---|---|---|---|---|---|---|\n')
        for c in classes:
            if 'error' in c:
                fh.write(f"| {c['label']} | {c['N']} | | | ERROR {c['error']} | | |\n")
                continue
            def short(reason):
                if reason.startswith('bad'):
                    return 'bad'
                if 'a_p != 0' in reason:
                    return 'ss3'
                if reason.startswith('precision'):
                    return 'prec'
                return reason[:20]
            nc = ', '.join(f'{p} ({short(r)})' for p, r in c['not_certified'])
            fh.write(f"| {c['label'][:-1]} | {c['N']} | {c['rank']} | {c['cm_disc'] or ''} | {','.join(map(str, c['certified_primes']))} | {nc} | {c['secs']} |\n")
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
