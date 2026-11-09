# STATE of research/rigorous-l (restart-safe log)

Task: bsd/rigorous_l.py (rigorous Arb enclosures of c_k = L^(k)(E,1)/k!), certification runs,
results.json, report.md, tests/test_rigorous_l.py.

## Status (update after each step)
- [x] 0. Studied package (pari_session, ecdata, analytic.exact_L1_ratio, curated.jsonl)
- [x] 1. Math derived & verified: Kummer form x^{-1-t}Gamma(1+t,x) = x^{-1-t}Gamma(1+t) - e^{-x} sum_m x^m/(1+t)_{m+1};
        validated vs closed forms (J0, J1=E1/x), Arb hypgeom, acb.integral; Fricke sign g(1/u)=w u^2 g(u) checked.
- [x] 2. bsd/rigorous_l.py: taylor_enclosures, certify_order_rigorous (+ external_zeros / external_lower_bound hook),
        modsym_L1, lambda_by_integration (independent rigorous method), ball_to_json, CLI (python3 -m bsd.rigorous_l 389a1)
- [x] 3. validation.jsonl DONE (25 curves, 300/300 PARI values inside balls, rad <= 2.8e-40 @128)
- [x] 3b. timing.jsonl DONE (N=19990: 0.27 s @128 bits, 4.1 s @512)
- [x] 4a. rank3.jsonl DONE (9 rank-3 curves N<20000; |L'(E,1)| <= ~3e-157 @512; c_3 != 0)
- [x] 4b. crosscheck.jsonl DONE (series vs acb.integral: 25/25 curves overlap)
- [~] 4c. beyond.jsonl (curated N>=20000 ranks 4,5,6; upper bounds): 234446a1, 19047851a1 done; rank-6 running
- [~] 4. certify.jsonl RUNNING (started 03:52; ~1000/2196 at 04:08; resumable)
        relaunch if killed:  cd repo && python3 research/rigorous-l/run_batch.py certify   (skips done labels)
- [ ] 5. FINAL: python3 research/rigorous-l/aggregate.py && python3 research/rigorous-l/make_report.py
        (both scripts work on partial data; re-run after certify completes)
- [x] 6. tests/test_rigorous_l.py: 15 tests, ~2 s (37 pass together with test_known_values.py)

## Notes
- cypari2: per-call precision is the `precision=` keyword in BITS (set_real_precision affects only GP strings).
- PARI Lambda (lfunlambda) = 2 x our Lambda (gamma factor Gamma_R(s)Gamma_R(s+1) = 2(2pi)^-s Gamma(s)).
- modsym (msfromell) at N~10^4: 2-9 s, needs > 1 GiB stack (batch uses 2 GiB = cap).
- Other agents' sweep processes use ~4.4 GB RSS + 3 cores; run only ONE heavy process.
- Open issue for parent: bsd/curated.py labels 67a1 'rank 1' but Cremona rank = 0 and ord = 0 is certified.
- eclib modular symbols NOT used as certificate (mod-p linear algebra + rational reconstruction).

## Integration (2026-11-09)
- Integrated into the tracked tree. The certify run stopped at 1001/2196 labels when the earlier session
  ended; its records are archived as `certify.jsonl.gz` (aggregate.py reads it). It is superseded by the
  complete run of bsd/bsd1.py over every isogeny class with N < 20000 (results/bsd1/, docs/06), which uses
  this module for every analytic certificate.
- `results.json` is derived: regenerate with `python3 research/rigorous-l/aggregate.py`.
