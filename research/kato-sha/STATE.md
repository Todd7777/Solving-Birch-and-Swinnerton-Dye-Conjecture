# STATE (kato-sha)

- 2026-11-07 03:33 : fresh start (directory was empty). Exploring code.
- 03:58 : literature: Stein-Wuthrich 2013 main thm (per WebSearch summary) = Sha(E/Q)[p]=0 for 1,534,422 pairs (E,p),
  E non-CM rank>=2 N<=30000, p good ordinary 5<=p<1000 with rhobar_{E,p} surjective. So our p>=5 ordinary surjective
  results are REPRODUCTIONS (weaker). New territory: p=3, non-surjective p, CM, supersingular.
  Skinner AWS 2018 notes downloaded (GitHub raw) to scratchpad dl_skinner/: confirm Kato Div-1 rational version holds
  even when E[p] reducible; Pollack L^+- in Lambda with interpolation formulas (sec 4.1.4); Kobayashi control (Prop 14).
- Plan: Sage eclib modular symbols (2.6s at N~1e4, ~15us/eval) + own Riemann sums with own precision analysis,
  denominators nu_p = min v_p over Manin-symbol values. Supersingular a_p=0: theta_1 => L^- mod p (Pollack/Kobayashi).
- 04:08 : kato_lib.py + run_sweep.py written and validated (ordinary coeffs == Sage series to stated precision on 11a1, 389a1;
  supersingular theta_1/theta_2 agree with Sage L_alpha digits on 1483a1, 1613a1, 1531a1). Pass-1 sweep launched:
  `python3 run_sweep.py --shard {0,1} --nshards 2 --pmax 100` -> data/sweep_shard{0,1}.jsonl (resumable: rerun same cmd).
  logs in logs/sweep_shard*.log

## Integration (2026-11-09)
- Requires Sage (eclib modular symbols). Pass-1 data in data/; scratch and logs dropped.
- Next: re-implement the ordinary-prime part without Sage (PARI msfromell/mseval are exact and fast,
  about 7 microseconds per symbol at N = 389); supersingular primes (Pollack L^+-) would remain here.
