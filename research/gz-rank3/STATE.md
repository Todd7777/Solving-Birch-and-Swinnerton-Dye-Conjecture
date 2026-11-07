# STATE (gz-rank3)
Started 2026-11-07. Scratch scripts in /tmp/claude-0/-home-user/49d47a5d-6794-5fd0-9790-0bd7cf4fd119/scratchpad/gz (may be lost; key results copied here).
## Done
- Rank-3 curves N<20000 (Cremona): 5077a1 11197a1 11642a1 12279a1 13766a1 16811a1 18097b1 18562c1 18745a1 (each class has 1 curve).
- PARI twisted modular symbol convention verified (7 curves, ranks 0-3):
  sum_{0<=a<|D|} (D|a) x^-(path oo -> a/|D|) = L(E^D,1) sqrt|D| / Omega^-   [PARI 2.17.2 help text omits sqrt|D|]
  (Omega^- = (-2/c_oo) imag(E.omega[2])); opposite path orientation flips sign.
- GZ normalisation check (PARI only): 120 pairs (rank-1 optimal curves N<=300, 3 odd Heegner D each):
  hK := sqrt|D| L'(E,1) L(E^D,1) / (2 E.area); hK/(2 Reg_Q) is always a perfect square (1,4,9,16,36,49,64,81) or 0 (11 cases, L(E^D,1)=0). => hhat_K = 2 hhat_Q (BSD/Cremona normalisation, PARI ellheight).
- Sage: heegner_point_height uses alpha = sqrt|D|/(2 complex_area) (same constant). Sage height_function().min(tol,n_max) = rigorous lower bound for all non-torsion points (CS2006/Tho2010; divides by tp^2).
- python-flint 0.9.0 available (arb.expint, arb.agm, fmpz_poly.complex_roots) -> use for interval arithmetic.
- bsd/gz_certificate.py written (certify_rank3); tests/test_gz_certificate.py written: 12 passed in 3 s.
- Arb expint loses up to ~60-110 bits -> +110 guard bits in lprime_ball.
- run_certificates.py (resumable, writes certificates.jsonl: primary + secondary D per curve) launched;
  5077a1, 11197a1 certified (margin 10^58.5, two D each).
## Next
- finish 9 curves; big curves via run_big.py (numpy a_n sieve, checkpointed): rank5 D=-7, rank7 D=-23 (4e7 terms).
- results.json + report.md.
