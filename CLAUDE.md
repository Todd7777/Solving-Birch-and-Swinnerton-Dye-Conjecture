# CLAUDE.md

Guidance for Claude Code sessions in this repository, local or cloud.

## What this is

A theorem ledger, obstruction analysis, research program, and computational verification toolkit
(PARI/GP through cypari2) for the Birch and Swinnerton-Dyer conjecture over Q. The conjecture is open, and
this repository does not prove it.

## Ground rules

- Do not claim or imply that the conjecture is proved. Keep **theorem**, **certified computation** and
  **numerical evidence** apart in prose and in code output, as the toolkit already does for every curve.
  Cite a theorem only with its hypotheses, and only after reading the source.
- A pattern across curves is evidence, not proof. Try to construct a counterexample before building on it.
- Run the checks before and after every change. A new failure is a regression until shown otherwise.
- Commit verified work in small steps and push before the session ends. Work left only inside a cloud
  session can be lost with it.

## Setup and checks

```bash
pip install -r requirements.txt && pip install -e .
python -m pytest -q          # 2026-11-13: 122 pass, 1 skipped (about 25 s)
```

The skipped test needs Sage (`pip install passagemath-standard`, a large download). `research/kato-sha/` and the
height bound of `bsd/gz_certificate.py` also need Sage; nothing else does.

## State on 2026-11-13

- The stranded agent work (`bsd/rigorous_l.py`, `bsd/gz_certificate.py`, `research/`) is integrated and tracked.
- `bsd/bsd1.py` certifies rank = analytic rank. It has been run on every isogeny class with N < 20000; all
  81,217 classes are certified (`results/bsd1/`, `docs/06-certified-results.md`). The raw JSONL is not tracked;
  rerun takes 16 CPU-minutes.
- `bsd/kato_bound.py` proves Selmer-corank bounds from exact p-adic L-function Riemann sums, for ordinary,
  multiplicative and supersingular (a_p = 0) p. Its theorem statements follow Stein-Wuthrich 2013, sections 3 and
  6-7, which were read; Kato's own paper and Pollack's Duke paper were not. The supersingular congruence
  L^- == -theta_1 mod omega_1 is proved in `SUPERSINGULAR_PROOF` there and checked exactly in the tests.
- `research/t2-389a1/` holds the finiteness of Sha(389a1)[p^oo] for every prime p < 5000 (`report.md`) and, stronger,
  Sha(389a1)[p] = 0 for every prime p < 5000, and with SW Thm 12.3 every p < 48859, and beyond it every
  p < 120000 by direct certificates (`beyond/run_beyond.py`, compiled level-2 sums) (`sha_zero_report.md`;
  `bsd/sha_zero.py` implements Stein-Wuthrich Algorithm 11.1 with the proofs of Lemma A, Lemma B and saturation
  in its docstring). Stein-Wuthrich Thm 12.3
  already proves Sha(389a1)[p] = 0 at p = 2 and every good ordinary p < 48859 except 16231; the primes new here are
  389, the supersingular ones, SW's exception 16231 and every prime from 48859 to 120000. 16231 needed [T^2]L_p mod p^2: the compiled level-3 sum in
  `research/t2-389a1/level3/` (C + AVX-512, exact; run_level3.py validates it and proves the Fricke symmetry it
  uses) gave ord_p [T^2]L_p = 1 in 12 CPU-hours. Read both SW Thm 1.1 and Thm 12.3 before claiming a range is new.
- PARI's `ellrank` uses GRH via `bnfinit` of the cubic 2-division field (read in PARI's `ellrank.c`). Certify
  that field with `bnfcertify` before calling an `ellrank` bound unconditional, as `bsd1.descent_bound` does.

## Where the next steps are

`docs/04-research-program.md` (section 4.5 is the progress log) and `docs/03-obstacles-and-strategies.md`.
Concrete open items left by the 2026-11-12 session:
- extending Sha(389a1)[p] = 0 beyond 120000: `beyond/run_beyond.py run --pmax X` is resumable. A prime costs about p^2/4
  symbol values (with the lean kernel of 2026-11-13, about 75 s at p = 100000 and 95 s at p = 120000; check
  a kernel change with `run_beyond.py validate-lean`). Supersingular primes above 75000 need more than
  1 GB of PARI stack for the D_p-valued regulator (`parisizemax` is set to 12 GB). A prime where c_2 = 0 mod p needs
  `research/t2-389a1/level3/` (about p^3/4 symbol values: 12 CPU-hours at 16231, so it grows quickly). PARI's
  overconvergent symbols are not an option at level 389: their memory grows like p^2;
- running `bsd1` on N < 500000 with the FASRC kit in `scripts/fasrc/`.
