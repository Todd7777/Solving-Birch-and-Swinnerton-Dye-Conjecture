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
python -m pytest -q          # 2026-11-09: 91 pass, 1 skipped (about 20 s)
```

The skipped test needs Sage (`pip install passagemath-standard`, a large download). `research/kato-sha/` and the
height bound of `bsd/gz_certificate.py` also need Sage; nothing else does.

## State on 2026-11-09

- The stranded agent work (`bsd/rigorous_l.py`, `bsd/gz_certificate.py`, `research/`) is integrated and tracked.
- `bsd/bsd1.py` certifies rank = analytic rank. It has been run on every isogeny class with N < 20000; all
  81,217 classes are certified (`results/bsd1/`, `docs/06-certified-results.md`). The raw JSONL is not tracked;
  rerun takes 16 CPU-minutes.
- `bsd/kato_bound.py` proves Selmer-corank bounds from exact p-adic L-function Riemann sums, for ordinary,
  multiplicative and supersingular (a_p = 0) p. Its theorem statements follow Stein-Wuthrich 2013, sections 3 and
  6-7, which were read; Kato's own paper and Pollack's Duke paper were not. The supersingular congruence
  L^- == -theta_1 mod omega_1 is proved in `SUPERSINGULAR_PROOF` there and checked exactly in the tests.
- `research/t2-389a1/` holds the finiteness of Sha(389a1)[p^oo] for every prime p < 5000; see `report.md` there.
- PARI's `ellrank` uses GRH via `bnfinit` of the cubic 2-division field (read in PARI's `ellrank.c`). Certify
  that field with `bnfcertify` before calling an `ellrank` bound unconditional, as `bsd1.descent_bound` does.

## Where the next steps are

`docs/04-research-program.md` (section 4.5 is the progress log) and `docs/03-obstacles-and-strategies.md`.
Concrete open items left by the 2026-11-09 session:
- extending T2 for 389a1 beyond p = 5000. Each prime costs about p^2 symbol values; supersingular primes cost the
  same, and a coefficient divisible by p needs the next level (p^3 or p^4);
- Sha[p] = 0, rather than finiteness, needs the integral form (SW Thm 7.3) and the p-adic regulator's valuation;
- running `bsd1` on N < 500000 with the FASRC kit in `scripts/fasrc/`.
