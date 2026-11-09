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
python -m pytest -q          # 2026-11-08: 21 pass, 1 skipped
```

The skipped test needs Sage (`pip install passagemath-standard`, a large download).

## State on 2026-11-08

- Commit `e984137` stopped tracking in-progress agent work until "integration": `research/`,
  `bsd/gz_certificate.py`, `bsd/rigorous_l.py` and their tests. Integration never reached GitHub. The last
  tracked versions are in `e984137^` (for example `git show e984137^:bsd/rigorous_l.py`). To integrate them,
  restore the files, remove their lines from `.gitignore`, and re-run the tests.

## Where the next steps are

`docs/04-research-program.md` and `docs/03-obstacles-and-strategies.md`.
