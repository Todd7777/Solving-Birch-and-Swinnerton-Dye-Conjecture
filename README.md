# The Birch and Swinnerton-Dyer conjecture — what is known, what is computed, what is open

**This repository does not contain a proof of the Birch and Swinnerton-Dyer conjecture.** Nobody has one.
It is a Millennium Prize problem and, as of November 2026, its status is exactly the one documented here.
What this repository contains instead is the most useful thing that can honestly be built around the problem:

1. **A precise, normalisation-safe statement** of the conjecture (`docs/01-statement.md`).
2. **A theorem ledger**: every major result about the conjecture, with its hypotheses, and an explicit list of
   what is *not* known (`docs/02-what-is-proven.md`).
3. **A structural analysis of why every known method stops at rank one** (`docs/03-obstacles-and-strategies.md`).
4. **An honest research program** with well-posed intermediate targets and a list of the standard ways attempted proofs fail
   (`docs/04-research-program.md`).
5. **A computational verification toolkit** (`bsd/`, built on PARI/GP via cypari2) that computes every quantity in the BSD
   formula for explicit curves, checks the formula numerically, certifies what can be certified (exact vanishing of $L(E,1)$
   by modular symbols, rank bounds by 2-descent, finiteness of $Ш[p^\infty]$ via Kato's theorem for rank $\ge2$ curves), and
   mechanically checks the hypotheses of the Skinner–Urban and Jetchev–Skinner–Wan theorems prime by prime.
6. **Results** of running the toolkit on Cremona's database and on the record-rank curves (`results/`, summarised in
   `docs/05-numerical-evidence.md`).

The distinction between *theorem*, *certified computation*, and *numerical evidence* is maintained everywhere; the code prints
it for every curve.

## Status of the conjecture in one table

| Situation | BSD (I): rank = analytic rank | $Ш$ finite | BSD (II): leading term |
|---|---|---|---|
| analytic rank 0 | **theorem** (Kolyvagin + Gross–Zagier; Kato) | **theorem** | $p$-part proven for most $p$ (Skinner–Urban + Kato; Rubin for CM); all $p$ only curve-by-curve |
| analytic rank 1 | **theorem** (Gross–Zagier + Kolyvagin) | **theorem** | $p$-part proven for $p\ge5$ good, semistable curves (Jetchev–Skinner–Wan); all $p$ only curve-by-curve |
| analytic rank 2 | verifiable curve by curve (e.g. 389a1), **no general theorem** | open for every curve (finite at individual primes only) | open |
| analytic rank $\ge3$ | not even verifiable for a single curve (vanishing of $L'$, $L''$ cannot be certified) | open | open |
| "most" curves | $\ge 66.48\%$ of all $E/\mathbb{Q}$ (Bhargava–Skinner–Zhang); 100 % of twists in suitable families (Smith) | same | — |
| function fields | **equivalent to finiteness of $Ш$** (Tate–Milne); many families proven | — | — |

## Quick start

```bash
pip install -r requirements.txt          # cypari2 (bundles PARI/GP 2.17), cysignals, matplotlib, pytest
python -m pytest -q                      # 19 regression tests against known values
python -m bsd verify 11a1 37a1 389a1 5077a1 571a1 --padic 13
python -m bsd curated --out results/curated --padic 97
python -m bsd sweep --min 1 --max 1000 --workers 4 --out results/sweep_0001_1000
python -m bsd report results/sweep_0001_1000.jsonl --md results/sweep_0001_1000.md
python -m bsd heuristic 11a1 37a1 389a1 5077a1 --X 300000 --plot results/heuristic.png
python -m bsd padic 389a1 --primes 3,5,7,11
```

Curves of conductor $< 20000$ can be given by Cremona label (the coefficients are fetched from
[JohnCremona/ecdata](https://github.com/JohnCremona/ecdata) and cached in `data/ecdata/`); any curve can be given by
a-invariants, e.g. `python -m bsd verify 0,0,1,-79,342` (the rank-5 curve of conductor 19047851).

## What `verify` computes for a curve

* minimal model, conductor, Kodaira symbols, Tamagawa numbers $c_p$, real period $\Omega_E$, torsion, root number, rational isogeny degrees;
* analytic rank and $L^{(r)}(E,1)/r!$ (PARI `ellanalyticrank`/`ellL1`), plus the **exact rational** $L(E,1)/\Omega^+$ from modular symbols;
* what is *rigorously* known about $\operatorname{ord}_{s=1}L(E,s)$ (e.g. "$\ge1$, odd" for a rank-3 curve — that is all anyone can prove);
* 2-descent rank bounds and the 2-rank of $Ш[2]/2Ш[4]$ (PARI `ellrank`), a Mordell–Weil basis saturated at all primes $<1000$, the regulator;
* the BSD quotient $\#Ш_{\mathrm{an}}$, rounded and tested for being a perfect square (Cassels);
* the status of the rank conjecture and of finiteness of $Ш$ *for this curve*, citing the theorems that apply;
* for analytic rank $\le1$: the prime-by-prime ledger of the $p$-part of the formula under the mechanised theorems;
* with `--padic B`: $p$-adic BSD data for good primes $3\le p\le B$ and, for rank $\ge2$ curves, the certificate that $Ш[p^\infty]$ is finite
  (Kato's divisibility + Mazur's control theorem + a certified upper bound for $\operatorname{ord}_{T=0}L_p$);
* a comparison with Cremona's independently computed table values.

## Layout

```
bsd/            the package (pari_session, ecdata, invariants, analytic, descent, theorems, padic, heuristic, verify, report, cli)
docs/           01 statement · 02 what is proven · 03 obstacles · 04 research program · 05 numerical evidence · references
results/        JSONL + Markdown output of the runs reported in docs/05, the heuristic plot
scripts/        batch runners
tests/          pytest regression suite
```

## Caveats, stated plainly

* Numerical non-vanishing ($L^{(k)}(E,1)\ne0$) is asserted when the value exceeds the working precision by many orders of magnitude; PARI does not do interval arithmetic.
* Saturation of Mordell–Weil bases is verified at all primes below 1000 and cross-checked against Cremona's independently saturated bases; a provable index bound (Cremona–Siksek) is not implemented.
* The theorem ledger mechanises the hypotheses of three theorems (Skinner–Urban 2014, Jetchev–Skinner–Wan 2017, Castella 2018) and Kato's rational divisibility. "Not covered" means not covered *by those*, not unknown in the literature.
* Bibliographic details in `docs/references.md` were compiled from memory and spot-checked; the mathematical statements were checked against the papers' stated hypotheses.

## License

No license file has been added; the repository owner should choose one. Cremona's database is used under the Artistic License 2.0.
