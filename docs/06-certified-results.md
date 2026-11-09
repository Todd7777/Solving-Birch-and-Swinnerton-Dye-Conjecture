# 6. Certified results

This file records statements that the repository's code **proves**, as opposed to the numerical evidence of
file 05. Each one is a *certified computation*: a theorem conditional only on the published theorems it cites
and on the correctness of the software it runs (PARI/GP 2.17.2, FLINT/Arb via python-flint 0.9.0). None of
them is a proof of the Birch and Swinnerton-Dyer conjecture or of any part of it beyond the finitely many curves
named. Section 6.4 says how they relate to what was already known.

## 6.1 Rank = analytic rank for every curve of conductor < 20000

**Certified computation (2026-11-09).** *For every elliptic curve $E/\mathbb{Q}$ of conductor $N<20000$,*
$$\operatorname{rank}E(\mathbb{Q}) \;=\; \operatorname{ord}_{s=1}L(E,s)\;\in\;\{0,1,2,3\}.$$

The run covers all 81,217 isogeny classes in Cremona's tables with $N<20000$ (132,535 curves). Rank and
$L(E,s)$ are isogeny invariants, so one curve per class (Cremona's curve 1) suffices.

| rank $r$ | classes | certified | how $\operatorname{ord}\ge r$ is proven | how $\operatorname{rank}\le r$ is proven |
|---|---|---|---|---|
| 0 | 33,551 | 33,551 | (nothing to prove) | Kolyvagin / Kato from $L(E,1)\ne0$ |
| 1 | 42,198 | 42,198 | root number $-1$ | Gross–Zagier–Kolyvagin from $\operatorname{ord}=1$ |
| 2 | 5,459 | 5,459 | a point of infinite order forces $L(E,1)=0$; root number $+1$ | 2-descent (all 5,459); Kato at $p\ge5$ as well (243 with $N\le3000$) |
| 3 | 9 | 9 | two independent points force $\operatorname{ord}\ne1$; root number $-1$ | 2-descent |

In every class $\operatorname{ord}\le r$ comes from an Arb ball for $L^{(r)}(E,1)/r!$ that excludes 0. In every
class of rank $\ge2$, $\operatorname{rank}\ge r$ comes from Cremona's generators, which are checked to lie on the
curve and certified independent. There were no failures and no disagreements with Cremona's ranks. Total time was
16 CPU-minutes (5.5 minutes on 3 cores).

**Consequences proven along the way.**
* *Rank $\le1$ (75,749 classes):* $Ш(E/\mathbb{Q})$ is finite (Kolyvagin; Kato). This was already a theorem once
  $\operatorname{ord}\le1$ is known. What is new here is only that the order is certified rather than observed.
* *Rank 2 and 3 (5,468 classes):* for curve 1 of each class, $\dim_{\mathbb{F}_2}\operatorname{Sel}_2 = r + \dim E(\mathbb{Q})[2]$.
  So $Ш[2]=0$, hence $Ш[2^\infty]=0$. For the other curves of the class, $Ш[2^\infty]$ is finite, because finiteness
  of $Ш[p^\infty]$ is an isogeny invariant (Selmer coranks are).
* *243 rank-2 classes with $N\le3000$:* at the good ordinary prime $p\ge5$ used for the Kato bound, $Ш[p^\infty]$
  is finite and the cyclotomic $p$-adic height pairing on $E(\mathbb{Q})$ is non-degenerate (Schneider's conjecture
  at that $p$; see 6.2(e)).

The finiteness of the whole of $Ш$ remains open for every one of the 5,468 classes of rank $\ge2$. No curve of
rank $\ge4$ has conductor $<20000$. For rank $\ge4$ no method can certify $\operatorname{ord}\ge4$ (section 6.5).

## 6.2 The proof chain

(a) **Parity** (functional equation, from modularity: Wiles, Taylor–Wiles, Breuil–Conrad–Diamond–Taylor):
$\Lambda(1+t)=w\Lambda(1-t)$, so $\operatorname{ord}_{s=1}L(E,s)\equiv(1-w)/2 \pmod 2$. $w$ comes from PARI `ellrootno`
(product of local root numbers, exact).

(b) **Upper bound for the order.** `bsd/rigorous_l.py` computes Arb balls that provably contain
$c_k=L^{(k)}(E,1)/k!$. It uses exact $a_n$ (PARI `ellan`), the Mellin integral split by the Fricke involution, a
cancellation-free form of the incomplete-gamma kernels, and proven tail bounds (B1)–(B4)
(`research/rigorous-l/report.md`, section 2). A ball for $c_r$ that excludes 0 proves $\operatorname{ord}\le r$. All
classes were certified at 64 bits, with radii near $10^{-20}$. The smallest leading coefficient among the rank-0/1
classes is $0.0478$ (6045a1).

(c) **Lower bound for the order from points.** These are two contrapositives of theorems:
* $L(E,1)\ne0 \Rightarrow E(\mathbb{Q})$ finite (Kolyvagin 1988 with Gross–Zagier 1986 and Bump–Friedberg–Hoffstein /
  Murty–Murty; independently Kato 2004). So one point of infinite order gives $L(E,1)=0$. With $w=+1$ this gives
  $\operatorname{ord}\ge2$.
* $\operatorname{ord}=1 \Rightarrow \operatorname{rank}=1$ (Gross–Zagier, Kolyvagin). So two independent points give
  $\operatorname{ord}\ne1$. With $w=-1$ this gives $\operatorname{ord}\ge3$.

  The second is the argument of Buhler–Gross–Zagier (Math. Comp. 44, 1985) for 5077a1. Kolyvagin's theorem stands in
  for their Heegner-point height bound. `bsd/gz_certificate.py` implements their original route, and
  `research/gz-rank3/certificates.jsonl` uses it to re-prove $L'(E,1)=0$ for 8 of the 9 rank-3 classes as a second proof.
  The new chain also covers 18097b1, which that run had not reached.

(d) **Independence of points**, exactly (`bsd/kato_bound.independence_certificate`). Suppose a primitive relation
$\sum n_iP_i$ were torsion. Then for a prime $\ell\nmid\#E(\mathbb{Q})_{\mathrm{tors}}$ the images in
$\prod_q E(\mathbb{F}_q)/\ell E(\mathbb{F}_q)$ would satisfy a non-trivial $\mathbb{F}_\ell$-relation. Here $q$ runs over
good primes with cyclic $\ell$-Sylow subgroup. A matrix of full rank $r$ over $\mathbb{F}_\ell$ rules this out. Used:
$\ell=2$ for 4,591 classes, $\ell=3$ for 871, $\ell=5$ for 6, with at most 9 primes $q$.

(e) **Upper bound for the rank.**
* *2-descent.* This is PARI `ellrank` applied to `ellrankinit(E)`. Reading the PARI 2.17.2 source (`src/basemath/ellrank.c`)
  shows where it depends on GRH. When the 2-division polynomial is irreducible, `makevbnf` calls `Buchall`
  (`bnfinit`) on the cubic field, and `bnfselmer` builds the 2-Selmer group from that field's class group and
  $S$-units. PARI guarantees these only under GRH. With a rational 2-torsion point, the Selmer group is built from
  Hilbert symbols over $\mathbb{Q}$ and a quadratic field, and no class group is used. Every cubic field that occurred
  (4,613 fields, $|d_K|\le 79{,}896$, class numbers up to 20) was certified with `bnfcertify` before `ellrank` was
  run on that same structure. The bound is therefore unconditional. In every class the plain Selmer bound
  $\dim\operatorname{Sel}_2-\dim E(\mathbb{Q})[2]$ already equals $r$, so the Cassels–Tate pairing step of `ellrank`
  was never needed.
* *Kato bound*, independent of descent (`bsd/kato_bound.py`). It uses exact Riemann sums of PARI's modular symbols
  for the $p$-adic $L$-function, with a proven congruence for each coefficient. For a good ordinary prime $p\ge5$ of a
  curve without CM, Kato's divisibility (in the form of Stein–Wuthrich 2013, Thm 7.3/7.4 with Prop. 7.2) and
  Schneider/Perrin-Riou (their Thm 6.1) give
  $r\le\operatorname{corank}\operatorname{Sel}_{p^\infty}\le\operatorname{ord}_T L_p$. Equality
  $\operatorname{ord}_T f_E=r$ holds iff $Ш[p^\infty]$ is finite and the $p$-adic height pairing is non-degenerate.
  The Kato bound was run for every rank-2 class with $N\le3000$ (243 classes) and certified
  $\operatorname{ord}_TL_p\le2$ in all of them: at $p=5$ (198 classes), $p=7$ (41) or $p=11$ (4); see
  `results/bsd1/classes.csv.gz`, column `kato_p`. Larger conductors were not run, because PARI's `msfromell` needs about 4 GB at $N\approx10^4$.

(f) **Rank.** For $r\le1$, $\operatorname{ord}=r$ gives $\operatorname{rank}=r$ by Gross–Zagier–Kolyvagin (Kato for
$r=0$). For $r\ge2$, (d) and (e) give $\operatorname{rank}=r$, and (b) and (c) give $\operatorname{ord}=r$.

## 6.3 Inputs and assumptions

* **Data.** Curves and generators come from Cremona's tables (`data/ecdata/`, fetched by `bsd/ecdata.py`). They are
  used only as *candidates*. Every generator is re-checked on the curve and certified independent, and Cremona's rank
  is the value being certified, not an input to any proof step.
* **Theorems.** Modularity; Gross–Zagier 1986; Kolyvagin 1988/1990 with Bump–Friedberg–Hoffstein 1990 /
  Murty–Murty 1991; Kato 2004 (Thm 17.4), as stated in Stein–Wuthrich 2013 §§6–7; Mazur's control theorem. The
  Stein–Wuthrich statements were read in the paper for this file; the others are as in file 02.
* **Software.** PARI/GP 2.17.2: `ellrootno`, `ellan`, `ellrank`, `bnfcertify`, `msfromell`/`mseval`, and point
  arithmetic over $\mathbb{F}_q$. FLINT/Arb: ball arithmetic.
* **Not assumed.** GRH (see 6.2(e)); BSD; finiteness of $Ш$; Cremona's ranks.

## 6.4 Relation to what was known

Cremona's tables already record these ranks, obtained with `mwrank` and analytic-rank computations. Analytic rank 3
was first proven for a curve (5077a1) by Buhler–Gross–Zagier in 1985. This run is therefore an **independent
re-derivation**, not a new theorem about any curve whose rank was unknown. Its use for the repository is that BSD (I)
for these 132,535 curves is now a certified statement produced and checkable here, with interval arithmetic in place
of floating-point non-vanishing. The pipeline also scales: the whole range took 16 CPU-minutes, so the FASRC kit in
`scripts/fasrc/` could run it to $N<500000$.

## 6.5 Where curve-by-curve certification stops

* **Rank $\ge4$.** For a curve of rank 4 ($w=+1$), $\operatorname{ord}\ge4$ needs $L''(E,1)=0$ exactly. Neither
  modular symbols (which see only $L(E,1)$) nor Gross–Zagier–Kolyvagin (which see only $\operatorname{ord}\le1$) can
  give that. For rank 5 one needs $L'''(E,1)=0$. So **no elliptic curve is known to have analytic rank $\ge4$**. For
  234446a1 (rank 4) and 19047851a1 (rank 5) the repository proves only $\operatorname{ord}\in\{2,4\}$ and
  $\operatorname{ord}\in\{3,5\}$, with rigorous bounds $|L''(E,1)|\le\varepsilon$ and $|L'''(E,1)|\le\varepsilon$
  (`research/rigorous-l/beyond.jsonl`).
* **Finiteness of $Ш$ for rank $\ge2$** is reachable only one prime at a time (6.1, 6.6). That is the subject of
  target T2 in file 04.

## 6.6 Target T2 for 389a1: $Ш[p^\infty]$ prime by prime

*Pending: filled in from `research/t2-389a1/` when the sweep over $p<3000$ completes.*

## 6.7 Reproduce

```bash
python -m bsd.bsd1 --min 1 --max 19999 --out results/bsd1/bsd1_raw.jsonl --workers 3 --kato-max-conductor 3000
python -c "from bsd import bsd1; bsd1.export('results/bsd1/bsd1_raw.jsonl', 'results/bsd1')"
python -m pytest -q tests/test_bsd1.py tests/test_kato_bound.py tests/test_rigorous_l.py
```

Tracked outputs: `results/bsd1/classes.csv.gz` (one row per class: leading-coefficient ball, independence primes,
descent bound, Kato prime), `results/bsd1/rank2plus.jsonl.gz` (full certificates for rank $\ge2$), and
`results/bsd1/summary.json`.
