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
$L(E,s)$ are isogeny invariants, so one curve per class (Cremona's curve 1) suffices. "Every curve" relies on the
completeness of those tables (see 6.3). Strictly, the proven statement is about the 132,535 curves listed there.

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

* **Completeness of the tables.** It is not re-derived here that Cremona's tables contain every isogeny class of
  conductor $<20000$. By modularity, isogeny classes correspond to rational weight-2 newforms, and Cremona enumerated
  those with modular symbols. This is the only input the statement "every curve of conductor $<20000$" takes on trust.
* **Data.** Curves and generators come from Cremona's tables (`data/ecdata/`, fetched by `bsd/ecdata.py`). They are
  used only as *candidates*. Every generator is re-checked on the curve and certified independent, and Cremona's rank
  is the value being certified, not an input to any proof step.
* **Theorems.** Modularity; Gross–Zagier 1986; Kolyvagin 1988/1990 with Bump–Friedberg–Hoffstein 1990 /
  Murty–Murty 1991; Kato 2004 (Thm 17.4), as stated in Stein–Wuthrich 2013 §§6–7; Mazur's control theorem. The
  Stein–Wuthrich statements were read in the paper for this file; the others are as in file 02.
* **Software.** PARI/GP 2.17.2: `ellrootno`, `ellan`, `ellrank`, `bnfcertify`, `msfromell`/`mseval`, and point
  arithmetic over $\mathbb{F}_q$. FLINT/Arb: ball arithmetic.
* **Not assumed.** GRH (see 6.2(e)); BSD; finiteness of $Ш$; Cremona's ranks, generators or $L$-values.

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
  give that. For rank 5 one needs $L'''(E,1)=0$. This is why, as is usually stated in the literature (not
  re-checked against a source here), no elliptic curve is known to have analytic rank $\ge4$. The two record curves
  show exactly where things stop (`research/rigorous-l/beyond_certified.json`, same chain as 6.2):
  * 234446a1 (rank 4): $\operatorname{rank}E(\mathbb{Q})=4$ is proven (4 independent points;
    certified 2-descent gives $\le4$), and $\operatorname{ord}\in\{2,4\}$ with $|L''(E,1)|\le5.5\cdot10^{-22}$ and
    $L^{(4)}(E,1)/4!\in[8.94384739590088904\pm2\cdot10^{-20}]$. **BSD (I) for 234446a1 is equivalent to $L''(E,1)=0$.**
  * 19047851a1 (rank 5): $\operatorname{rank}=5$ is proven, and $\operatorname{ord}\in\{3,5\}$ with
    $|L'''(E,1)|\le5.0\cdot10^{-22}$. **BSD (I) for 19047851a1 is equivalent to $L'''(E,1)=0$.**
* **Finiteness of $Ш$ for rank $\ge2$** is reachable only one prime at a time (6.1, 6.6). That is the subject of
  target T2 in file 04.

## 6.6 Target T2 for 389a1: $Ш[p^\infty]$ prime by prime

Target T2 of file 04 asks for the finiteness of $Ш(E/\mathbb{Q})$ for one curve of rank 2 at *all* primes at once.
This section records how far a prime-by-prime computation gets for $E=$ 389a1 ($y^2+y=x^3+x^2-2x$, rank 2, the
smallest conductor of rank 2). It does not reach "all primes", and no computation can.

**Certified computation (2026-11-09; extended 2026-11-10).** *$Ш(389\mathrm{a}1)[p^\infty]$ is finite for every prime $p<5000$*: all 669
of them, with $Ш[2^\infty]=0$. At every odd $p<5000$ the $p$-adic regulator of $E(\mathbb{Q})$ is moreover non-zero:
for $p\ne107,599,1049,2957$ it is the cyclotomic height pairing that is non-degenerate, and at those four supersingular
primes the $D_p$-valued regulator is non-zero.

| primes | how |
|---|---|
| $p=2$ | $\dim\operatorname{Sel}_2=2$, $E(\mathbb{Q})[2]=0$, so $Ш[2]=0$ (certified 2-descent, 6.2(e)) |
| good ordinary $5\le p<5000$ (662 primes) | exact level-2 Riemann sum: $c_2\not\equiv0\pmod p$, so $\operatorname{ord}_TL_p\le2$; Stein–Wuthrich Thm 7.3/7.4 (Kato) give $\operatorname{ord}_Tf_E\le2$; their Thm 6.1 (Schneider, Perrin-Riou) gives $\operatorname{ord}_Tf_E\ge2$, with equality iff $Ш[p^\infty]$ is finite and the $p$-adic height is non-degenerate |
| $p=3$ (good ordinary, $a_3=-2$) | the same, with Thm 7.3, which needs $\rho_3$ surjective. That is proven in `research/t2-389a1/rho3_surjective.py`; see below |
| $p=389$ (split multiplicative) | `kato_bound.multiplicative_bound`: $c_0=c_1=c_2=0$, $c_3\equiv321\pmod{389}$, so $\operatorname{ord}_TL_p\le3$. Thm 7.3 in the split case ($T$ divides $L_p/f_E$) gives $\operatorname{ord}_Tf_E\le2$. $\rho_{389}$ is surjective because $E$ has no 389-isogeny (Serre, as quoted in SW §7; SW Prop. 7.2) |
| supersingular ($a_p=0$): 107, 599, 1049, 2957 | `kato_bound.supersingular_bound`: $[T^2]L^-\equiv-[T^2]\theta_1\not\equiv0\pmod p$ (residues $21, 56, 847, 2547$), so $\operatorname{ord}_TL_p\le2$ for the $D_p$-valued $L_p$; SW Thm 7.3/7.4 and Thm 6.2 (Perrin-Riou) as above. At 107 and 599, $L^+$ from $\theta_2$ gives a second, independent certificate. See below |

**$\rho_3$ is surjective** for 389a1. The proof has three parts.
1. $\bar\rho_3$ is surjective: $\psi_3$ is irreducible with Galois group $S_4$, and $SL_2(\mathbb{F}_3)$ has no subgroup of index 2.
2. Surjectivity mod 9. Elkies (arXiv:math/0612734) shows that the curves with $\rho_3$ surjective mod 3 but not mod 9
   are those with $j=f(x)$ for an explicit $f\in\mathbb{Q}(x)$ of degree 27. Brute force over $GL_2(\mathbb{Z}/9)$
   confirms the group theory behind this: all 27 lifts of $SL_2(\mathbb{Z}/3)$ are conjugate, $\mathfrak{sl}_2(\mathbb{F}_3)$
   is irreducible, and the normaliser of the lift is Elkies' $G'$. $j(389\mathrm{a}1)=1404928/389$ is not a value of $f$
   on $\mathbb{P}^1(\mathbb{Q})$; this is checked by exact factorisation. The transcription of $f$ was checked against
   both of Elkies' forms and his table of values.
3. The standard mod-9 lifting lemma, whose proof is in the script's docstring, then gives surjectivity onto $GL_2(\mathbb{Z}_3)$.

The one input taken on trust here is the moduli interpretation of $X_9=X(9)/G'$.

**The supersingular primes.** At a supersingular $p$ with $a_p=0$ the measure $\mu_\alpha$ is unbounded,
$v_p(\alpha)=\tfrac12$. Riemann sums for $L_\alpha$ itself would need about $p^5$ symbol values before they could
resolve the $T^2$ coefficient. Stein–Wuthrich (Remark 12.5) say their algorithms apply here but their implementation
was too slow, and they did not run it at scale. The way around is Pollack's splitting of $L_\alpha$ into two
*bounded* functions. Pollack's paper could not be obtained in this session, so `bsd.kato_bound.SUPERSINGULAR_PROOF`
proves what is used from scratch, in four steps.
1. From the Hecke relation with $a_p=0$, the Mazur–Tate elements satisfy $\theta_n\equiv-\Phi_{p^{n-1}}(1+T)\,\theta_{n-2}$
   modulo $(1+T)^{p^{n-1}}-1$. So $\theta_{2k+1}$ vanishes at every $\zeta-1$ with $\zeta$ of order $p^t$, $t$ even, $2\le t\le2k$.
2. Hence $L^-:=G/\log^+$ is the limit of the *integral* polynomials $(-1)^{k+1}\theta_{2k+1}/\prod_{t\text{ even}\le2k}\Phi_{p^t}(1+T)$.
   Here $L_\alpha=G+\alpha H$ is the decomposition SW §3.5 uses to define the $D_p$-valued $L_p$, and
   $\log^+=\tfrac1p\prod_{t\text{ even}}\Phi_{p^t}(1+T)/p$ is Pollack's half-logarithm, with constant term $1/p$. So
   $\operatorname{ord}_TL_p=\min(\operatorname{ord}_TL^-,\operatorname{ord}_TL^+)$ with $L^\pm\in\mathbb{Z}_p[[T]]$.
3. Interpolation at roots of unity and Weierstrass division give $L^-\equiv-\theta_1$ modulo $((1+T)^p-1)\,\mathbb{Z}_p[[T]]$.
   So $[T^j]L^-\equiv-[T^j]\theta_1\pmod p$ for $j<p$. That costs $p^2$ symbol values, the same as an ordinary prime.
4. The algebraic side is SW Thm 7.3/7.4 (stated for every semistable $p$, so including good supersingular $p$, with
   $L_p$ and $f_E$ $D_p$-valued) and Thm 6.2 (Perrin-Riou).

For 389a1 the $T^0$ coefficient vanishes exactly and the $T^1$ coefficient vanishes mod $p$ at all four primes, as
rank 2 forces, while the $T^2$ coefficient is a unit. The method was checked in five ways:
- exactly on test cases: the predicted divisibilities of $\theta_2,\theta_3$ and the level-4 congruence hold for
  37a1 at $p=17,19$ and 11a1 at $p=19,29$;
- by controls: rank 0 (11a1) gives order 0 and rank 1 (37a1) gives order 1;
- by the constant terms: $L^-(0)=(p-1)m'(0)$ and $L^+(0)=2m'(0)$, matching Pollack's interpolation formulas as quoted
  in Skinner's 2018 Arizona Winter School notes;
- against a published example: 1483a1 at $p=5$ (rank 2, SW §12.2). The order is 2 there, and all four $T^2$, $T^3$
  coefficients they display are reproduced, reading their display as $(1-\varphi)^{-2}L_p$ in their period normalisation
  (which is ours divided by 2, as found for their ordinary example 446d1);
- by a second certificate from $L^+$ via $\theta_2$ at $p=107$ and $p=599$.

**Regulator anomalies.** At level 2 the coefficient $c_2$ is known modulo $p$. It was a unit at all 663 good ordinary
primes. If the normalised leading coefficient were a random $p$-adic integer, about $\sum1/p\approx1.56$ exceptions
would be expected in $[5,5000)$. Zero is unremarkable (probability about $e^{-1.56}\approx0.21$), but it means no
prime needed a higher level. At the Mazur-anomalous primes $3, 859, 1213, 1493, 1951, 2521$ ($a_p\equiv1$), the factor
$(1-1/\alpha)^2$ is divisible by $p^2$, yet $c_2$ is still a unit. This is consistent with (not a proof of) the $p$-adic
BSD formula, in which the $p$-adic regulator then carries a compensating $p^{-2}$. As an internal check,
$c_0=0$ exactly and $c_1\equiv0\pmod p$ at all 663 primes, as Kato's theorem and rank 2 force. A random value would
vanish mod $p$ only once in $p$ tries. The GP loop that computes these coefficients agrees exactly (modulo $p^3$) with the
independent Python implementation `kato_bound.riemann_coefficients`, checked at $p\le211$ and at both ends of
the extension, $p=3001$ and $p=4999$.

**Relation to known results, and what is left.** Stein–Wuthrich prove much more than finiteness for this curve.
Their Theorem 1.1 gives $Ш[p]=0$ for every non-CM curve of conductor $\le30000$ and rank $\ge2$ at good ordinary
$5\le p<1000$ with $\bar\rho_p$ surjective. Their Theorem 12.3, specific to 389a1, gives $Ш(389\mathrm{a}1)[p]=0$ for
$p=2$ and for all 5,005 good ordinary primes $p<48859$ except $p=16231$. Recounted here, 5,005 is exactly the number of
good ordinary odd primes below 48859, so $p=3$ is included. *Correction (2026-11-10):* earlier versions of this
paragraph, of file 04 and of the README cited only Theorem 1.1. They said the range $1000<p<5000$ went beyond the
published range. That was wrong: at every good ordinary prime the finiteness proved here is weaker than a published
theorem. Theorem 12.3 excludes exactly $p=389$ and the supersingular primes, which their Remark 12.5 calls too slow for
their implementation. Below 5000 those are 107, 389, 599, 1049 and 2957; 6.7 proves $Ш[p]=0$ there. The cost is about
$p^2$ symbol evaluations per prime (318 CPU-minutes
in total, including the $L^+$ cross-checks; the extension from 3000 to 5000 alone took about four CPU-hours), so the method never reaches all $p$. T2 asks for exactly that, and it needs a mechanism linking the $p$-adic $L$-functions for different
$p$, which nobody has (file 03, §3.4).

Reproduce: `python research/t2-389a1/run_t2.py --pmax 3000`, `python research/t2-389a1/run_t2.py --pmin 3001 --pmax 4999 --tag ext5000` (shard with `--shard i --nshards k`), `python research/t2-389a1/run_supersingular.py`,
`python research/t2-389a1/rho3_surjective.py`, `python research/t2-389a1/aggregate.py`. Per-prime records: `research/t2-389a1/data/`; summary: `report.md` and `summary.json` there.

## 6.7 $Ш(389\mathrm{a}1)[p]=0$ for every prime $p<5000$

**Certified computation (2026-11-10).** *$Ш(389\mathrm{a}1)[p]=0$ for every prime $p<5000$*: all 669 primes, with
$b_p=0$ at every odd one. Combined with Stein–Wuthrich's Theorem 12.3 (cited, see 6.6), the supersingular primes
certified beyond 5000 extend this: the paragraph *Beyond 5000* below states how far.

*What is new.* Relative to the literature checked (Stein–Wuthrich 2013, which lists what it excludes), the new primes
are $107, 389, 599, 1049, 2957$: the split multiplicative prime and the supersingular primes below 5000. At $p=2$ and
at the 663 good ordinary odd primes, the result reproduces their Theorem 12.3 with an independent implementation and
partly different proofs. Proving $Ш[p]=0$ rather than finiteness also proves the cyclotomic Iwasawa main conjecture
for $(389\mathrm{a}1,p)$ at each odd $p<5000$ (corollary below). It is still a statement about one curve and finitely
many primes. It is not progress on T2, which asks for all $p$ at once.

**Method.** This is Stein–Wuthrich's Algorithm 11.1, implemented in `bsd/sha_zero.py` and run by
`research/t2-389a1/run_sha_zero.py`. For odd $p$ it computes $b_p$ with $\#Ш(p)\le p^{b_p}$. Since $\#Ш(p)$ is a
square (Cassels), $b_p\le1$ already forces $Ш(p)=0$. The theorem chain, all as stated in SW and read there:
- Kato (SW Thm 7.3; Prop. 7.2 for $p\ge5$): if $\rho_p$ is surjective and $p$ semistable, then $L_p=f_E\,d$ with
  $d\in\Lambda$, and $T\mid d$ in the split case. So $\operatorname{ord}_p f_E^*\le\operatorname{ord}_pL_p^*$.
- Schneider, Perrin-Riou, Jones (SW Thm 6.1; ordinary and multiplicative $p$) and Perrin-Riou (SW Thm 6.2;
  supersingular $p$, $D_p$-valued). Given $\operatorname{ord}_Tf_E=r$, which T2 (6.6) proved, the leading term of
  $f_E$ has the valuation of $\varepsilon_p\prod c_v\,\#Ш(p)\,\mathrm{Reg}_\gamma/\#E(\mathbb{Q})(p)^2$. Here
  $\varepsilon_p=(1-1/\alpha)^2$, with $\operatorname{ord}_p\varepsilon_p=2\operatorname{ord}_p\#\tilde E(\mathbb{F}_p)$; at
  split $p$ it is replaced by $\mathcal{L}_p/\log_p\kappa(\gamma)$; and
  $\mathrm{Reg}_\gamma=\mathrm{Reg}_p/\log_p(1+p)^r$.

Three steps of the algorithm are replaced by short proofs, given in full in the module docstring:
- **Lemma A** (replaces Serre's and Dickson's criteria for surjectivity of $\bar\rho_p$). A subgroup of
  $GL_2(\mathbb{F}_\ell)$ containing a transvection fixes a line or contains $SL_2(\mathbb{F}_\ell)$. By Tate's
  uniformisation, inertia at the multiplicative prime 389 gives a transvection, because $\ell\nmid\operatorname{ord}_{389}\Delta=1$.
  At $\ell=389$ itself the transvection comes from the Kummer extension $\mathbb{Q}_{389}^{ur}(\zeta_{389},q_E^{1/389})$.
  A Frobenius whose characteristic polynomial is irreducible mod $\ell$ rules out a Borel subgroup. The determinant is
  cyclotomic, so the image is $GL_2(\mathbb{F}_\ell)$. The tests check the group theory by brute force for $\ell=3,5,7$.
  At $p=3$ the 3-adic statement is needed, and comes from `rho3_surjective.py` (6.6).
- **Lemma B** (replaces the regulator computation and the saturation proof at good ordinary $p\nmid\prod c_v$).
  SW (4.1) give $h(mP)=2\log_p(e/\sigma_p(t))$, with $\sigma_p=t+O(t^2)$ in $\mathbb{Z}_p[[t]]$ (SW §4.1, citing
  Mazur–Tate). For a point that lies in the formal group and reduces into the identity component everywhere,
  $e/\sigma_p(t)$ is a $p$-adic unit, so the height is in $p\mathbb{Z}_p$. If $p\nmid\#\tilde E(\mathbb{F}_p)$, this gives
  $\operatorname{ord}_p\mathrm{Reg}_\gamma\ge0$. If $\operatorname{ord}_p\#\tilde E(\mathbb{F}_p)=1$ (the anomalous primes
  3, 859, 1213, 1493, 1951, 2521), a basis adapted to the kernel of reduction gives $\ge-2$, exactly offsetting
  $\operatorname{ord}_p\varepsilon_p=2$. Both bounds hold for *any* basis, so $b_p\le\operatorname{ord}_pL_p^*=0$ without
  knowing $E(\mathbb{Q})$ exactly.
- **Saturation**, needed where PARI's regulator enters the bound: the supersingular primes and 389. If $p$ divided
  $[E(\mathbb{Q}):\langle P_1,P_2\rangle]$, some $a\not\equiv0\pmod p$ would have $a_1P_1+a_2P_2\in pE(\mathbb{Q})$.
  Reduction at a prime $\ell$ with $\operatorname{ord}_p\#E(\mathbb{F}_\ell)=1$ gives a linear form on $\mathbb{F}_p^2$
  that kills $a$. Two independent forms exclude every $a$. If the forms stay degenerate, the common kernel is tested
  with `ellisdivisible` and returned as a witness, so the tests' non-saturated controls fail visibly.

| primes | count | $b_p$ | leading term | regulator |
|---|---|---|---|---|
| 2 | 1 | – | – | 2-descent: $\dim\operatorname{Sel}_2=2=\operatorname{rank}$, $E(\mathbb{Q})[2]=0$, so $Ш[2]=0$ |
| good ordinary, non-anomalous | 657 | 0 | $[T^2]L_p$ a unit (6.6) | Lemma B(i) |
| good ordinary, anomalous | 6 | 0 | $[T^2]L_p$ a unit, $\operatorname{ord}_p\varepsilon_p=2$ | Lemma B(ii) |
| supersingular 107, 599, 1049, 2957 | 4 | 0 | $[T^2]L^-$ a unit, so $(1-\varphi)^{-2}L_p^*$ has valuation 0 | PARI $D_p$-valued, $\operatorname{ord}=2$; generators $p$-saturated |
| 389, split multiplicative | 1 | 0 | $[T^3]L_{389}\equiv321$ | PARI's modified height, $\operatorname{ord}\mathrm{Reg}_\gamma=0$, $\operatorname{ord}\mathcal{L}_{389}=1$; 389-saturated |

At a supersingular $p$, $L_p=G\,\omega-pH\,\varphi(\omega)$ with $G=L^-\log^+$ and $H=L^+\log^-$ (6.6). With
$\varphi^2=-1/p$ one gets
$(1-\varphi)^{-2}L_p^*=\big(((p-1)L_2^-+2pL_2^+)\,\omega+p(2L_2^--(p-1)L_2^+)\,\varphi(\omega)\big)/(p+1)^2$, where
$L_2^\pm=[T^2]L^\pm$. At $T=0$ this formula reproduces SW's interpolation $(1-\varphi)^{-2}L_p(0)=(L(E,1)/\Omega)\,\omega$.
Valuations are taken in the lattice $\mathbb{Z}_p\omega+\mathbb{Z}_p\varphi(\omega)$, which here equals
$\mathbb{Z}_p\omega+\mathbb{Z}_p\eta$ (checked with PARI's Frobenius matrix).

**Conventions.** The normalisations are those of the theorems: $\mathrm{Reg}_\gamma$ at supersingular $p$ (as in SW's
own example 12.2) and $\mathcal{L}_p/\log_p\kappa(\gamma)$ at split $p$ (SW Thm 6.1). Read literally, Algorithm 11.1
prints $\mathrm{Reg}_p$ and $\mathcal{L}_p$. That reading would give the impossible $b_p=-2$ (supersingular) and $-1$
(at 389), and $-1$ for SW's own example 446d1 at $p=223$. The theorem versions also give the larger, weaker bounds.

**Checks against Stein–Wuthrich's worked examples** (`tests/test_sha_zero.py`):
- PARI reproduces their $\mathrm{Reg}_5(446\mathrm{d}1)$ (4.2) digit for digit.
- At 446d1, $p=223$ (§12.1), PARI reproduces $\mathcal{L}_{223}$ and the modified $\mathrm{Reg}_{223}$ digit for digit, and $b_{223}=0$.
- At 1483a1, $p=5$ (§12.2, supersingular), $b_5=0$. The $D_p$-valued regulator matches theirs to $O(5^6)$ in the
  $\varphi(\omega)$-coordinate and to $O(5^4)$ in the $\omega$-coordinate. There PARI has $4\cdot5^4$ at every precision
  tried, where they print $O(5^5)$.
- 446d1 at $p=5$ (§7.1) is anomalous ($\#\tilde E(\mathbb{F}_5)=10$), so $\operatorname{ord}_5\varepsilon_5=2$. SW's remark
  that $\varepsilon_5$ is prime to 5 is therefore inaccurate. Their conclusion still stands, because with
  $\mathrm{Reg}_\gamma$ (valuation $-1$) the count is the same. Lemma B gives $b_5=1$ there, and the regulator gives 0.

These examples matter because two of PARI's regulators first disagreed with SW: `ellrank`'s points for 446d1 have
index 3 in $E(\mathbb{Q})$, which is exactly what the saturation certificate is for.

**Cross-checks, not used in the proof.**
- The generators are $p$-saturated at all 663 ordinary primes as well. So PARI's regulator gives the same $b_p=0$ by a
  second, independent route. Lemma B predicts its valuation exactly (2, or 0 at anomalous $p$), and no prime deviates.
- Direction check. Kato's divisibility and Perrin-Riou's formula make $(1-\varphi)^{-2}L_p^*$ a scalar multiple of the
  $D_p$-valued regulator. Its coordinate ratio, $-2-L_2^+/L_2^-\bmod p$, matches PARI's at $p=107$ (81) and $p=599$
  (413), the two primes where $L^+$ is known. Swapping $L^\pm$ gives 47 and 151, so the check is not vacuous; it
  confirms the conventions of 6.6 independently of SW.
- PARI's `ellsaturation(E, gens, 5000)` returns a basis of the same regulator.

**Corollary (main conjecture).** By the remark after SW Algorithm 11.1, $b_p=0$ means $f_E$ and $L_p$ have the same
order and leading-term valuation, so $L_p/f_E\in\Lambda^\times$ (in $T\Lambda^\times$ at 389). So the cyclotomic
main conjecture holds for $(389\mathrm{a}1,p)$ at every odd $p<5000$. At good ordinary $p$ this is also a case of
Skinner–Urban (SW Thm 7.5, with $q=389$). The supersingular and multiplicative literature (Wan, Sprung, Skinner) was
not checked.

**Beyond 5000.** SW Thm 12.3 leaves out the supersingular primes, the first of which above 5000 are 6661, 8263, 9397,
9551, 14633 and 15101. The same certificate, applied to them, is in `research/t2-389a1/sha_zero_report.md`. If all six are
certified, then $Ш(389\mathrm{a}1)[p]=0$ for every prime $p<16231$: our certificates for $p<5000$ and at those six,
together with SW Thm 12.3 (cited) for the good ordinary primes between. 16231 is the ordinary prime SW could not
finish: $\operatorname{ord}_p\mathrm{Reg}_p=3$ there (their Remark 12.4), so $[T^2]L_p$ is needed modulo $p^2$. That means
level-3 Riemann sums at a cost of $p^3\approx4\cdot10^{12}$, or overconvergent symbols. The remaining supersingular
primes below 48859 are 28591, 30671, 30869, 31799, 34781, 36263 and 45161, about 7 CPU-hours at level 2.

*Trusted:* PARI/GP 2.17.2 (`ellpadicregulator`, `ellpadicheightmatrix`, `ellpadics2`, `ellpadicfrobenius`, Tate
parametrisation, `ellap`, `elllog`, `ellisdivisible`, `bnfcertify`), and the inputs of 6.6. *Cited:* SW Thms 6.1, 6.2,
7.3, Prop. 7.2 and 12.3; the integrality of $\sigma_p$ (Mazur–Tate, as stated by SW); Cassels; Tate's uniformisation.
Reproduce: `python research/t2-389a1/run_sha_zero.py --shard i --nshards 4` for $i=0,\dots,3$ (11 CPU-minutes), then
`--aggregate`. Per-prime records: `data/sha_zero_shard*.jsonl`; report: `sha_zero_report.md`, `sha_zero_summary.json`.

## 6.8 Reproduce

```bash
python -m bsd.bsd1 --min 1 --max 19999 --out results/bsd1/bsd1_raw.jsonl --workers 3 --kato-max-conductor 3000
python -c "from bsd import bsd1; bsd1.export('results/bsd1/bsd1_raw.jsonl', 'results/bsd1')"
python -m pytest -q tests/test_bsd1.py tests/test_kato_bound.py tests/test_rigorous_l.py
```

Tracked outputs: `results/bsd1/classes.csv.gz` (one row per class: leading-coefficient ball, independence primes,
descent bound, Kato prime), `results/bsd1/rank2plus.jsonl.gz` (full certificates for rank $\ge2$), and
`results/bsd1/summary.json`.
