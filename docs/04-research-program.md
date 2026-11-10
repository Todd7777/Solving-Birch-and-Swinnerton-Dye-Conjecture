# 4. An honest research program

This file is what I would hand to someone who wants to *work on* BSD rather than read about it. It
contains no claim of progress. It (i) decomposes the conjecture into sub-statements with their current
status, (ii) lists concrete, falsifiable intermediate targets ordered by plausibility, (iii) records the
standard failure modes of amateur "proofs", and (iv) says what the computational part of this repository can
and cannot contribute.

## 4.1 Decomposition into sub-conjectures

For a fixed $E/\mathbb{Q}$ and prime $p$, with $r=\operatorname{rank}E(\mathbb{Q})$, $r_{\mathrm{an}}=\operatorname{ord}_{s=1}L(E,s)$, $s_p=\operatorname{corank}\operatorname{Sel}_{p^\infty}$:

| | Statement | Status |
|---|---|---|
| **S1** | $Ш[p^\infty]$ finite (⇔ $s_p = r$) | theorem if $r_{\mathrm{an}}\le1$; for $r\ge2$ provable *prime by prime* via Kato + $p$-adic $L$ (file 02 F4), never for all $p$ |
| **S2** | $s_p \equiv r_{\mathrm{an}} \pmod 2$ | **theorem** for all $E,p$ (Dokchitser–Dokchitser) |
| **S3** | $s_p \le r_{\mathrm{an}}$ | theorem for $r_{\mathrm{an}}\le1$; open otherwise (needs a rank-$r$ Euler system) |
| **S4** | $s_p \ge r_{\mathrm{an}}$ | theorem for $s_p\le1$ under hypotheses (converse theorems); open otherwise |
| **S5** | $\operatorname{ord}_{T=0}L_p(E,T) = s_p$ (+1 if split mult.) | ⇐ main conjecture (theorem in great generality) + Schneider non-degeneracy (open) |
| **S6** | $\operatorname{ord}_{T=0}L_p(E,T) = r_{\mathrm{an}}$ | open for $r_{\mathrm{an}}\ge2$; no mechanism |
| **S7** | $p$-part of formula (II) | theorem for $r_{\mathrm{an}}\le1$ under hypotheses (B5, B7, B8), open for $r\ge2$ |

BSD (I) for $E$ ⇔ S3 ∧ S4 ∧ S1 (for one $p$). BSD (II) ⇔ S7 for all $p$ (with S1 for all $p$).

**The minimal missing piece is S3 for $r_{\mathrm{an}}=2$**: *"if $L(E,1)=L'(E,1)=0\ne L''(E,1)$ then $\operatorname{corank}\operatorname{Sel}_{p^\infty}(E/\mathbb{Q})\le2$."*
Proving this for a single prime $p$ and a single curve (e.g. 389a1) by a method that does not amount to computing $\operatorname{Sel}_{p^\infty}$ directly
would be the first genuinely new bridge since Kolyvagin.

## 4.2 Intermediate targets, in increasing order of ambition

Each of these is a well-posed theorem to aim at; each is open; none is known to be out of reach.

**T1. Non-vanishing of $p$-adic heights (Schneider's conjecture) for non-CM curves of rank 1.**
Equivalent to a $p$-adic transcendence statement for a single point. Would close the gap between the main conjecture and the $p$-adic rank conjecture in rank 1 for all good ordinary $p$. Tools: $p$-adic Gross–Zagier (the height of the Heegner point is $L_p'(E,0)$, so non-vanishing of the $p$-adic height is equivalent to $L_p'(E,0)\ne0$ — a statement about a $p$-adic $L$-function that might be attacked by $p$-adic analytic methods or by the "Iwasawa-theoretic non-vanishing" results of the kind Kato/Rohrlich prove for values).

**T2. Finiteness of $Ш(E/\mathbb{Q})$ for one curve of rank 2 (e.g. 389a1) at all primes simultaneously.**
Needs: Kato's divisibility at all good ordinary $p$ (known), supersingular analogue (Kobayashi/Sprung — known in many cases), *and* a uniform way to certify $\operatorname{ord}_{T=0}L_p(E,T)\le 2$ for all $p$ at once. The last point is the obstacle: it is a non-vanishing statement for infinitely many $p$-adic $L$-functions (of the second derivative). An analogue of Rohrlich's theorem "for derivatives" would do it. **This is the single most concrete open target that a determined research effort could plausibly attack**, and success would give the first rank-2 curve with finite $Ш$ — and, by Dokchitser parity, S2, nothing more; BSD (I) for that curve already holds by direct verification.

**T3. Gross–Zagier for Kolyvagin classes.**
Formulate and prove a formula expressing $L^{(M+1)}(E/K,1)$, where $M+1$ is the Selmer corank, in terms of Kolyvagin's derived classes $\kappa_n$ with $\nu(n)=M$ (file 03, §3.3). Even a *conjectural* formula with numerical support (the $\kappa_n$ are computable: Jetchev–Lauter–Stein 2009) would be a real step; computing enough examples to guess the shape of the constant is a feasible project with the tools in this repository plus Sage's Heegner-point code.

**T4. A number-field analogue of shtukas with two legs**, or any geometric object over $\mathbb{Z}$ whose "intersection numbers" are second derivatives of $L$-functions. Programmatic; this is where plectic cohomology (Nekovář–Scholl) and arithmetic-shtuka ideas live.

**T5. A rank-2 Euler system for $E$** (a non-zero element of $\bigwedge^2H^1(\mathbb{Q},T_pE)$ with Euler-system norm relations whose regulator is $L''(E,1)$). Equivalent, by the Mazur–Rubin/Burns–Sano formalism, to S3 for $r_{\mathrm{an}}=2$ together with explicit bounds on $Ш$.

## 4.3 Standard failure modes (read before writing a proof)

Every one of the following has appeared in attempted proofs. If an argument uses one of them, it is wrong.

1. **Evaluating the Euler product at $s=1$.** $\prod_p(1-a_pp^{-1}+p^{-1})^{-1}$ does not converge; the heuristic $\prod_{p\le X}N_p/p\sim C(\log X)^r$ is a conjecture (file 01 §1.6), and even if it held it gives $r_{\mathrm{an}}$, not $r$.
2. **Assuming $Ш$ is finite** (explicitly, or implicitly by writing $\#Ш$ in an identity, or by assuming $\operatorname{Sel}_n$ has the "expected" size).
3. **Assuming the rank is computable**, or that a descent that fails to find points proves there are none (it does not: that is exactly the $Ш[2]\ne0$ situation of 571a1).
4. **Confusing $r_{\mathrm{an}}$ with a numerically observed vanishing.** $L'(E,1)=10^{-40}$ proves nothing. Only $L(E,1)=0$ can be certified (modular symbols).
5. **Normalisation errors**: $\Omega^+$ vs $\Omega_E$; Silverman's $\hat h$ vs Cremona's; a non-minimal model; forgetting $c_\infty$; using $L^{(r)}(E,1)$ instead of $L^{(r)}(E,1)/r!$. The code in this repository exists partly to make these impossible to get wrong.
6. **Transferring Tate's function-field theorem to $\mathbb{Q}$** ("BSD ⇔ finiteness of $Ш$"): no such equivalence is known over number fields (file 03 §3.6).
7. **Using the parity conjecture as if it were the rank conjecture.** $w(E)=(-1)^r$ is implied by BSD and is itself unproven in general (only the $p$-parity statement for Selmer coranks is known).
8. **Arguing from the functional equation alone** (e.g. "the order of vanishing is forced by symmetry"). The functional equation controls parity only.
9. **Applying Gross–Zagier/Kolyvagin outside rank 1**, e.g. "take the Heegner point for a rank-2 curve": it is torsion, and Kolyvagin's theorem then says nothing (the derived classes do, but there is no formula — file 03 §3.3).
10. **Treating "$\ge 66.48\%$ of curves satisfy BSD" as evidence of a method.** Those results are statistical consequences of the rank-$\le1$ theorems; they cannot be refined to density $1$ by the same methods, let alone to every curve.

## 4.4 What this repository's computations contribute, and what they do not

**They do:**
* verify formula (II) numerically to $\sim10^{-9}$ relative precision for every curve of conductor $\le 10\,000$ (and more) in Cremona's database, independently of Cremona's own software (PARI vs eclib), with $\#Ш_{\mathrm{an}}$ always an integer square — see `docs/05-numerical-evidence.md`;
* certify, for each curve, exactly which part of the conclusion is a *theorem* (file 02) and which is *numerics*: e.g. for 389a1, $r=r_{\mathrm{an}}=2$ is proven, $Ш[p^\infty]$ is finite for every odd $p<100$, and $Ш$ itself is of unknown finiteness;
* reproduce the 1965 heuristic, the $p$-adic BSD checks, and the hypothesis checks for the Skinner–Urban / Jetchev–Skinner–Wan theorems at each prime;
* with Sage, re-run the Kolyvagin/Kato-based proof machinery (`prove_BSD`) on every curve of conductor $\le1000$, so that for those curves of analytic rank $\le1$ the *full* formula (II) is a theorem whenever the routine returns no exceptional primes.

**They do not:**
* prove anything for a curve of rank $\ge2$ beyond prime-by-prime finiteness of $Ш[p^\infty]$;
* prove formula (II) at any prime for rank $\ge2$;
* substitute for a provable saturation bound unless Sage is available (`--sage` uses Sage's provable height lower bound and index bound; without it the bases are saturated at all primes $<1000$ and cross-checked against Cremona's independently saturated bases);
* *(superseded 2026-11-09)* ~~use interval arithmetic~~: `bsd/rigorous_l.py` now gives Arb enclosures of $L^{(k)}(E,1)/k!$ with proven tail
  bounds, and `bsd/bsd1.py` uses them to certify rank $=$ analytic rank for every curve of conductor $<20000$ (file 06).
  The PARI-based sweeps of file 05 still use floating-point non-vanishing.

## 4.5 Progress log

**2026-11-09.** Nothing here moves the conjecture itself. The computational frontier moved as follows (details in file 06):

* **BSD (I) is certified for every curve of conductor $<20000$** (ranks 0–3, 81,217 isogeny classes), with interval
  arithmetic for the analytic side and an unconditional 2-descent for the rank. This re-derives Cremona's ranks
  independently; it is not a new rank.
* **T2 for 389a1.** $Ш[p^\infty]$ is now proven finite at every good ordinary prime $5\le p<5000$, at $p=389$
  (split multiplicative, exceptional zero), at $p=2$ (where it is 0), and at $p=3$. For $p=3$ one first has to prove
  $\rho_3$ surjective, which uses Elkies' classification of curves whose 3-adic image is surjective mod 3 but not
  mod 9. The four supersingular primes $107, 599, 1049, 2957$ are covered too, so the statement holds for
  **every prime $p<5000$** (669 primes; extended from 3000 on 2026-11-10). They use a self-contained proof of the congruence $L^-\equiv-\theta_1\bmod((1+T)^p-1)$
  for Pollack's bounded $L^-$, which makes supersingular primes as cheap as ordinary ones (section 6.6). At the good
  ordinary primes this is weaker than a published theorem: Stein–Wuthrich's Theorem 12.3 gives $Ш(389\mathrm{a}1)[p]=0$
  at $p=2$ and at every good ordinary $p<48859$ except 16231. *Correction (2026-11-10):* this item used to say the
  work extended Stein–Wuthrich beyond $p=1000$; only their Theorem 1.1 had been taken into account. The per-prime cost
  grows like $p^2$ at level 2, so a computation can never reach "all $p$". The obstruction of T2 is untouched: no
  mechanism links the $p$-adic $L$-functions for different $p$.
* **$Ш(389\mathrm{a}1)[p]=0$ for every prime $p<5000$** (2026-11-10, section 6.7). This is Stein–Wuthrich's
  Algorithm 11.1 at every odd $p<5000$, with $b_p=0$ throughout. Beyond their Theorem 12.3, which excludes them, are 389
  and the supersingular primes 107, 599, 1049, 2957 (other literature was not checked); at the ordinary primes it
  reproduces that theorem independently.
  Three parts of the algorithm are replaced by short proofs. Lemma A proves surjectivity of $\bar\rho_p$ from a
  transvection in inertia at 389 plus one Frobenius. Lemma B uses integrality of $p$-adic heights, which makes the
  regulator and the Mordell–Weil basis unnecessary at ordinary $p$. A saturation certificate by reduction maps covers
  the primes where PARI's regulator is used. Consequences: the cyclotomic main conjecture for $(389\mathrm{a}1,p)$ at
  every odd $p<5000$ (at ordinary $p$ also Skinner–Urban). And the supersingular primes above 5000 can be certified the
  same way; with Theorem 12.3, all six below 16231 would give $Ш[p]=0$ for every $p<16231$ (status in 6.7). Checked against
  all four of SW's worked examples. Two of PARI's regulators first disagreed with SW because `ellrank`'s points for
  446d1 have index 3, which the saturation certificate detects.
* **T1-type statements as a by-product.** By Schneider/Perrin-Riou (Stein–Wuthrich Thm 6.1),
  $\operatorname{ord}_TL_p=r$ also proves that the cyclotomic $p$-adic height pairing on $E(\mathbb{Q})$ is
  non-degenerate. So the computations above prove Schneider's conjecture for the rank-2 curve 389a1 at each of
  those primes, and for 243 rank-2 curves at one prime each. They are instances, not the general statement T1 asks for.
* **Corrected claim.** Earlier versions of files 01, 02 and the README said that analytic rank 3 cannot be certified
  for any curve. That was wrong: Buhler–Gross–Zagier certified it for 5077a1 in 1985. The true frontier is analytic
  rank $\ge4$, where $L''(E,1)=0$ (root number $+1$) or $L'''(E,1)=0$ (root number $-1$) would be needed.

## 4.6 If I had to bet

The conjecture will not be settled by a trick. Among the targets above, **T2** (full finiteness of $Ш$ for one rank-2 curve) and **T1** (Schneider for rank 1) look like decade-scale problems that existing technology could conceivably reach with a new non-vanishing idea for derivatives of $p$-adic $L$-functions. **T3–T5** need a new object. The function-field theory (Yun–Zhang) says such objects *can* exist in principle; the question of whether anything like "two legs" can live over $\operatorname{Spec}\mathbb{Z}$ is, in my judgement, the real content of the Birch and Swinnerton-Dyer conjecture beyond rank one.
