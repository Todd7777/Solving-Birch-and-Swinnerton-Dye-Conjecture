# 2. What is actually proven (a theorem ledger)

Each entry states a theorem with its hypotheses, as precisely as I can vouch for. Where a hypothesis
list is long and technical I say so and point to the paper rather than paraphrase loosely. The entries
are grouped by what they tell us about the two halves of the conjecture: **(I)** rank equality,
**(II)** the leading-term formula, and the prerequisite **finiteness of $Ш$**.

Notation: $E/\mathbb{Q}$ of conductor $N$, $r = \operatorname{rank}E(\mathbb{Q})$, $r_{\mathrm{an}} = \operatorname{ord}_{s=1}L(E,s)$,
$\bar\rho_{E,p}: G_{\mathbb{Q}}\to \mathrm{GL}_2(\mathbb{F}_p)$ the mod-$p$ Galois representation on $E[p]$,
$\operatorname{Sel}_{p^\infty}(E/\mathbb{Q})$ the $p^\infty$-Selmer group, which sits in
$0\to E(\mathbb{Q})\otimes\mathbb{Q}_p/\mathbb{Z}_p\to\operatorname{Sel}_{p^\infty}\to Ш[p^\infty]\to0$, so that
$\operatorname{corank}_{\mathbb{Z}_p}\operatorname{Sel}_{p^\infty} = r + \operatorname{corank}Ш[p^\infty]$ and
**$Ш[p^\infty]$ is finite iff $\operatorname{corank}\operatorname{Sel}_{p^\infty} = r$.**

---

## A. Unconditional structure

**A1 (Mordell–Weil).** $E(\mathbb{Q})$ is finitely generated. **(Mazur 1977)** The torsion is one of 15 groups.

**A2 (Modularity: Wiles, Taylor–Wiles 1995; Breuil–Conrad–Diamond–Taylor 2001).** $L(E,s)$ is entire with
functional equation, so $r_{\mathrm{an}}$ exists and the conjecture is meaningful for every $E/\mathbb{Q}$. Before 2001 the
conjecture could only be *stated* unconditionally for modular curves.

**A3 (Cassels 1962).** If $Ш(E/\mathbb{Q})$ is finite then $\#Ш$ is a perfect square. $Ш[n]$ is finite for every $n$.

**A4 (Isogeny invariance: Cassels 1965, Tate 1966, Milne).** For $\mathbb{Q}$-isogenous $E, E'$ with finite $Ш$:
(II) holds for $E$ iff it holds for $E'$.

**A5 ($p$-parity: Dokchitser–Dokchitser 2010; Monsky 1996 for $p=2$; Nekovář 2006 for ordinary $p$).** For every
$E/\mathbb{Q}$ and every prime $p$: $\operatorname{corank}_{\mathbb{Z}_p}\operatorname{Sel}_{p^\infty}(E/\mathbb{Q})\equiv r_{\mathrm{an}}\pmod 2$.
*Consequence:* if $Ш[p^\infty]$ is finite for a single $p$, then $r\equiv r_{\mathrm{an}}\pmod 2$ — the parity conjecture holds for $E$.

**A6 (Manin 1971).** If $Ш(E/\mathbb{Q})$ is finite then $E(\mathbb{Q})$ is computable by a terminating algorithm
(descent for an upper bound, point search for a lower bound). Unconditionally, no algorithm is known to compute $r$.

**A7 (Kato 2004 + Rohrlich 1984: Mazur's conjecture).** For every $E/\mathbb{Q}$ and every prime $p$, $E(\mathbb{Q}(\mu_{p^\infty}))$
is finitely generated. (Rohrlich: $L(E,\chi,1)\ne0$ for all but finitely many characters $\chi$ of $p$-power conductor;
Kato's Euler system converts non-vanishing into finiteness of the $\chi$-parts of the Selmer group.)

---

## B. Analytic rank 0 and 1 — where the conjecture is a theorem

**B1 (Coates–Wiles 1977).** If $E/\mathbb{Q}$ has complex multiplication and $L(E,1)\ne0$ then $E(\mathbb{Q})$ is finite.
(Stated for CM by the maximal order of a class-number-one field; every CM curve over $\mathbb{Q}$ is $\mathbb{Q}$-isogenous to such a curve, and rank is an isogeny invariant.)

**B2 (Gross–Zagier 1986).** Let $K$ be an imaginary quadratic field in which every prime dividing $N$ splits
(*Heegner hypothesis*), and $y_K\in E(K)$ the Heegner point obtained from a modular parametrisation $X_0(N)\to E$. Then
$$L'(E/K,1) = c_{E,K}\cdot\hat h(y_K)$$
with an explicit *non-zero* constant $c_{E,K}$ (involving the Petersson norm of $f$, $\sqrt{|d_K|}$, the number of units of $K$ and the degree of the parametrisation).
In particular $\operatorname{ord}_{s=1}L(E/K,s)=1 \iff y_K$ has infinite order.

**B3 (Kolyvagin 1988, 1990).** If $y_K$ has infinite order, then $\operatorname{rank}E(K) = 1$ and $Ш(E/K)$ is finite.
Moreover, for all primes $p$ outside an explicit finite set (roughly: $p\nmid 2\,\#E(K)_{\mathrm{tors}}\prod c_\ell$ with
$\bar\rho_{E,p}$ surjective), $\#Ш(E/K)[p^\infty]$ divides $[E(K):\mathbb{Z}y_K]_p^2$ — an explicit upper bound that
is the engine of all "exact" verifications (section B9).

**B4 (Non-vanishing of twists: Waldspurger 1985; Bump–Friedberg–Hoffstein 1990; Murty–Murty 1991).** For any $E/\mathbb{Q}$
there are infinitely many fundamental discriminants $D<0$ satisfying the Heegner hypothesis with
$\operatorname{ord}_{s=1}L(E/\mathbb{Q}(\sqrt D),s) = \operatorname{ord}_{s=1}L(E,s)+\operatorname{ord}_{s=1}L(E^D,s)$ equal to $1$,
whenever $r_{\mathrm{an}}\le 1$.

**Theorem (B2 + B3 + B4).** *If $r_{\mathrm{an}}\le1$ then $r = r_{\mathrm{an}}$ and $Ш(E/\mathbb{Q})$ is finite.*
Equivalently, for every $E/\mathbb{Q}$: $\ r_{\mathrm{an}}\ge\min(r,2)$ — the only general inequality between $r$ and $r_{\mathrm{an}}$ known.

**B5 (Rubin 1987, 1991; Perrin-Riou 1987; later Pollack–Rubin 2004 for supersingular $p$).** For CM curves with
$r_{\mathrm{an}}\le1$: $Ш$ is finite and the $p$-part of (II) holds for all primes $p\nmid 6$ (ordinary and supersingular cases
by the two-variable main conjecture for imaginary quadratic fields; some further exclusions in the rank-1 case).

**B6 (Kato 2004).** Using the Euler system of Beilinson–Kato elements:
(i) if $L(E,1)\ne0$ then $E(\mathbb{Q})$ and $Ш(E/\mathbb{Q})$ are finite (a second proof of the rank-0 case, independent of Heegner points);
(ii) for good ordinary $p$, $X_\infty = \operatorname{Sel}_{p^\infty}(E/\mathbb{Q}_\infty)^\vee$ is $\Lambda$-torsion and
$\operatorname{char}_\Lambda(X_\infty)$ divides the $p$-adic $L$-function $L_p(E,T)$ in $\Lambda\otimes\mathbb{Q}_p$, integrally when
$\bar\rho_{E,p}$ is surjective (Thm 17.4); (iii) for $L(E,1)\ne0$ this gives the explicit bound
$\#Ш[p^\infty]\le\bigl|L(E,1)/\Omega_E\bigr|_p^{-1}\cdot(\text{Tamagawa/torsion correction})$ for such $p$.

**B7 (Skinner–Urban 2014, Thm 2).** Let $p\ge3$ with $E$ good ordinary at $p$, $\bar\rho_{E,p}$ irreducible, and assume
there is a prime $q\ne p$ with $q\,\|\,N$ at which $\bar\rho_{E,p}$ is ramified. Then the cyclotomic Iwasawa main conjecture
holds for $(E,p)$ (the divisibility opposite to Kato's), and consequently **if $L(E,1)\ne0$ the $p$-part of formula (II) holds**:
$\operatorname{ord}_p\bigl(L(E,1)/\Omega_E\bigr) = \operatorname{ord}_p\bigl(\#Ш\prod c_\ell/\#E(\mathbb{Q})_{\mathrm{tors}}^2\bigr)$.
Further primes: Skinner 2016 (multiplicative $p$), Wan / Sprung / Castella–Çiperiani–Skinner–Sprung (supersingular $p$),
Castella–Grossi–Lee–Skinner 2022 (Eisenstein primes, where $\bar\rho_{E,p}$ is reducible).

**B8 (Jetchev–Skinner–Wan 2017; Castella 2018).** Let $E$ be semistable with $r_{\mathrm{an}}=1$. Then the $p$-part of (II)
holds for every prime $p\ge5$ of good reduction with $E[p]$ irreducible (also $p=3$ if $E$ is ordinary at 3, or supersingular
with $a_3=0$); Castella extends this to primes $p>3$ of multiplicative reduction with $E[p]$ irreducible.

**B9 (Converse theorems: Skinner 2020; W. Zhang 2014; Burungale–Tian 2020 for CM; Castella–Grossi–Lee–Skinner 2022.)**
Under hypotheses on $p$ (good ordinary, $\bar\rho_{E,p}$ irreducible/surjective, ramification conditions at primes dividing $N$):
if $\operatorname{corank}_{\mathbb{Z}_p}\operatorname{Sel}_{p^\infty}(E/\mathbb{Q}) = 1$ then $r_{\mathrm{an}} = 1$; and (Skinner–Urban) if the corank is 0 then $L(E,1)\ne0$.
So for a curve of algebraic rank $\le1$, *finiteness of $Ш[p^\infty]$ for one suitable $p$ implies BSD (I) and finiteness of all of $Ш$*.
Zhang's paper also proves **Kolyvagin's conjecture** (some derived Heegner class $\kappa_n$ is non-zero) under similar hypotheses — see file 03.

**B10 (Exact verification for individual curves: Grigorov–Jorza–Patrikis–Stein–Tarniţă 2009; Miller 2011; Creutz–Miller 2012).**
Combining B3, B6, B7 with explicit 2-, 3- and 4-descents and Heegner-index computations, the **full** formula (II) — every prime —
is a theorem for the great majority of curves of conductor $\le 5000$ with $r_{\mathrm{an}}\le1$; the remaining $(E,p)$ pairs are listed in those papers.
(Sage's `prove_BSD` implements this machinery.) This is the only sense in which (II) is known "in full" for any curve.

---

## C. Statistical results (what is known for "most" curves)

**C1 (Bhargava–Shankar 2015 ff.).** Ordering $E: y^2 = x^3+Ax+B$ by height, the average size of $\operatorname{Sel}_n$ is $\sigma(n)$
for $n=2,3,4,5$; hence the average rank is $\le 0.885$, and a positive proportion of curves have rank 0 (and, with Dokchitser parity, rank 1 for a positive proportion of the rest).

**C2 (Bhargava–Skinner 2014; Bhargava–Skinner–Zhang 2014).** A positive proportion of $E/\mathbb{Q}$ have rank 1 and satisfy BSD (I);
in total **at least 66.48 % of elliptic curves over $\mathbb{Q}$ satisfy BSD (I) and have finite $Ш$** (ordered by naive height).

**C3 (A. Smith 2017, 2022).** For $E$ with full rational 2-torsion and no rational cyclic 4-isogeny, $100\%$ of quadratic twists have
$2^\infty$-Selmer rank $\le 1$; with parity and the converse theorems this gives Goldfeld's conjecture and BSD (I) for $100\%$ of the twists
(e.g. the congruent number curves $y^2=x^3-n^2x$). Earlier: Heath-Brown 1994, Kane 2013, Klagsbrun–Mazur–Rubin 2013 on Selmer distributions.

None of these results touches any *single* curve of rank $\ge2$; they prove BSD on a set of density $<1$ (or $=1$ within a family).

---

## D. Function fields — the case that is understood

**D1 (Tate 1966, Milne 1975; Kato–Trihan 2003 for the $p$-part).** Let $E$ be an elliptic curve over $K=\mathbb{F}_q(C)$. The following are equivalent:
(i) BSD (I) and (II) for $E/K$; (ii) $Ш(E/K)$ is finite; (iii) $Ш(E/K)[\ell^\infty]$ is finite for *one* prime $\ell$ (including $\ell = p = \operatorname{char}K$);
(iv) the Tate conjecture for the elliptic surface $\mathcal{E}\to C$ (in the form $\operatorname{rank}\mathrm{NS} = \operatorname{ord}_{s=1}\zeta$-pole).

**D2 (Known cases).** BSD holds for isotrivial curves, for curves whose elliptic surface is dominated by a product of curves or by a Fermat surface
(Shioda; **Ulmer 2002**: $y^2+xy=x^3-t^d$ over $\mathbb{F}_p(t)$ has unbounded rank and satisfies BSD), and for elliptic K3 surfaces (Tate conjecture
for K3 surfaces over finite fields: Charles 2013, Madapusi Pera 2015, Kim–Madapusi Pera).

**D3 (Yun–Zhang 2017, 2019: higher Gross–Zagier).** For an unramified cuspidal automorphic representation $\pi$ of $\mathrm{PGL}_2$ over $\mathbb{F}_q(X)$
and a quadratic extension, the $r$-th central derivative $L^{(r)}(\pi_F,1/2)$ is, for every $r\ge0$, an explicit multiple of the self-intersection
number of the $\pi$-isotypic Heegner–Drinfeld cycle on the moduli stack of shtukas with $r$ legs. **This is the only known formula anywhere relating
$r$-th derivatives of $L$-functions to algebraic cycles for arbitrary $r$**; its mechanism (shtukas with $r$ legs) has no number-field analogue.

---

## E. Other number fields

**E1 (Zhang 2001; Yuan–Zhang–Zhang 2013; Nekovář 2007).** For a modular elliptic curve $E$ over a totally real field $F$ with $r_{\mathrm{an}}\le1$,
under a mild hypothesis guaranteeing a Shimura-curve parametrisation ($[F:\mathbb{Q}]$ odd, or $E$ multiplicative at some prime), $\operatorname{rank}E(F)=r_{\mathrm{an}}$
and $Ш(E/F)$ is finite. Modularity is known for all $E$ over real quadratic and totally real cubic fields (Freitas–Le Hung–Siksek 2015, Derickx–Najman–Siksek 2020).

**E2.** Over CM fields only potential modularity is known (meromorphic continuation of $L$); BSD-type results exist for special classes (e.g. base change).

---

## F. $p$-adic BSD and Iwasawa theory

**F1 (Greenberg–Stevens 1993).** The Mazur–Tate–Teitelbaum *exceptional zero conjecture*: if $p$ is split multiplicative,
$L_p'(E,0) = \dfrac{\log_p q_E}{\operatorname{ord}_p q_E}\cdot\dfrac{L(E,1)}{\Omega_E}$, where $q_E$ is the Tate period. A complete theorem.

**F2 (Perrin-Riou 1987; Kobayashi 2013).** $p$-adic Gross–Zagier: $L_p'(E,0)$ (ordinary) resp. the $\pm$ versions (supersingular) equal an explicit
multiple of the *$p$-adic* height of the Heegner point.

**F3 (Main conjectures).** Kato (one divisibility, all good ordinary $p$, rationally), Skinner–Urban (equality, ordinary $p$, hypotheses as in B7),
Rubin (CM), Wan, Sprung, Castella–Çiperiani–Skinner–Sprung (supersingular), Skinner (multiplicative), Castella–Grossi–Lee–Skinner (Eisenstein).
*Consequence (Perrin-Riou, Schneider):* for such $(E,p)$, $\operatorname{ord}_{T=0}L_p(E,T)\ \ge\ \operatorname{corank}\operatorname{Sel}_{p^\infty}(E/\mathbb{Q})$, with equality
iff the cyclotomic $p$-adic height pairing is non-degenerate on the Selmer group (**Schneider's conjecture** — open in general, known for CM curves by Bertrand).
Note that even Kato's one-sided divisibility already gives the inequality $\operatorname{corank}\operatorname{Sel}_{p^\infty}(E/\mathbb{Q}) \le \operatorname{ord}_{T=0}L_p(E,T)$ for good ordinary $p$ (via Mazur's control theorem and the structure theorem for $\Lambda$-modules).

**F4 (Stein–Wuthrich 2013).** Turning F3 into an algorithm: for good ordinary $p\ge3$, computing $L_p(E,T)$ to finite precision certifies an upper bound
$\operatorname{ord}_{T=0}L_p\le r_p$; if $r_p$ equals the number of independent rational points found, then $Ш(E/\mathbb{Q})[p^\infty]$ is finite, and with the
integral form of Kato's divisibility one gets an explicit upper bound for $\#Ш[p^\infty]$. **This works for curves of any rank** and is implemented in `bsd/theorems.py`
(`sha_p_finite_via_kato`); e.g. $Ш(389\mathrm{a}1)[p^\infty]$ is finite for every odd prime $p<100$ (see `results/curated.md`).

---

## G. The frontier: precise statements of what is *not* known

1. **No elliptic curve over $\mathbb{Q}$ of rank $\ge2$ is known to have finite $Ш$.** For individual odd primes $p$ one can prove $Ш[p^\infty]$ finite (F4),
   but no argument handles all $p$ at once, so finiteness of the full group is open for every such curve (389a1 included).
2. **No general theorem gives $r = r_{\mathrm{an}}$ when $r_{\mathrm{an}}\ge2$.** For a *specific* curve with $r_{\mathrm{an}}=2$ one can certify both sides
   (file 01, remark 4): 389a1 has $r=2$ by 2-descent and $r_{\mathrm{an}}=2$ by an exact vanishing plus a numerically non-zero $L''$. But for $r_{\mathrm{an}}\ge3$
   not even individual verification is possible: 5077a1 has $r=3$, and all we can *prove* about $L$ is $\operatorname{ord}\ge1$, odd.
3. **No inequality $r\le r_{\mathrm{an}}$ or $r_{\mathrm{an}}\le r$ is known in general**; the only known relation is $r_{\mathrm{an}}\ge\min(r,2)$ (from B2–B4).
4. **Formula (II) is not a theorem for any curve of rank $\ge2$**, and for rank $\le1$ it is a theorem only curve-by-curve (B10) or prime-by-prime (B5, B7, B8); there is no infinite family for which (II) is known at every prime.
5. **Boundedness of ranks is open.** Elkies (2006) found a curve of rank $\ge28$; Elkies–Klagsbrun (2024) one of rank $\ge29$. The heuristic of
   Park–Poonen–Voight–Wood (2019) predicts that only finitely many $E/\mathbb{Q}$ have rank $>21$.
6. **Computability of the rank is open** (A6): it follows from BSD, and is not known otherwise.
7. **Non-degeneracy of $p$-adic heights (Schneider)** is open even for a single non-torsion point in general; it is the missing link between the
   main conjecture (a theorem in large generality) and the $p$-adic rank statement.

The next file explains *why* every known method stops exactly at these points.
