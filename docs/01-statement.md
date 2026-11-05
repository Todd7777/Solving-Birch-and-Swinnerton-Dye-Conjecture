# 1. The conjecture, stated precisely

This file fixes notation and normalisations. Almost every "counterexample" to BSD that a
newcomer produces is a normalisation error (period convention, height convention, non-minimal
model, Tamagawa factor at infinity), so the conventions below are spelled out exactly as the
code in `bsd/` uses them.

## 1.1 Elliptic curves over $\mathbb{Q}$

An elliptic curve $E/\mathbb{Q}$ has a **global minimal Weierstrass model**

$$y^2 + a_1 xy + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6,\qquad a_i\in\mathbb{Z},$$

unique up to the substitutions $x = u^2x' + r,\ y = u^3y' + u^2 s x' + t$ with $u=\pm1$.
Its discriminant $\Delta_{\min}$ is the **minimal discriminant**, and
$\omega = dx/(2y + a_1x + a_3)$ is the **Néron differential** (well defined up to sign).
All invariants below are computed on this model; using a non-minimal model rescales the period
and silently breaks the formula.

**Mordell–Weil (Mordell 1922, Weil 1928).** $E(\mathbb{Q}) \cong \mathbb{Z}^r \oplus E(\mathbb{Q})_{\mathrm{tors}}$.
The integer $r = \operatorname{rank} E(\mathbb{Q})$ is the **(algebraic) rank**.
**Mazur (1977):** $E(\mathbb{Q})_{\mathrm{tors}}$ is one of $\mathbb{Z}/n$ ($1\le n\le 10$ or $n=12$)
or $\mathbb{Z}/2\times\mathbb{Z}/2n$ ($1\le n\le 4$).

## 1.2 Reduction, conductor, $L$-function

For a prime $p$, reduce the minimal model modulo $p$. Either the reduction is an elliptic curve over
$\mathbb{F}_p$ (**good** reduction, $p\nmid\Delta_{\min}$), or it has a node (**multiplicative**,
*split* if the tangent slopes are in $\mathbb{F}_p$, *non-split* otherwise), or a cusp (**additive**).
Put

$$a_p = p + 1 - \#\tilde E(\mathbb{F}_p)\ (p \text{ good}),\qquad a_p = 1,\ -1,\ 0\ \text{(split, non-split, additive)}.$$

Hasse: $|a_p| \le 2\sqrt p$ for good $p$. The **conductor** is $N = \prod_p p^{f_p}$ with $f_p = 0,\,1,\,\ge 2$
for good, multiplicative, additive reduction ($f_p = 2$ for additive $p\ge 5$; Tate's algorithm computes
everything, including the Kodaira symbol and the Tamagawa number $c_p$ below).

$$L(E,s) = \prod_{p\nmid N}\bigl(1 - a_p p^{-s} + p^{1-2s}\bigr)^{-1}\prod_{p\mid N}\bigl(1 - a_p p^{-s}\bigr)^{-1} = \sum_{n\ge1} a_n n^{-s},$$

absolutely convergent for $\operatorname{Re} s > 3/2$. The point $s=1$ is *outside* the region of convergence:
the conjecture is about the analytic continuation.

**Modularity (Wiles 1995, Taylor–Wiles 1995, Breuil–Conrad–Diamond–Taylor 2001).** There is a
newform $f = \sum a_n q^n \in S_2(\Gamma_0(N))$ with $L(E,s) = L(f,s)$. Consequently
$\Lambda(E,s) := N^{s/2}(2\pi)^{-s}\Gamma(s)L(E,s)$ is entire and

$$\Lambda(E,s) = w(E)\,\Lambda(E,2-s),\qquad w(E)\in\{\pm1\}$$

the **root number**, a product of computable local root numbers. Thus $r_{\mathrm{an}} := \operatorname{ord}_{s=1}L(E,s)$,
the **analytic rank**, is well defined for every $E/\mathbb{Q}$, and $r_{\mathrm{an}}$ is even iff $w(E)=+1$.

## 1.3 The invariants on the right-hand side

* **Real period.** Let $\Lambda_E = \{\int_\gamma \omega\}\subset\mathbb{C}$ be the period lattice and
  $\Omega^+$ its least positive real period. $E(\mathbb{R})$ has $c_\infty = 2$ connected components if
  $\Delta_{\min} > 0$ and $c_\infty = 1$ if $\Delta_{\min} < 0$. The period in BSD is

  $$\Omega_E := \int_{E(\mathbb{R})}|\omega| = c_\infty\,\Omega^+ .$$

  (PARI's `E.omega[1]` is $\Omega^+$; PARI's `elltamagawa` multiplies the finite Tamagawa product by
  $c_\infty$; Cremona's tables list $\Omega_E$. The code uses $\Omega_E$ and the *finite* Tamagawa product.)
* **Tamagawa numbers.** $c_p := [E(\mathbb{Q}_p) : E_0(\mathbb{Q}_p)]$, where $E_0(\mathbb{Q}_p)$ is the set of points
  with non-singular reduction. $c_p = 1$ for good $p$; $c_p = \operatorname{ord}_p\Delta_{\min}$ for split
  multiplicative $p$; $c_p\in\{1,2\}$ for non-split; $c_p\in\{1,2,3,4\}$ for additive.
* **Néron–Tate height and regulator.** For $P\in E(\mathbb{Q})$ with $x(P) = a/b$ in lowest terms let
  $h(P) = \log\max(|a|,|b|)$ and

  $$\hat h(P) := \lim_{n\to\infty} 4^{-n} h(2^nP).$$

  This is the normalisation of Cremona's book, PARI, LMFDB and Sage; **Silverman's book uses $\hat h/2$**.
  $\hat h$ is a positive definite quadratic form on $E(\mathbb{Q})/\mathrm{tors}$ with bilinear form
  $\langle P,Q\rangle = \tfrac12(\hat h(P+Q) - \hat h(P) - \hat h(Q))$. The **regulator** is
  $\operatorname{Reg}(E/\mathbb{Q}) = \det(\langle P_i,P_j\rangle)$ for a $\mathbb{Z}$-basis $P_1,\dots,P_r$ of
  $E(\mathbb{Q})/\mathrm{tors}$ (and $\operatorname{Reg}=1$ if $r=0$). With Silverman's normalisation the regulator
  is $2^{-r}$ times this, and the formula below must be modified accordingly.
* **Tate–Shafarevich group.**

  $$Ш(E/\mathbb{Q}) := \ker\Bigl(H^1(\mathbb{Q},E) \to \prod_v H^1(\mathbb{Q}_v,E)\Bigr),$$

  the group of isomorphism classes of torsors under $E$ that have points over every completion
  $\mathbb{Q}_v$ (so they are "everywhere locally trivial" but may have no rational point). It is a torsion
  group and $Ш[n]$ is finite for every $n$, because of the exact sequence
  $0\to E(\mathbb{Q})/nE(\mathbb{Q})\to \operatorname{Sel}_n(E/\mathbb{Q})\to Ш(E/\mathbb{Q})[n]\to 0$ and the finiteness of Selmer
  groups. **It is conjectured, and not known in general, that $Ш(E/\mathbb{Q})$ is finite.**
  **Cassels (1962):** if $Ш$ is finite, the Cassels–Tate pairing is non-degenerate and alternating,
  hence $\#Ш$ is a perfect square. (For abelian varieties of higher dimension Poonen–Stoll 1999 showed
  $\#Ш$ can be twice a square; for elliptic curves it is a square.)
* **Torsion.** $\#E(\mathbb{Q})_{\mathrm{tors}}$, computable by Nagell–Lutz or reduction modulo small primes.

## 1.4 The conjecture

**Conjecture (Birch–Swinnerton-Dyer, 1963–65).** Let $E/\mathbb{Q}$ be an elliptic curve.

**(I) Rank.** $\operatorname{ord}_{s=1}L(E,s) = \operatorname{rank} E(\mathbb{Q})$.

**(II) Leading term.** With $r = \operatorname{rank}E(\mathbb{Q})$,

$$\frac{L^{(r)}(E,1)}{r!} = \frac{\Omega_E\cdot\operatorname{Reg}(E/\mathbb{Q})\cdot\prod_{p\mid N}c_p\cdot\#Ш(E/\mathbb{Q})}{(\#E(\mathbb{Q})_{\mathrm{tors}})^2}.$$

Part (II) presupposes that $Ш(E/\mathbb{Q})$ is finite; the *Clay Millennium problem* statement
(Wiles's official description) asks for (I), and the full conjecture is (I) + (II).

Given everything except $Ш$ one defines the **analytic order of $Ш$**

$$\#Ш_{\mathrm{an}} := \frac{L^{(r)}(E,1)/r!}{\Omega_E\operatorname{Reg}(E/\mathbb{Q})\prod c_p/(\#E(\mathbb{Q})_{\mathrm{tors}})^2},$$

which the conjecture predicts to be $\#Ш$. A computation "verifies BSD numerically" when $\#Ш_{\mathrm{an}}$ is
an integer to high precision, is a perfect square, and its $p$-parts agree with whatever descents can
say about $Ш[p]$. This is what `python -m bsd verify` does; it is evidence, not proof.

### Remarks

1. **Isogeny invariance (Cassels 1965; Tate 1966; Milne).** If $E\to E'$ is a $\mathbb{Q}$-isogeny then
   $L(E,s) = L(E',s)$, the ranks agree, and the right-hand side of (II) is the same for $E$ and $E'$
   (each factor changes, the product does not), assuming $Ш$ finite. The conjecture is really about isogeny classes.
2. **Parity.** (I) implies the **parity conjecture** $w(E) = (-1)^{r}$. The weaker **$p$-parity conjecture**
   $\operatorname{corank}_{\mathbb{Z}_p}\operatorname{Sel}_{p^\infty}(E/\mathbb{Q}) \equiv r_{\mathrm{an}} \pmod 2$ is a theorem for all
   $E/\mathbb{Q}$ and all $p$ (Dokchitser–Dokchitser 2010; Monsky 1996 for $p=2$; Nekovář for ordinary $p$).
3. **Why one cannot even compute the rank.** The $n$-Selmer group is computable, so one can bound
   $r \le \dim\operatorname{Sel}_n - \dim E(\mathbb{Q})[n]$, and search for points to bound $r$ from below.
   If $Ш$ is finite these two bounds eventually meet (Manin 1971), so **finiteness of $Ш$ (hence BSD) implies
   that the rank is computable**. No algorithm is *known* to compute the rank of an arbitrary elliptic curve.
4. **What can be computed exactly.** $L(E,1)/\Omega^+$ is a rational number computable exactly by modular symbols
   (Manin, Cremona); this is how one *certifies* $L(E,1)=0$ (as opposed to observing a value $10^{-20}$).
   Derivatives $L^{(k)}(E,1)$ are computed numerically from rapidly convergent series (about $\sqrt N$ terms).
   Hence: $r_{\mathrm{an}} = 0$ can be certified; $r_{\mathrm{an}} = 1$ can be certified (root number $-1$ gives
   $\operatorname{ord}\ge1$, and $L'(E,1)\neq0$ is a verifiable numerical inequality); $r_{\mathrm{an}}=2$ can be
   certified ($L(E,1) = 0$ exactly, root number $+1$ gives $\operatorname{ord}\ge 2$, $L''(E,1)\ne 0$ numerically);
   **$r_{\mathrm{an}}\ge 3$ can never be certified** by present methods, because vanishing of $L'(E,1)$ or
   $L''(E,1)$ cannot be established numerically. The code records exactly this distinction for every curve.

## 1.5 Variants and refinements

* **Over a number field $K$.** Same shape, with $\Omega_E$ replaced by the product of the local periods at the
  infinite places, an extra factor $|d_K|^{-1/2}$, heights relative to $K$, and $L(E/K,s)$. The analytic
  continuation is known for totally real $K$ only when $E$ is modular (proved for real quadratic fields by
  Freitas–Le Hung–Siksek 2015 and totally real cubic fields by Derickx–Najman–Siksek 2020), and only
  meromorphically over CM fields (potential modularity, Allen–Calegari–Caraiani–Gee–Helm–Le Hung–Newton–Scholze–Taylor–Thorne 2018).
* **Function fields.** For $E$ over $K = \mathbb{F}_q(C)$, $L(E,s)$ is a polynomial in $q^{-s}$ and the same
  conjecture is formulated. Here BSD is *equivalent* to finiteness of $Ш$ (Tate 1966, Milne 1975; see file 02).
* **$p$-adic BSD (Mazur–Tate–Teitelbaum 1986).** The cyclotomic $p$-adic $L$-function $L_p(E,T)$ of
  Mazur–Swinnerton-Dyer should vanish at $T=0$ to order $r$ ($r+1$ if $p$ is split multiplicative) with leading
  term $(1-\alpha^{-1})^2\operatorname{Reg}_p\prod c_p\#Ш/\#\mathrm{tors}^2$, $\operatorname{Reg}_p$ the regulator of the
  $p$-adic height pairing. `python -m bsd padic` checks this numerically.
* **Bloch–Kato.** BSD is the Tamagawa number conjecture for the motive $h^1(E)(1)$; the equivariant version
  (Burns–Flach) also encodes the Mazur–Tate refined conjectures.
* **Goldfeld's conjecture (1979).** Among quadratic twists $E^D$, half have rank 0 and half rank 1
  (average analytic rank $1/2$). **Minimalist conjecture:** the average rank of all $E/\mathbb{Q}$ ordered by height
  is $1/2$ (Bhargava–Shankar prove $\le 0.885$).

## 1.6 The original heuristic (1965)

Birch and Swinnerton-Dyer computed on EDSAC-2 the partial products
$f_E(X) = \prod_{p\le X,\ p\nmid N} N_p/p$ ($N_p = \#\tilde E(\mathbb{F}_p)$) and observed
$f_E(X)\approx C(\log X)^{r}$. Heuristically, $N_p/p$ is the $p$-th Euler factor of $L(E,s)$ at $s=1$,
so the product "is" $1/L(E,1)$ and should blow up like $(\log X)^r$ when $L$ vanishes to order $r$.
This is not a theorem (the Euler product does not converge at $s=1$); Goldfeld (1982) proved that **if**
$f_E(X)\sim C(\log X)^r$ then $r = r_{\mathrm{an}}$, and $C$ equals $L^{(r)}(E,1)/r!$ times $\sqrt2 e^{r\gamma}$,
the $\sqrt2$ being explained by Conrad (2005) as a general feature of partial Euler products on the critical line.
`python -m bsd heuristic` reproduces the experiment; `results/heuristic.png` shows slopes $\approx 0,1,2,3$ for
the curves 11a1, 37a1, 389a1, 5077a1.
