# 3. Why it is hard: where each method stops, and why

The conjecture equates an **archimedean-analytic** quantity (the order of vanishing and leading term of a
Dirichlet series at a point *outside* its half-plane of convergence) with **Diophantine / Galois-cohomological**
quantities (rational points, Selmer groups, $Ш$). Nothing in the definitions connects them. Every theorem in
file 02 rests on one of a very small number of *bridges*, each of which was built by an explicit construction
from modular forms, and each of which is structurally limited to vanishing order $\le 1$.

## 3.1 The bridges, and the rank at which each one breaks

| Bridge | Input | What it produces | Detects | Why it stops there |
|---|---|---|---|---|
| **Beilinson–Kato elements** (Kato 2004) | $K_2$ of modular curves, Siegel units | classes $z\in H^1(\mathbb{Q},T_pE)$ whose dual-exponential is $L(E,1)/\Omega$ (and twisted values) | $L(E,1)\ne0$ | the class is a single element of $H^1$: an Euler system of **rank 1** in the sense of Mazur–Rubin. It bounds a Selmer group of core rank 1 by the *value* $L(E,1)$, i.e. the rank-0 case. Nothing in the construction sees derivatives. |
| **Heegner points** (Gross–Zagier, Kolyvagin) | CM points on $X_0(N)$ (or Shimura curves) | a point $y_K\in E(K)$; derived classes $\kappa_n\in H^1(K,E[p^M])$ | $L'(E/K,1)\ne0$ | Gross–Zagier is a formula for the **first** derivative because the Heegner divisor has degree 0 and its height is a *single* quadratic quantity. Higher derivatives would need "cycles with $r$ Heegner legs" — these exist over function fields (shtukas, Yun–Zhang) and have no number-field analogue. |
| **Diagonal cycles** (Gross–Kudla–Schoen; Darmon–Rotger 2014–17) | cycles on triple products of modular curves | classes in $H^1(\mathbb{Q},V_pE\otimes\rho_1\otimes\rho_2)$ | $L(E,\rho_1\otimes\rho_2,1)\ne0$ | again rank-0 type for the *twisted* $L$-function; produces no rational points on $E$ itself beyond special cases. |
| **$p$-adic $L$-functions / Iwasawa theory** (Mazur–Swinnerton-Dyer, Kato, Skinner–Urban, …) | modular symbols, Eisenstein congruences | equality $\operatorname{char}_\Lambda(X_\infty) = (L_p(E,T))$ | $p$-adic order of vanishing ↔ Selmer corank (modulo Schneider) | $L_p$ is tied to the **complex** $L$-function only through *values* $L(E,\chi,1)$ at finite-order characters. The derivatives $L_p^{(k)}(0)$ interpolate nothing complex. Only in rank 1 is there a comparison, and it goes through the Heegner point ($p$-adic Gross–Zagier vs complex Gross–Zagier). |
| **Analytic non-vanishing** (Waldspurger, Bump–Friedberg–Hoffstein, Iwaniec–Sarnak, Ono–Skinner, …) | Rankin–Selberg, mollifiers, half-integral weight forms | many twists with $r_{\mathrm{an}}\le1$ | — | produces no points and no Selmer information; useful only in combination with the bridges above. |
| **Selmer statistics** (Bhargava–Shankar, Smith) | geometry of numbers, class-group counting | average Selmer sizes | — | proves things for a *proportion* of curves; blind to any individual curve of rank $\ge2$. |

The pattern: every bridge is a **rank-one Euler system**, and rank-one Euler systems see exactly the vanishing orders $0$ and $1$.

## 3.2 The rank-2 wall, stated as three separate missing theorems

Write $s_p = \operatorname{corank}_{\mathbb{Z}_p}\operatorname{Sel}_{p^\infty}(E/\mathbb{Q}) = r + \operatorname{corank}Ш[p^\infty]$. BSD for a curve of rank $r\ge2$ decomposes into:

**(W1) Upper bound on Selmer from $L$:** $\ L^{(r)}(E,1)\ne0 \implies s_p\le r$.
Known for $r=0$ (Kato, Kolyvagin) and $r=1$ (Kolyvagin). For $r\ge2$ this is what an Euler system **of rank $r$** would give
(Mazur–Rubin 2016 and Burns–Sano 2021 develop the formalism of "Kolyvagin systems of core rank $r$"), but no such system has been constructed for any
elliptic curve: its bottom class would be an element of $\bigwedge^r H^1(\mathbb{Q},T)$ whose "regulator" is $L^{(r)}(E,1)$, i.e. a *rank-$r$ zeta element* in the sense of
Perrin-Riou's conjectures (Astérisque 229). Constructing it is essentially equivalent to the problem.

**(W2) Lower bound on vanishing from points:** $\ r\ge k \implies L^{(j)}(E,1)=0$ for $j<k$.
Known for $k\le2$ (contrapositive of Kolyvagin and Gross–Zagier–Kolyvagin). For $k=3$ one would need: *three independent points force $L''(E,1)=0$.* No construction turns rational points into vanishing of $L$ beyond the first derivative. (Parity stretches $k\le2$ by one step: if the root number is $-1$, then $r\ge2$ already forces $\operatorname{ord}\ge3$. This is the Buhler–Gross–Zagier situation, and file 06 uses it for every rank-3 curve of conductor $<20000$. The first case parity cannot reach is root number $+1$ with $r\ge4$, where $\operatorname{ord}\ge4$ needs $L''(E,1)=0$.)

**(W3) Finiteness of $Ш$.** All known proofs of finiteness go through (W1): they bound $s_p$ from above by an $L$-value or $L$-derivative, and compare with explicit points. There is no "intrinsic" finiteness proof, not even for one curve of rank 2.

Everything proven in file 02 for ranks $\le1$ is (W1)+(W2)+(W3) in the only two cases where a rank-one Euler system suffices.

## 3.3 What Kolyvagin's derived classes *do* see (and what they do not)

Kolyvagin's construction attaches to each square-free $n$ (product of "Kolyvagin primes" $\ell$, inert in $K$ with $p^M \mid \ell+1$ and $p^M\mid a_\ell$)
a class $\kappa_n\in H^1(K,E[p^M])$ derived from the Heegner point of conductor $n$. Let $\nu(n)$ be the number of prime factors.

* **Kolyvagin's conjecture:** some $\kappa_n$ is non-zero (for $M$ large). **Theorem (W. Zhang 2014)** under hypotheses on $p$.
* **Kolyvagin's structure theorem (1991):** if the conjecture holds, let $M_\infty = \min\{\nu(n):\kappa_n\ne0\}$. Then
  $\operatorname{corank}\operatorname{Sel}_{p^\infty}(E/K) = M_\infty+1$, and the structure of $Ш(E/K)[p^\infty]$ (which is finite!) is read off from the orders of the $\kappa_n$.

So **Heegner-derived classes control Selmer groups in every rank** — the algebraic side of BSD is "known" in terms of Kolyvagin classes. What is completely
missing is the *analytic* input: a formula of the shape

$$L^{(M_\infty+1)}(E/K,1)\ \overset{?}{\sim}\ \text{"size" of the first non-vanishing } \kappa_n ,$$

a Gross–Zagier formula for derived Heegner classes. Nobody knows what the right-hand side should even be as a real number: the $\kappa_n$ are torsion classes
(elements of $H^1(K,E[p^M])$), not points with heights. This is, in my view, the sharpest way to say where the number-field theory ends: *the algebraic side has
a rank-$r$ theory (Kolyvagin classes), the analytic side has a rank-$\le1$ theory (Gross–Zagier), and no object lives on both sides.*

Over function fields the two sides **do** meet: Yun–Zhang's cycles on $\operatorname{Sht}^r$ (shtukas with $r$ legs) are geometric objects whose intersection numbers are $r$-th derivatives,
and whose cohomology classes play the role of the $\kappa_n$. The dream of "shtukas over number fields" (the motivation behind a great deal of recent arithmetic
geometry — Fargues–Scholze is the *local* version) is precisely the dream of a number-field object with $r$ legs.

## 3.4 Why Iwasawa theory cannot finish the job on its own

Suppose the main conjecture holds for $(E,p)$ and $p$ is good ordinary. Then (file 02, F3)

$$\operatorname{corank}\operatorname{Sel}_{p^\infty}(E/\mathbb{Q})\ \le\ \operatorname{ord}_{T=0}L_p(E,T),$$

with equality iff the $p$-adic height is non-degenerate. Two gaps remain:

1. **Schneider's conjecture** (non-degeneracy of the cyclotomic $p$-adic height pairing). Open in general; proved for CM curves (Bertrand). Without it, $L_p$ may vanish to higher order than the Selmer corank. It is not even known in general that $\hat h_p(P)\ne0$ for a single point of infinite order.
2. **No $p$-adic-to-complex comparison for derivatives.** $L_p(E,T)$ interpolates $L(E,\chi,1)$ for characters $\chi$ of $p$-power conductor; at $T=0$ it gives $(1-\alpha^{-1})^2L(E,1)/\Omega$. There is no reason, and no theorem, for $\operatorname{ord}_{T=0}L_p$ to equal $\operatorname{ord}_{s=1}L$ when both are $\ge2$. (For rank 1 the equality goes through the Heegner point: both derivatives are heights of $y_K$, one complex, one $p$-adic.)

So the best Iwasawa theory could deliver, even with Schneider, is the **$p$-adic** BSD rank statement $\operatorname{ord}_{T=0}L_p = r$ together with finiteness of $Ш[p^\infty]$ — for one $p$ at a time. That is already valuable (it is how `results/curated.md` certifies $Ш[p^\infty]$ finite for rank-2 and rank-3 curves), but it is a different conjecture from the complex one.

## 3.5 Why pure analysis cannot do it

Analytic number theory can prove non-vanishing (of $L$ or $L'$) for many twists in a family, and can compute $L^{(k)}(E,1)$ to any precision for a given curve. It cannot:
* prove that $L^{(k)}(E,1)=0$ for $k\ge1$ for a given curve (a real number equal to $0$ cannot be certified numerically; $L(E,1)=0$ is certified *algebraically*, by modular symbols);
* produce a rational point;
* say anything about $Ш$.

Any "proof" of BSD that uses only properties of $L(E,s)$ (functional equation, Euler product, growth, zeros) is therefore structurally incapable of success; it must at some point manufacture an element of $E(\mathbb{Q})$ or of a Galois cohomology group.

## 3.6 Why the function-field proof does not transfer

Over $\mathbb{F}_q(C)$, $E$ *is* a surface $\mathcal{E}\to C$ over a finite field, $L(E,s)$ is a factor of the zeta function of $\mathcal{E}$ (a polynomial, by Grothendieck–Deligne), rational points are sections, and $Ш$ is the Brauer group. BSD becomes a statement about algebraic cycles on a surface (the Tate conjecture), and finiteness of $Ш$ ⇔ BSD (Tate–Milne). The three tools that make this work — étale cohomology of varieties over finite fields, the Weil conjectures, and the geometric interpretation of $Ш$ — have no counterpart over $\mathbb{Q}$: there is no surface, no Frobenius, and $Ш$ is not a Brauer group of anything (the "arithmetic surface" $\mathcal{E}\to\operatorname{Spec}\mathbb{Z}$ has no good cohomological zeta-function theory at $s=1$).

## 3.7 Strategies currently being pursued (and how far they get)

1. **Higher-rank Euler systems** (Mazur–Rubin, Burns–Sano–Sakamoto, Kings–Loeffler–Zerbes "Euler systems for $\mathrm{GSp}_4$, $\mathrm{GU}(2,1)$…"). The *formalism* exists; constructions exist for other Galois representations where the relevant automorphic period is a value rather than a derivative. For $E/\mathbb{Q}$ itself nothing of rank $\ge2$ is constructed.
2. **Plectic cohomology** (Nekovář–Scholl 2016; Fornea–Gehrmann 2023, Fornea–Guitart–Masdeu). A conjectural theory in which products of Hilbert modular varieties carry "plectic" structure of weight equal to the number of real places, producing classes that should detect $r$-th derivatives for curves over totally real fields of degree $r$. Plectic Stark–Heegner points are constructed $p$-adically and tested numerically; the theory itself is conjectural.
3. **Darmon's program** (Stark–Heegner points, $p$-adic uniformisation, "elliptic Stark conjectures" of Bertolini–Darmon–Rotger). Conjectural $p$-adic constructions of points; proven cases again in rank $\le1$ situations.
4. **Converse theorems and statistics** (Skinner, Zhang, Bhargava–Skinner–Zhang, Smith). Push the proportion of curves with proven BSD toward $100\%$; cannot reach any individual rank-$\ge2$ curve.
5. **Computation** (Stein–Wuthrich, Jetchev–Lauter–Stein, Cremona, LMFDB). Tests refined conjectures (Kolyvagin's conjecture, Mazur–Tate refined conjectures, Kurihara's conjectures on Selmer structure via modular symbols) in situations where a structural insight might appear; proves $Ш[p^\infty]$ finite prime-by-prime for rank-$\ge2$ curves.
6. **Function-field → number-field dictionary** (Yun–Zhang; the search for arithmetic shtukas). The only setting with a higher-rank Gross–Zagier formula; the transfer problem is wide open.

None of these is close to a proof of (W1)–(W3) for a single curve of rank 2. A fair summary of 60 years of work: **the conjecture is proven exactly as far as rank-one Euler systems reach, and no one has built anything that reaches further.**
