# Rigorous enclosures of L^(k)(E,1)/k! and certified analytic ranks

*Generated 2026-11-07 04:08 by `research/rigorous-l/make_report.py` from `results.json`. Code: `bsd/rigorous_l.py`; tests: `tests/test_rigorous_l.py`; driver: `research/rigorous-l/run_batch.py`.*

**Scope.** Nobody has proved the Birch and Swinnerton-Dyer conjecture, and nothing here does. This component proves statements about the *analytic* side only. It gives Arb balls that provably contain c_k = L^(k)(E,1)/k! (k = 0..4), and from them, plus exact side information (parity, modular symbols), the order of vanishing ord_{s=1} L(E,s) for each curve. Every number labelled *rigorous* is an interval enclosure with explicit, proven truncation bounds. It is conditional only on the theorems and the software listed in section 7. Numbers labelled *numerical* (PARI's `lfun`, Cremona's table values) are floating-point references used for validation only.

## Summary

* **Validation**: 300 reference values (c_k and lambda_k = Lambda^(k)(1)/k!, k = 0..4, 25 curves including 11a1, 37a1, 389a1, 5077a1 and 571a1, at 128 and 256 bits). PARI's 135-digit value lies inside our ball in 300/300 cases. Largest radius: 2.8e-40 at 128 bits, 5.9e-79 at 256 bits.
* **Certification** of 997 curves (all curated curves with N < 20000, one curve for each isogeny class of Cremona rank 2 or 3 with N <= 10000, and 200 random rank-0/1 curves with N <= 10000): ord_{s=1} L(E,s) is **certified exactly for 995** curves. In every one of them it equals Cremona's rank (0 disagreements). For the 2 rank-3 curves the result is **ord ∈ {1, 3}** with c_3 ≠ 0 certified, plus a rigorous bound |L'(E,1)| ≤ ε. Exact vanishing of L'(E,1) cannot come from a ball; it needs e.g. a Gross–Zagier certificate. Errors: 0.
* **All 9 rank-3 curves with N < 20000**: c_3 ≠ 0 certified. |L'(E,1)| ≤ 1.6e-41 (128 bits) and ≤ 4.0e-157 (512 bits); see section 5.3.
* **Cost** (one core, kmax = 4, machine shared with other jobs): the enclosures for the 997 curves of the certification run took median 0.07 s and at most 0.30 s at 128 bits. Conductor 19990 takes 0.27 s at 128 bits and 4.1 s at 512 bits. Whole certificates, including the exact modular symbols needed for L(E,1) = 0 when w = +1, took median 0.50 s and at most 9.1 s.

## 1. Method

### 1.1 Functional equation and the integral at s = 1

By the modularity theorem, L(E,s) = L(f,s) for the weight-2 newform f = Σ a_n q^n of level N (the conductor). Put Λ(s) = N^{s/2} (2π)^{-s} Γ(s) L(E,s) and g(u) = f(iu/√N) = Σ_n a_n e^{-x_n u} with x_n = 2πn/√N. For Re s > 2 the series and integral converge absolutely, and termwise Mellin transformation gives

    Λ(s) = ∫_0^∞ g(u) u^{s-1} du .

f is an eigenform of the Fricke involution: f(-1/(Nz)) = ε N z² f(z) with ε = ±1. At z = iu/√N we have -1/(Nz) = i/(u√N) and N z² = -u², so g(1/u) = -ε u² g(u). Splitting the integral at u = 1 and substituting u -> 1/u on (0,1) gives

    Λ(s) = ∫_1^∞ g(u) (u^{s-1} - ε u^{1-s}) du ,

which is entire and satisfies Λ(2-s) = -ε Λ(s). Hence the root number is w = -ε, and **g(1/u) = w u² g(u)** as stated in the task. (Numerical check: g(1/u)/(u² g(u)) = +1.00000000000000 for 11a1 and 389a1 (w = +1) and -1.00000000000000 for 37a1 and 5077a1 (w = -1), at u = 2 and 3/2.) With s = 1 + t:

    Λ(1+t) = ∫_1^∞ g(u) (u^t + w u^{-t}) du ,     λ_k := Λ^(k)(1)/k! = (1 + w(-1)^k) Σ_n a_n K_k(x_n),
    K_k(x) := J_k(x)/k!,   J_k(x) = ∫_1^∞ e^{-xu} (log u)^k du .

The task's formulas are correct as given; no sign or normalisation had to be changed. Expanding u^{±t} = Σ (±log u)^k t^k/k! under the integral is justified by dominated convergence (|Σ_k| ≤ 2u^{|t|}, and g decays like e^{-2πu/√N}). Exchanging Σ_n and ∫ is Fubini, using Σ_n |a_n| J_k(x_n) < ∞ by (B1), (B2) below. For (-1)^k = -w the factor is exactly 0. These parity zeros are built into the enclosures as exact zeros; they are a theorem (the functional equation), not an assumption.

### 1.2 The special functions K_k(x), without cancellation

F(t) := ∫_1^∞ e^{-xu} u^t du = x^{-1-t} Γ(1+t, x) has Taylor coefficients [t^k] F = K_k(x). Using Γ(1+t,x) = Γ(1+t) - γ(1+t,x) and the series γ(a,x) = x^a e^{-x} Σ_{m≥0} x^m/(a)_{m+1} (DLMF 8.7.1 together with γ(a,x) = Γ(a) x^a γ*(a,x); this is Kummer's transformation of (x^a/a) 1F1(a; a+1; -x)), with a = 1 + t:

    x^{-1-t} Γ(1+t, x) = x^{-1-t} Γ(1+t)  -  e^{-x} Σ_{m≥0} x^m / (1+t)_{m+1} .

Here (-1)^k [t^k] 1/(1+t)_{m+1} = q_{k,m} := h_k(1, 1/2, …, 1/(m+1))/(m+1)! > 0, where h_k is the complete homogeneous symmetric polynomial; the series in m converges uniformly for |t| ≤ 1/2, so coefficients can be extracted termwise. Therefore

    K_k(x) = (1/x) Σ_{i=0}^k γ_i (-log x)^{k-i}/(k-i)!  -  (-1)^k e^{-x} Q_k(x),   Q_k(x) = Σ_{m≥0} q_{k,m} x^m,

with Γ(1+t) = Σ γ_i t^i (Arb's rigorous `arb_series.gamma`). Q_k has positive coefficients, so Horner evaluation at x > 0 does not cancel. The final subtraction loses only *relative* accuracy, and only where K_k(x) ≈ e^{-x}/x^{k+1} is itself negligible; the *absolute* error stays ~2^{-prec}. The textbook alternating series for γ(1+t,x) would lose ≈ 1.44 x bits (x reaches ≈ 96 at 128 bits). The q_{k,m} come from the exact recurrence q_{k,0} = 1, q_{k,m+1} = Σ_{i≤k} q_{k-i,m}/(m+2)^{i+1} in ball arithmetic. Checks: K_0 = e^{-x}/x and K_1 = E_1(x)/x (closed forms), Arb's rigorous `acb_series.hypgeom` for 1F1(1+t; 2+t; -x), and Arb's rigorous integrator `acb.integral` (scratch runs and `tests/test_rigorous_l.py`).

### 1.3 From Λ to L

L(E,1+t) = C(t) Λ(1+t), with C(t) = (2π/√N)^{1+t}/Γ(1+t) = (2π/√N) exp( t(log(2π/√N) + γ) - Σ_{m≥2} (-1)^m ζ(m) t^m/m ). Since log Γ(1+t) = -γt + Σ_{m≥2} (-1)^m ζ(m) t^m/m for |t| < 1, this is exact. The coefficients C_0..C_kmax depend only on the exponent up to t^kmax, so truncating the exponent there is exact (no remainder term). They are computed from Arb balls for π, γ, ζ(m) and log, and cross-checked against Arb's `rgamma` power series. Then c_k = Σ_{j≤k} λ_j C_{k-j}, with no assumption that lower derivatives vanish. Because C_0 = 2π/√N ≠ 0, the triangular system gives c_0 = … = c_{r-1} = 0 ⟺ λ_0 = … = λ_{r-1} = 0. Note that L itself is not even or odd about s = 1; for w = -1 only c_0 vanishes by parity. For example c_2 = λ_1 C_1 ≠ 0 for 37a1.

### 1.4 Implementation

`taylor_enclosures(ainvs, kmax, prec)` does the following:
1. Takes N (`ellglobalred`), w (`ellrootno`) and a_1..a_M (`ellan`, exact integers) from PARI.
2. Chooses M so that 2·(B3) ≤ 2^-(prec+4), and works at prec + 24 + log2 M + 4(kmax+1) bits.
3. For each n with a_n ≠ 0, evaluates K_k(x_n) for the k of the right parity. The degree D of Q_k is chosen per x_n, in steps of 8, so that the tail bound (B4) is negligible, and the rigorous (B4) is added to the ball.
4. Adds 2·(B3) to each λ_k and convolves with C(t).

Every choice (M, D, working precision) affects only the radius, never validity. `certify_order_rigorous` adds the exact information (section 3) and returns a JSON-ready dict. For every ball it gives Arb's enclosing decimal string, a midpoint, and rigorous decimal upper bounds for the radius and for |x|.

## 2. Proofs of the truncation bounds

**(B1) |a_n| ≤ d(n)√n ≤ 2n for all n ≥ 1.** a_n is multiplicative. For a good prime p, a_{p^{j+1}} = a_p a_{p^j} - p a_{p^{j-1}} with |a_p| ≤ 2√p (Hasse), so a_{p^j} = Σ_{i=0}^j α^i β^{j-i} with |α| = |β| = √p and |a_{p^j}| ≤ (j+1) p^{j/2} = d(p^j)(p^j)^{1/2}. For a bad prime, a_{p^j} = a_p^j with a_p ∈ {0, ±1}, so |a_{p^j}| ≤ 1. Multiplying gives |a_n| ≤ d(n)√n. Finally d(n) ≤ 2√n, because divisors pair as (d, n/d) and the smaller element of each pair is ≤ √n.

**(B2) 0 ≤ K_k(x) ≤ e^{-x}/x^{k+1} for x > 0.** On [1, ∞), 0 ≤ log u = ∫_1^u dv/v ≤ u - 1. Hence J_k(x) ≤ ∫_1^∞ e^{-xu}(u-1)^k du = e^{-x} ∫_0^∞ e^{-xv} v^k dv = e^{-x} k!/x^{k+1}.

**(B3) n-tail.** With a = 2π/√N and x_n = an:

    |Σ_{n>M} a_n K_k(x_n)| ≤ Σ_{n>M} 2n · e^{-an}/(an)^{k+1} = (2/a^{k+1}) Σ_{n>M} n^{-k} e^{-an}
                          ≤ 2 e^{-a(M+1)} / ( a^{k+1} (M+1)^k (1 - e^{-a}) ) .

This is evaluated in ball arithmetic; 2·(B3) is added as a radius to λ_k = 2 Σ a_n K_k(x_n).

**(B4) m-tail of Q_k.** Each monomial of h_k(y_1..y_r) occurs in (y_1+…+y_r)^k with coefficient ≥ 1, so for y_i ≥ 0 we have h_k ≤ (Σy_i)^k and q_{k,m} ≤ H_{m+1}^k/(m+1)! =: b_m/x^m (H = harmonic number). Since H_{m+2} = H_{m+1} + 1/(m+2) and H_{m+1} ≥ 1, b_{m+1}/b_m ≤ (1 + 1/(m+2))^k x/(m+2) ≤ ρ := (1 + 1/(D+3))^k x/(D+3) for m ≥ D+1. If ρ < 1 (checked in ball arithmetic), then

    0 ≤ Σ_{m>D} q_{k,m} x^m ≤ H_{D+2}^k x^{D+1} / ( (D+2)! (1 - ρ) ) .

**Rounding.** All other operations are Arb ball operations (FLINT 3 / python-flint 0.9.0). Exact zeros (parity) are exact balls. Validity therefore rests only on (B1)–(B4) and the correctness of Arb. `tests/test_rigorous_l.py` also checks (B3) against 400 explicit tail terms for 37a1.

## 3. Certification logic

* Exact zeros of λ_k come from *parity* (λ_k = 0 when (-1)^k = -w), from *modular symbols* (when w = +1 and the ball for λ_0 contains 0, `bsd.analytic.exact_L1_ratio` computes L(E,1)/Ω⁺ exactly as a rational via PARI `msfromell`/`mseval`, an exact linear-algebra computation), or from optional *external* certificates passed by another component (`external_zeros`, `external_lower_bound`). External inputs are recorded and named as external in the conclusion, and they are rejected if a ball contradicts them.
* A coefficient is *certified non-zero* when its ball excludes 0.
* ord ≥ r_lower := the number of leading coefficients proven zero. ord ≤ r_upper := the first index whose c- or λ-ball excludes 0, or 0 if modular symbols give L(E,1) ≠ 0. Parity of ord is that of (1-w)/2. The candidates are the r in [r_lower, r_upper] of the right parity with λ_r not proven zero. **ord is certified when exactly one candidate remains.**
* A ball around 0 never proves vanishing. For every candidate below the largest, the output lists rigorous bounds |c_k| ≤ ε, |λ_k| ≤ ε' and |L^(k)(E,1)| ≤ k!ε (`undetermined`). For w = -1 and analytic rank 3 the field `L_prime_abs_upper` is the bound on |L'(E,1)| that a Gross–Zagier certificate must turn into exact vanishing, for instance by showing that the Heegner point is torsion (Buhler–Gross–Zagier 1985 for 5077a). Alternatively, Gross–Zagier–Kolyvagin (ord ≤ 1 ⇒ rank = ord) together with a *proof* that rank E(Q) ≥ 2 gives ord ≥ 2, hence ord ≥ 3 by parity, and our certified c_3 ≠ 0 gives ord = 3. Such an input can be passed through `external_lower_bound`. Neither route is part of this component's results below.

## 4. Validation against PARI (floating point, 450 bits ≈ 135 digits)

Reference: c_k^PARI = `lfun(E,1,k)`/k! and λ_k^PARI = `lfunlambda(E,1,k)`/(2·k!). PARI's Λ carries the gamma factor Γ_R(s)Γ_R(s+1) = 2(2π)^{-s}Γ(s), twice ours. 'inside' means PARI's value lies in our ball. Where our coefficient is an exact zero (parity), PARI's value must be < 1e-60 in absolute value; it is 0.E-135.

| curve | N | w | k | our ball for c_k (128 bits) | PARI c_k (first digits) | inside |
|---|---|---|---|---|---|---|
| 11a1 | 11 | +1 | 0 | `[0.253841860855910684337758923350909461044 +/- 3.27…` | 0.25384186085591068433775892 | yes |
| 11a1 | 11 | +1 | 1 | `[0.308708533963172285620043118807185401970 +/- 3.77…` | 0.30870853396317228562004311 | yes |
| 11a1 | 11 | +1 | 2 | `[0.0113280542179201936918730233038576078627 +/- 5.6…` | 0.01132805421792019369187302 | yes |
| 11a1 | 11 | +1 | 3 | `[-0.0367068776204353155007383441280255754575 +/- 7.…` | -0.03670687762043531550073834 | yes |
| 11a1 | 11 | +1 | 4 | `[0.0088561798530780549620969310031189119623 +/- 2.6…` | 0.00885617985307805496209693 | yes |
| 37a1 | 37 | -1 | 0 | `0` | 0 | exact 0 (parity); PARI 0.E-135 |
| 37a1 | 37 | -1 | 1 | `[0.305999773834052301820483683321676474453 +/- 5.40…` | 0.30599977383405230182048368 | yes |
| 37a1 | 37 | -1 | 2 | `[0.186547797268161964173817368779507591454 +/- 1.99…` | 0.18654779726816196417381736 | yes |
| 37a1 | 37 | -1 | 3 | `[-0.136791463097187666302582216428158570769 +/- 3.1…` | -0.13679146309718766630258221 | yes |
| 37a1 | 37 | -1 | 4 | `[0.0161066468496400535050977294565270785182 +/- 1.6…` | 0.01610664684964005350509772 | yes |
| 389a1 | 389 | +1 | 0 | `[+/- 4.66e-41]` | 4.2337403813834066e-138 | yes |
| 389a1 | 389 | +1 | 1 | `[+/- 2.64e-41]` | -2.399247720994272e-138 | yes |
| 389a1 | 389 | +1 | 2 | `[0.7593165002884267702301926078947220190781 +/- 3.3…` | 0.75931650028842677023019260 | yes |
| 389a1 | 389 | +1 | 3 | `[-0.4303023375833619992903517750600442361904 +/- 5.…` | -0.43030233758336199929035177 | yes |
| 389a1 | 389 | +1 | 4 | `[-0.1935093138299810209672215497507966216008 +/- 5.…` | -0.19350931382998102096722154 | yes |
| 5077a1 | 5077 | -1 | 0 | `0` | 0 | exact 0 (parity); PARI 0.E-135 |
| 5077a1 | 5077 | -1 | 1 | `[+/- 1.56e-41]` | -2.779251361973017e-138 | yes |
| 5077a1 | 5077 | -1 | 2 | `[+/- 2.88e-41]` | 1.6052162830136721e-138 | yes |
| 5077a1 | 5077 | -1 | 3 | `[1.7318499001193006897919750850601528449544 +/- 2.1…` | 1.73184990011930068979197508 | yes |
| 5077a1 | 5077 | -1 | 4 | `[-3.2059055884439030547878000598586718161245 +/- 4.…` | -3.20590558844390305478780005 | yes |
| 571a1 | 571 | +1 | 0 | `[1.7293650250874438809412704044171868580140 +/- 8.1…` | 1.72936502508744388094127040 | yes |
| 571a1 | 571 | +1 | 1 | `[-1.3118995480263278564047737532927042624693 +/- 5.…` | -1.31189954802632785640477375 | yes |
| 571a1 | 571 | +1 | 2 | `[0.8883617697453767133392084279047020943713 +/- 7.1…` | 0.88836176974537671333920842 | yes |
| 571a1 | 571 | +1 | 3 | `[0.2706747383079536477639837917878724312867 +/- 7.2…` | 0.27067473830795364776398379 | yes |
| 571a1 | 571 | +1 | 4 | `[-1.1562788650365833146749061664165432072403 +/- 3.…` | -1.15627886503658331467490616 | yes |

All runs (every k = 0..4, both c_k and λ_k):

| curve | N | w | Cremona rank | bits | M | time (s) | max radius | all inside |
|---|---|---|---|---|---|---|---|---|
| 11a1 | 11 | +1 | 0 | 128 | 48 | 0.015 | 2.8e-40 | yes |
| 37a1 | 37 | -1 | 1 | 128 | 85 | 0.013 | 1.8e-40 | yes |
| 389a1 | 389 | +1 | 2 | 128 | 299 | 0.101 | 4.6e-41 | yes |
| 571a1 | 571 | +1 | 0 | 128 | 363 | 0.060 | 4.7e-41 | yes |
| 1638n2 | 1638 | +1 | 0 | 128 | 622 | 0.048 | 3.5e-41 | yes |
| 2856f3 | 2856 | -1 | 1 | 128 | 788 | 0.026 | 3.1e-41 | yes |
| 3146a1 | 3146 | -1 | 1 | 128 | 828 | 0.045 | 3.1e-41 | yes |
| 3276d2 | 3276 | -1 | 1 | 128 | 845 | 0.015 | 3.1e-41 | yes |
| 5077a1 | 5077 | -1 | 3 | 128 | 1057 | 0.117 | 2.9e-41 | yes |
| 6237c1 | 6237 | +1 | 0 | 128 | 1231 | 0.086 | 2.8e-41 | yes |
| 6258c1 | 6258 | +1 | 0 | 128 | 1233 | 0.131 | 2.8e-41 | yes |
| 6864x2 | 6864 | +1 | 0 | 128 | 1293 | 0.053 | 2.6e-41 | yes |
| 7644a1 | 7644 | -1 | 1 | 128 | 1303 | 0.029 | 2.5e-41 | yes |
| 8676c1 | 8676 | +1 | 0 | 128 | 1457 | 0.045 | 2.5e-41 | yes |
| 8850m1 | 8850 | -1 | 1 | 128 | 1404 | 0.080 | 2.5e-41 | yes |
| 9576h2 | 9576 | +1 | 0 | 128 | 1532 | 0.053 | 2.5e-41 | yes |
| 9888j1 | 9888 | +1 | 2 | 128 | 1557 | 0.077 | 2.5e-41 | yes |
| 11400a3 | 11400 | -1 | 1 | 128 | 1597 | 0.041 | 2.4e-41 | yes |
| 11601b1 | 11601 | +1 | 0 | 128 | 1690 | 0.107 | 2.3e-41 | yes |
| 12492b2 | 12492 | -1 | 1 | 128 | 1674 | 0.043 | 2.3e-41 | yes |
| 12540a1 | 12540 | +1 | 0 | 128 | 1758 | 0.090 | 2.3e-41 | yes |
| 13266b1 | 13266 | -1 | 1 | 128 | 1726 | 0.075 | 2.3e-41 | yes |
| 15158c2 | 15158 | +1 | 0 | 128 | 1937 | 0.301 | 2.1e-41 | yes |
| 15950j2 | 15950 | +1 | 2 | 128 | 1988 | 0.127 | 2.1e-41 | yes |
| 17550bn1 | 17550 | -1 | 1 | 128 | 1991 | 0.077 | 2.1e-41 | yes |
| 11a1 | 11 | +1 | 0 | 256 | 95 | 0.038 | 5.9e-79 | yes |
| 37a1 | 37 | -1 | 1 | 256 | 171 | 0.045 | 2.4e-79 | yes |
| 389a1 | 389 | +1 | 2 | 256 | 577 | 0.176 | 1.6e-79 | yes |
| 571a1 | 571 | +1 | 0 | 256 | 701 | 0.076 | 1.2e-79 | yes |
| 5077a1 | 5077 | -1 | 3 | 256 | 2056 | 0.399 | 8.1e-80 | yes |

Total: 300 of 300 reference values inside the balls. The midpoints differ from PARI by about 1e-42 at 128 bits, about 100 times less than our radii.

### 4.1 Two independent rigorous algorithms

`lambda_by_integration` computes λ_k with none of the K_k machinery. It integrates ∫_1^U g_M(u)(log u)^k du, where g_M(u) = P(e^{-au}) for the integer polynomial P = Σ_{n≤M} a_n X^n, with Arb's certified Gauss–Legendre integrator `acb.integral` (Petras' error bounds, log taken with the analytic flag). Two proven tails are added: |∫_U^∞ g_M (log u)^k| ≤ Σ_{n≤M} |a_n| Γ(k+1, x_n U)/x_n^{k+1} (using log u ≤ u), and k!·(B3) for n > M. On the same 25 curves the two balls overlap for every k (all_overlap = True). PARI's floating-point agreement checks against an outside reference; this checks the formula and the code against an unrelated rigorous method.

| curve | N | w | overlap (k=0..4) | series s | integral s | max rad series | max rad integral |
|---|---|---|---|---|---|---|---|
| 11a1 | 11 | +1 | yes | 0.01 | 0.01 | 1.2e-40 | 1.2e-40 |
| 37a1 | 37 | -1 | yes | 0.01 | 0.03 | 1.8e-40 | 1.8e-40 |
| 389a1 | 389 | +1 | yes | 0.07 | 0.18 | 1.4e-40 | 1.4e-40 |
| 571a1 | 571 | +1 | yes | 0.03 | 0.14 | 1.8e-40 | 1.8e-40 |
| 1638n2 | 1638 | +1 | yes | 0.10 | 0.24 | 1.8e-40 | 1.8e-40 |
| 2856f3 | 2856 | -1 | yes | 0.04 | 0.22 | 1.7e-40 | 1.7e-40 |
| 3146a1 | 3146 | -1 | yes | 0.07 | 0.18 | 1.7e-40 | 1.7e-40 |
| 3276d2 | 3276 | -1 | yes | 0.02 | 0.25 | 1.8e-40 | 1.8e-40 |
| 5077a1 | 5077 | -1 | yes | 0.14 | 0.31 | 1.7e-40 | 1.7e-40 |
| 6237c1 | 6237 | +1 | yes | 0.16 | 0.59 | 1.8e-40 | 6.0e-40 |
| 6258c1 | 6258 | +1 | yes | 0.18 | 0.60 | 1.8e-40 | 6.3e-40 |
| 6864x2 | 6864 | +1 | yes | 0.06 | 0.38 | 1.7e-40 | 8.1e-40 |
| 7644a1 | 7644 | -1 | yes | 0.04 | 0.20 | 1.7e-40 | 1.7e-40 |
| 8676c1 | 8676 | +1 | yes | 0.09 | 0.69 | 1.8e-40 | 1.4e-39 |
| 8850m1 | 8850 | -1 | yes | 0.11 | 0.39 | 1.7e-40 | 1.7e-40 |
| 9576h2 | 9576 | +1 | yes | 0.11 | 0.49 | 1.8e-40 | 1.8e-39 |
| 9888j1 | 9888 | +1 | yes | 0.07 | 0.55 | 1.8e-40 | 2.0e-39 |
| 11400a3 | 11400 | -1 | yes | 0.05 | 0.33 | 1.8e-40 | 1.8e-40 |
| 11601b1 | 11601 | +1 | yes | 0.12 | 0.68 | 1.7e-40 | 3.0e-39 |
| 12492b2 | 12492 | -1 | yes | 0.06 | 0.37 | 1.8e-40 | 1.8e-40 |
| 12540a1 | 12540 | +1 | yes | 0.12 | 0.79 | 1.8e-40 | 3.5e-39 |
| 13266b1 | 13266 | -1 | yes | 0.17 | 0.39 | 1.8e-40 | 1.8e-40 |
| 15158c2 | 15158 | +1 | yes | 0.35 | 1.05 | 1.8e-40 | 6.0e-39 |
| 15950j2 | 15950 | +1 | yes | 0.25 | 0.99 | 1.8e-40 | 6.6e-39 |
| 17550bn1 | 17550 | -1 | yes | 0.15 | 0.73 | 1.8e-40 | 1.8e-40 |

## 5. Certification results (128 bits, kmax = 4)

### 5.1 Counts

| Cremona rank | curves | ord certified exactly | ord ∈ {1,3} (rank-3 gap) | other / not certified |
|---|---|---|---|---|
| 0 | 126 | 126 | 0 | 0 |
| 1 | 99 | 99 | 0 | 0 |
| 2 | 770 | 770 | 0 | 0 |
| 3 | 2 | 0 | 2 | 0 |

| set | curves | certified | ord ∈ {1,3} |
|---|---|---|---|
| curated | 31 | 29 | 2 |
| random_rank01 | 200 | 200 | 0 |
| rank2_class_rep | 770 | 770 | 0 |
| rank3_class_rep | 1 | 0 | 1 |

* How ord is certified: rank 0, c_0's ball excludes 0. Rank 1 (w = -1), c_0 = 0 by parity and c_1's ball excludes 0. Rank 2 (w = +1), c_0 = 0 by modular symbols, c_1 = 0 by modular symbols plus parity, and c_2's ball excludes 0. Modular symbols were computed for 770 curves (823 s in total; slowest: 5586c1 8.9 s, 5712c1 8.8 s, 5760r1 8.4 s).
* Certified orders that disagree with Cremona's rank: 0. The largest relative difference between our leading coefficient and Cremona's tabulated L^(r)(1)/r! (printed to 15 significant digits) is 3.7e-15. Cremona's values were computed with different software (eclib/gp), so this is an independent numerical check.

### 5.2 Selected certificates

| curve | N | w | conclusion | leading coefficient (rigorous, conditioned on exact zeros) |
|---|---|---|---|---|
| 11a1 | 11 | +1 | ord_(s=1) L(E,s) = 0 rigorously: c_0 != 0 (ball for c_0 excludes 0). | `[0.253841860855910684337758923350909461044 +/- 3.…` |
| 37a1 | 37 | -1 | ord_(s=1) L(E,s) = 1 rigorously: c_0 = 0 (parity); c_1 != 0 (ball for c_1 excludes 0). | `[0.305999773834052301820483683321676474453 +/- 5.…` |
| 389a1 | 389 | +1 | ord_(s=1) L(E,s) = 2 rigorously: c_0 = 0 (modular symbols); c_1 = 0 (modular symbols; parity); c_2 != 0 (ball for c_2 excludes 0). | `[0.7593165002884267702301926078947220190781 +/- 2…` |
| 433a1 | 433 | +1 | ord_(s=1) L(E,s) = 2 rigorously: c_0 = 0 (modular symbols); c_1 = 0 (modular symbols; parity); c_2 != 0 (ball for c_2 excludes 0). | `[0.9470207808658145334897264009807512532043 +/- 4…` |
| 571a1 | 571 | +1 | ord_(s=1) L(E,s) = 0 rigorously: c_0 != 0 (ball for c_0 excludes 0). | `[1.729365025087443880941270404417186858014 +/- 8.…` |
| 571b1 | 571 | +1 | ord_(s=1) L(E,s) = 2 rigorously: c_0 = 0 (modular symbols); c_1 = 0 (modular symbols; parity); c_2 != 0 (ball for c_2 excludes 0). | `[0.9031605169716644699311486252572393786392 +/- 2…` |
| 681b1 | 681 | +1 | ord_(s=1) L(E,s) = 0 rigorously: c_0 != 0 (ball for c_0 excludes 0). | `[1.844815206126820716928527729159800755354 +/- 3.…` |
| 3364c1 | 3364 | +1 | ord_(s=1) L(E,s) = 0 rigorously: c_0 != 0 (ball for c_0 excludes 0). | `[5.210595695541168848304965225586067651633 +/- 1.…` |
| 5077a1 | 5077 | -1 | rigorous: 1 <= ord <= 3 and ord = 1 mod 2, so ord in [1, 3]. \|c_1\| <= 1.5532861e-41 (not proven zero). ord = 3 iff L'(E,1) = 0; exact vanishing needs another argument (e.g. Gross-Zagier: the Heegner | `c_3 = [1.731849900119300689791975085060152844954 …` |
| 11197a1 | 11197 | -1 | rigorous: 1 <= ord <= 3 and ord = 1 mod 2, so ord in [1, 3]. \|c_1\| <= 1.0493051e-41 (not proven zero). ord = 3 iff L'(E,1) = 0; exact vanishing needs another argument (e.g. Gross-Zagier: the Heegner | `c_3 = [2.871585754190565016616006139530539041317 …` |

### 5.3 Rank-3 curves (all with N < 20000): the bound a Gross–Zagier certificate must close

For each curve, c_0 = 0 exactly (parity), c_2 = λ_1 C_1 vanishes iff c_1 does, and c_3 ≠ 0 is certified. So ord ∈ {1,3}, with ord = 3 ⟺ L'(E,1) = 0. Rigorous upper bounds for |L'(E,1)| = |c_1|:

| curve | N | \|L'(E,1)\| ≤ (128 bits) | (256 bits) | (512 bits) | c_3 = L'''(E,1)/6 (512 bits, truncated) | time 512 bits (s) |
|---|---|---|---|---|---|---|
| 5077a1 | 5077 | 1.5532861e-41 | 4.4075627e-80 | 3.9758563e-157 | `[1.731849900119300689791975085060152844…` | 1.15 |
| 11197a1 | 11197 | 1.0493051e-41 | 3.0906627e-80 | 2.6635826e-157 | `[2.871585754190565016616006139530539041…` | 1.89 |
| 11642a1 | 11642 | 1.0213654e-41 | 3.0386751e-80 | 2.7094016e-157 | `[3.165929870472220684418157806176209517…` | 2.11 |
| 12279a1 | 12279 | 1.0244588e-41 | 3.0127639e-80 | 2.5334602e-157 | `[2.919032476745083165894859203432197165…` | 1.89 |
| 13766a1 | 13766 | 9.5121043e-42 | 2.8595215e-80 | 2.4514587e-157 | `[3.620272457612264809363321817968910493…` | 2.04 |
| 16811a1 | 16811 | 8.5882371e-42 | 2.5602603e-80 | 2.1788158e-157 | `[2.705328413403246815776206619728878029…` | 2.06 |
| 18097b1 | 18097 | 8.5554539e-42 | 2.4548476e-80 | 2.0874669e-157 | `[4.106943857290598778081931101453522876…` | 2.34 |
| 18562c1 | 18562 | 8.2001086e-42 | 2.4747192e-80 | 2.1352307e-157 | `[4.097177493674534133245676081586504479…` | 2.48 |
| 18745a1 | 18745 | 8.2759885e-42 | 2.4253689e-80 | 2.0790950e-157 | `[3.246308889404806198667895622520851251…` | 2.32 |

A precision of p bits gives |L'(E,1)| ≤ ~2^-p. Any finite precision leaves a gap; closing it needs an algebraic argument, not more digits.

### 5.4 Beyond N = 20000: rigorous upper bounds for the record-rank curves

The method has no conductor limit; only M ∝ √N grows. For the curated record curves, c_r ≠ 0 is certified, so ord ≤ r rigorously. Modular symbols were not run at these levels. Known exact methods (parity, modular symbols, Gross–Zagier–Kolyvagin with proven points) give at most ord ≥ 2 for w = +1 and ord ≥ 3 for w = -1, and *no* curve is known for which ord ≥ 4 can be proven. So these analytic ranks remain unproven, although the floating-point evidence for them is overwhelming.

| curve | N | w | PARI numerical rank | M | time (s) | certified non-zero c_k | candidates for ord | bounds for the other candidates |
|---|---|---|---|---|---|---|---|---|
| 234446a1 | 234446 | +1 | 4 | 7827 | 1.1 | [4] | [0, 2, 4] | \|c_0\| ≤ 2.3782879e-42, \|c_2\| ≤ 1.4922025e-41 |
| 19047851a1 | 19047851 | -1 | 5 | 70399 | 8.8 | [5] | [1, 3, 5] | \|c_1\| ≤ 2.6425751e-43, \|c_3\| ≤ 4.4857823e-42 |

The rank-7 curve (N ≈ 3.8·10^11, about 1.1·10^7 terms) is feasible but was not run. It needs about 0.8 GB and about 40 min with the current non-streaming handling of the a_n.

## 6. Timing (one core, kmax = 4; wall-clock seconds, machine shared with other jobs)

| curve | N | 128 bits: M / s | 256 bits: M / s | 512 bits: M / s |
|---|---|---|---|---|
| 11a1 | 11 | 48 / 0.014 | 95 / 0.039 | 189 / 0.138 |
| 37a1 | 37 | 85 / 0.011 | 171 / 0.034 | 342 / 0.141 |
| 100a1 | 100 | 149 / 0.013 | 290 / 0.038 | 573 / 0.159 |
| 389a1 | 389 | 299 / 0.049 | 577 / 0.192 | 1134 / 0.594 |
| 1001a1 | 1001 | 461 / 0.027 | 905 / 0.088 | 1795 / 0.498 |
| 2001a1 | 2001 | 689 / 0.121 | 1321 / 0.298 | 2584 / 1.178 |
| 5077a1 | 5077 | 1057 / 0.071 | 2056 / 0.245 | 4060 / 1.200 |
| 7501a1 | 7501 | 1353 / 0.163 | 2576 / 0.508 | 5022 / 2.374 |
| 10001a1 | 10001 | 1494 / 0.070 | 2896 / 0.237 | 5709 / 1.218 |
| 12502a1 | 12502 | 1674 / 0.086 | 3242 / 0.306 | 6387 / 1.405 |
| 15002a1 | 15002 | 1838 / 0.130 | 3554 / 0.432 | 7000 / 2.112 |
| 17502a1 | 17502 | 1988 / 0.166 | 3842 / 0.488 | 7564 / 2.525 |
| 19990a1 | 19990 | 2230 / 0.267 | 4227 / 0.814 | 8220 / 4.095 |

M grows like √N·prec, and the cost of one term grows with the degree D (∝ prec) and the cost of a multiplication, so the total cost grows roughly like √N·prec^2, a little faster at high precision. Conductor 20000 at 128 bits takes about 0.2–0.3 s, far below the one-minute target. For w = +1 curves with L(E,1) = 0, the exact modular-symbol computation (seconds near N = 10^4, and more than 1 GiB of PARI stack at some highly composite levels; this run used 2 GiB) costs much more than the enclosure.

## 7. What the certificates rest on (stated explicitly)

* Modularity theorem (Wiles; Taylor-Wiles; Breuil-Conrad-Diamond-Taylor): L(E,s) = L(f,s) for a weight-2 newform f of level N, so Lambda(s) is entire and Lambda(s) = w Lambda(2-s) with w the global root number; equivalently g(1/u) = w u^2 g(u) for g(u) = f(iu/sqrt N).
* a_n from PARI ellan (exact integers); conductor N from PARI ellglobalred (Tate's algorithm, exact); root number w from PARI ellrootno (exact algorithm: product of local root numbers).
* Hasse's bound |a_p| <= 2 sqrt p (theorem), giving |a_n| <= d(n) sqrt n <= 2n for all n >= 1.
* Correctness of the Arb ball arithmetic (FLINT 3 via python-flint) and of the PARI routines used.
* For L(E,1) = 0 when w = +1: PARI's `msfromell`/`mseval` (exact arithmetic over Q). The plus modular symbol of E evaluated on {0, ∞} is a non-zero rational multiple of L(E,1)/Ω⁺, so it vanishes iff L(E,1) does.
* Floating-point PARI and Cremona values are used for validation only, never as inputs to a conclusion.
* Limitations: certified non-vanishing is limited to k ≤ kmax = 4 here (the code accepts any kmax). A rank-3 certificate needs exact vanishing of L'(E,1), which this component does not provide. Certified ord = 2 is a statement about L(E,s) only. The rank part of BSD for these curves also needs rank E(Q) = 2, which is a different computation (descent plus independent points).

## 8. Interface for other components (e.g. a Gross–Zagier certificate)

```python
from bsd.rigorous_l import certify_order_rigorous, taylor_enclosures
res = certify_order_rigorous(ainvs)                # JSON-ready dict; res['conclusion'] is a sentence
eps = res['L_prime_abs_upper']                     # w = -1, rank-3 gap: rigorous |L'(E,1)| <= eps (decimal string)
enc = taylor_enclosures(ainvs, kmax=3, prec=512)   # Arb balls: enc.c[1] = L'(E,1), enc.c[3] = L^(3)(E,1)/6
twist = taylor_enclosures(twist_ainvs, kmax=0)     # any curve, e.g. a quadratic twist: ball for L(E^D,1)
# once L'(E,1) = 0 has been proven exactly elsewhere:
res3 = certify_order_rigorous(ainvs, external_zeros={1: 'Gross-Zagier: Heegner point is torsion'})
# -> res3['order'] == 3 and res3['order_certified'] is True; the source is quoted in res3['conclusion']
```

## 9. Reproduce

```
python3 -m pytest -q tests/test_rigorous_l.py          # ~2 s
python3 -m bsd.rigorous_l 389a1 --kmax 4 --prec 128     # one certificate as JSON
python3 research/rigorous-l/run_batch.py all            # resumable; skips finished items
python3 research/rigorous-l/aggregate.py && python3 research/rigorous-l/make_report.py
```
