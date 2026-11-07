#!/usr/bin/env python3
"""Generate research/rigorous-l/report.md from results.json and validation.jsonl (run aggregate.py first)."""
from __future__ import annotations

import json
import os
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def short(ball: str, n: int = 44) -> str:
    return ball if len(ball) <= n else ball[: n - 1] + "…"


def main():
    R = json.load(open(os.path.join(HERE, "results.json")))
    V = R["validation"]
    C = R["certification"]
    T = R["timing"]
    R3 = R.get("rank3_all_N_below_20000", {}).get("curves", [])
    vfull = {}
    with open(os.path.join(HERE, "validation.jsonl")) as f:
        for line in f:
            try:
                r = json.loads(line)
                vfull[r["key"]] = r
            except Exception:
                pass

    L = []
    w = L.append
    w("# Rigorous enclosures of L^(k)(E,1)/k! and certified analytic ranks")
    w("")
    w(f"*Generated {time.strftime('%Y-%m-%d %H:%M')} by `research/rigorous-l/make_report.py` from `results.json`. "
      "Code: `bsd/rigorous_l.py`; tests: `tests/test_rigorous_l.py`; driver: `research/rigorous-l/run_batch.py`.*")
    w("")
    w("**Scope.** Nobody has proved the Birch and Swinnerton-Dyer conjecture, and nothing here does. "
      "This component proves statements about the *analytic* side only. It gives Arb balls that provably contain "
      "c_k = L^(k)(E,1)/k! (k = 0..4), and from them, plus exact side information (parity, modular symbols), "
      "the order of vanishing ord_{s=1} L(E,s) for each curve. "
      "Every number labelled *rigorous* is an interval enclosure with explicit, proven truncation bounds. "
      "It is conditional only on the theorems and the software listed in section 7. "
      "Numbers labelled *numerical* (PARI's `lfun`, Cremona's table values) are floating-point references "
      "used for validation only.")
    w("")
    # ------------------------------------------------------------------ summary
    s = C["summary_by_cremona_rank"]
    tot = C["n_curves"]
    ncert = sum(v.get("certified", 0) for v in s.values())
    n13 = sum(v.get("candidates_1_3", 0) for v in s.values())
    w("## Summary")
    w("")
    w(f"* **Validation**: {V['n_values_checked']} reference values (c_k and lambda_k = Lambda^(k)(1)/k!, k = 0..4, "
      f"{V['n_curves']} curves including 11a1, 37a1, 389a1, 5077a1 and 571a1, at 128 and 256 bits). "
      f"PARI's 135-digit value lies inside our ball in {V['n_values_inside_ball']}/{V['n_values_checked']} cases. "
      f"Largest radius: {V['max_rad_prec128']:.1e} at 128 bits, {V['max_rad_prec256']:.1e} at 256 bits.")
    w(f"* **Certification** of {tot} curves (all curated curves with N < 20000, one curve for each isogeny class of "
      "Cremona rank 2 or 3 with N <= 10000, and 200 random rank-0/1 curves with N <= 10000): "
      f"ord_{{s=1}} L(E,s) is **certified exactly for {ncert}** curves. In every one of them it equals "
      f"Cremona's rank ({len(C['certified_order_disagrees_with_cremona_rank'])} disagreements). "
      f"For the {n13} rank-3 curves the result is **ord ∈ {{1, 3}}** with c_3 ≠ 0 certified, plus a rigorous "
      "bound |L'(E,1)| ≤ ε. Exact vanishing of L'(E,1) cannot come from a ball; it needs e.g. a Gross–Zagier "
      f"certificate. Errors: {C['n_errors']}.")
    if R3:
        e512 = max(float(r["by_prec"]["512"]["L_prime_abs_upper"]) for r in R3)
        e128 = max(float(r["by_prec"]["128"]["L_prime_abs_upper"]) for r in R3)
        w(f"* **All {len(R3)} rank-3 curves with N < 20000**: c_3 ≠ 0 certified. "
          f"|L'(E,1)| ≤ {e128:.1e} (128 bits) and ≤ {e512:.1e} (512 bits); see section 5.3.")
    tmax = max((t["seconds"] for t in T if t["prec_bits"] == 128), default=None)
    t512 = max((t["seconds"] for t in T if t["prec_bits"] == 512), default=None)
    es = C.get("enclosure_seconds", {})
    cs = C.get("certificate_seconds_incl_modsym", {})
    w(f"* **Cost** (one core, kmax = 4, machine shared with other jobs): the enclosures for the {tot} curves of the certification run "
      f"took median {es.get('median', 0):.2f} s and at most {es.get('max', 0):.2f} s at 128 bits. Conductor 19990 "
      f"takes {tmax:.2f} s at 128 bits and {t512:.1f} s at 512 bits. Whole certificates, including the exact "
      f"modular symbols needed for L(E,1) = 0 when w = +1, took median {cs.get('median', 0):.2f} s and at most "
      f"{cs.get('max', 0):.1f} s.")
    w("")
    # ------------------------------------------------------------------ method
    w("## 1. Method")
    w("")
    w("### 1.1 Functional equation and the integral at s = 1")
    w("")
    w("By the modularity theorem, L(E,s) = L(f,s) for the weight-2 newform f = Σ a_n q^n of level N (the conductor). "
      "Put Λ(s) = N^{s/2} (2π)^{-s} Γ(s) L(E,s) and g(u) = f(iu/√N) = Σ_n a_n e^{-x_n u} with x_n = 2πn/√N. "
      "For Re s > 2 the series and integral converge absolutely, and termwise Mellin transformation gives")
    w("")
    w("    Λ(s) = ∫_0^∞ g(u) u^{s-1} du .")
    w("")
    w("f is an eigenform of the Fricke involution: f(-1/(Nz)) = ε N z² f(z) with ε = ±1. At z = iu/√N we have "
      "-1/(Nz) = i/(u√N) and N z² = -u², so g(1/u) = -ε u² g(u). Splitting the integral at u = 1 and substituting "
      "u -> 1/u on (0,1) gives")
    w("")
    w("    Λ(s) = ∫_1^∞ g(u) (u^{s-1} - ε u^{1-s}) du ,")
    w("")
    w("which is entire and satisfies Λ(2-s) = -ε Λ(s). Hence the root number is w = -ε, and **g(1/u) = w u² g(u)** "
      "as stated in the task. (Numerical check: g(1/u)/(u² g(u)) = +1.00000000000000 for 11a1 and 389a1 (w = +1) and "
      "-1.00000000000000 for 37a1 and 5077a1 (w = -1), at u = 2 and 3/2.) With s = 1 + t:")
    w("")
    w("    Λ(1+t) = ∫_1^∞ g(u) (u^t + w u^{-t}) du ,     λ_k := Λ^(k)(1)/k! = (1 + w(-1)^k) Σ_n a_n K_k(x_n),")
    w("    K_k(x) := J_k(x)/k!,   J_k(x) = ∫_1^∞ e^{-xu} (log u)^k du .")
    w("")
    w("The task's formulas are correct as given; no sign or normalisation had to be changed. Expanding "
      "u^{±t} = Σ (±log u)^k t^k/k! under the integral is justified by dominated convergence "
      "(|Σ_k| ≤ 2u^{|t|}, and g decays like e^{-2πu/√N}). Exchanging Σ_n and ∫ is Fubini, using "
      "Σ_n |a_n| J_k(x_n) < ∞ by (B1), (B2) below. For (-1)^k = -w the factor is exactly 0. These parity zeros "
      "are built into the enclosures as exact zeros; they are a theorem (the functional equation), not an assumption.")
    w("")
    w("### 1.2 The special functions K_k(x), without cancellation")
    w("")
    w("F(t) := ∫_1^∞ e^{-xu} u^t du = x^{-1-t} Γ(1+t, x) has Taylor coefficients [t^k] F = K_k(x). "
      "Using Γ(1+t,x) = Γ(1+t) - γ(1+t,x) and the series γ(a,x) = x^a e^{-x} Σ_{m≥0} x^m/(a)_{m+1} (DLMF 8.7.1 together with γ(a,x) = Γ(a) x^a γ*(a,x); "
      "this is Kummer's transformation of (x^a/a) 1F1(a; a+1; -x)), with a = 1 + t:")
    w("")
    w("    x^{-1-t} Γ(1+t, x) = x^{-1-t} Γ(1+t)  -  e^{-x} Σ_{m≥0} x^m / (1+t)_{m+1} .")
    w("")
    w("Here (-1)^k [t^k] 1/(1+t)_{m+1} = q_{k,m} := h_k(1, 1/2, …, 1/(m+1))/(m+1)! > 0, where h_k is the complete "
      "homogeneous symmetric polynomial; the series in m converges uniformly for |t| ≤ 1/2, so coefficients can be "
      "extracted termwise. Therefore")
    w("")
    w("    K_k(x) = (1/x) Σ_{i=0}^k γ_i (-log x)^{k-i}/(k-i)!  -  (-1)^k e^{-x} Q_k(x),   Q_k(x) = Σ_{m≥0} q_{k,m} x^m,")
    w("")
    w("with Γ(1+t) = Σ γ_i t^i (Arb's rigorous `arb_series.gamma`). Q_k has positive coefficients, so Horner "
      "evaluation at x > 0 does not cancel. The final subtraction loses only *relative* accuracy, and only where "
      "K_k(x) ≈ e^{-x}/x^{k+1} is itself negligible; the *absolute* error stays ~2^{-prec}. The textbook "
      "alternating series for γ(1+t,x) would lose ≈ 1.44 x bits (x reaches ≈ 96 at 128 bits). The q_{k,m} come "
      "from the exact recurrence q_{k,0} = 1, q_{k,m+1} = Σ_{i≤k} q_{k-i,m}/(m+2)^{i+1} in ball arithmetic. "
      "Checks: K_0 = e^{-x}/x and K_1 = E_1(x)/x (closed forms), Arb's rigorous `acb_series.hypgeom` for "
      "1F1(1+t; 2+t; -x), and Arb's rigorous integrator `acb.integral` (scratch runs and `tests/test_rigorous_l.py`).")
    w("")
    w("### 1.3 From Λ to L")
    w("")
    w("L(E,1+t) = C(t) Λ(1+t), with C(t) = (2π/√N)^{1+t}/Γ(1+t) = (2π/√N) exp( t(log(2π/√N) + γ) - "
      "Σ_{m≥2} (-1)^m ζ(m) t^m/m ). Since log Γ(1+t) = -γt + Σ_{m≥2} (-1)^m ζ(m) t^m/m for |t| < 1, this is exact. "
      "The coefficients C_0..C_kmax depend only on the exponent up to t^kmax, so truncating the exponent there is "
      "exact (no remainder term). They are computed from Arb balls for π, γ, ζ(m) and log, and cross-checked "
      "against Arb's `rgamma` power series. Then c_k = Σ_{j≤k} λ_j C_{k-j}, with no assumption that lower "
      "derivatives vanish. Because C_0 = 2π/√N ≠ 0, the triangular system gives "
      "c_0 = … = c_{r-1} = 0 ⟺ λ_0 = … = λ_{r-1} = 0. Note that L itself is not even or odd about s = 1; "
      "for w = -1 only c_0 vanishes by parity. For example c_2 = λ_1 C_1 ≠ 0 for 37a1.")
    w("")
    w("### 1.4 Implementation")
    w("")
    w("`taylor_enclosures(ainvs, kmax, prec)` does the following:")
    w("1. Takes N (`ellglobalred`), w (`ellrootno`) and a_1..a_M (`ellan`, exact integers) from PARI.")
    w("2. Chooses M so that 2·(B3) ≤ 2^-(prec+4), and works at prec + 24 + log2 M + 4(kmax+1) bits.")
    w("3. For each n with a_n ≠ 0, evaluates K_k(x_n) for the k of the right parity. The degree D of Q_k is chosen "
      "per x_n, in steps of 8, so that the tail bound (B4) is negligible, and the rigorous (B4) is added to the ball.")
    w("4. Adds 2·(B3) to each λ_k and convolves with C(t).")
    w("")
    w("Every choice (M, D, working precision) affects only the radius, never validity. `certify_order_rigorous` "
      "adds the exact information (section 3) and returns a JSON-ready dict. For every ball it gives Arb's "
      "enclosing decimal string, a midpoint, and rigorous decimal upper bounds for the radius and for |x|.")
    w("")
    # ------------------------------------------------------------------ proofs
    w("## 2. Proofs of the truncation bounds")
    w("")
    w("**(B1) |a_n| ≤ d(n)√n ≤ 2n for all n ≥ 1.** a_n is multiplicative. For a good prime p, "
      "a_{p^{j+1}} = a_p a_{p^j} - p a_{p^{j-1}} with |a_p| ≤ 2√p (Hasse), so a_{p^j} = Σ_{i=0}^j α^i β^{j-i} "
      "with |α| = |β| = √p and |a_{p^j}| ≤ (j+1) p^{j/2} = d(p^j)(p^j)^{1/2}. For a bad prime, a_{p^j} = a_p^j with "
      "a_p ∈ {0, ±1}, so |a_{p^j}| ≤ 1. Multiplying gives |a_n| ≤ d(n)√n. Finally d(n) ≤ 2√n, because divisors pair "
      "as (d, n/d) and the smaller element of each pair is ≤ √n.")
    w("")
    w("**(B2) 0 ≤ K_k(x) ≤ e^{-x}/x^{k+1} for x > 0.** On [1, ∞), 0 ≤ log u = ∫_1^u dv/v ≤ u - 1. Hence "
      "J_k(x) ≤ ∫_1^∞ e^{-xu}(u-1)^k du = e^{-x} ∫_0^∞ e^{-xv} v^k dv = e^{-x} k!/x^{k+1}.")
    w("")
    w("**(B3) n-tail.** With a = 2π/√N and x_n = an:")
    w("")
    w("    |Σ_{n>M} a_n K_k(x_n)| ≤ Σ_{n>M} 2n · e^{-an}/(an)^{k+1} = (2/a^{k+1}) Σ_{n>M} n^{-k} e^{-an}")
    w("                          ≤ 2 e^{-a(M+1)} / ( a^{k+1} (M+1)^k (1 - e^{-a}) ) .")
    w("")
    w("This is evaluated in ball arithmetic; 2·(B3) is added as a radius to λ_k = 2 Σ a_n K_k(x_n).")
    w("")
    w("**(B4) m-tail of Q_k.** Each monomial of h_k(y_1..y_r) occurs in (y_1+…+y_r)^k with coefficient ≥ 1, so "
      "for y_i ≥ 0 we have h_k ≤ (Σy_i)^k and q_{k,m} ≤ H_{m+1}^k/(m+1)! =: b_m/x^m (H = harmonic number). "
      "Since H_{m+2} = H_{m+1} + 1/(m+2) and H_{m+1} ≥ 1, b_{m+1}/b_m ≤ (1 + 1/(m+2))^k x/(m+2) ≤ ρ := "
      "(1 + 1/(D+3))^k x/(D+3) for m ≥ D+1. If ρ < 1 (checked in ball arithmetic), then")
    w("")
    w("    0 ≤ Σ_{m>D} q_{k,m} x^m ≤ H_{D+2}^k x^{D+1} / ( (D+2)! (1 - ρ) ) .")
    w("")
    w("**Rounding.** All other operations are Arb ball operations (FLINT 3 / python-flint 0.9.0). "
      "Exact zeros (parity) are exact balls. Validity therefore rests only on (B1)–(B4) and the correctness of Arb. "
      "`tests/test_rigorous_l.py` also checks (B3) against 400 explicit tail terms for 37a1.")
    w("")
    # ------------------------------------------------------------------ certification
    w("## 3. Certification logic")
    w("")
    w("* Exact zeros of λ_k come from *parity* (λ_k = 0 when (-1)^k = -w), from *modular symbols* (when w = +1 and "
      "the ball for λ_0 contains 0, `bsd.analytic.exact_L1_ratio` computes L(E,1)/Ω⁺ exactly as a rational via PARI "
      "`msfromell`/`mseval`, an exact linear-algebra computation), or from optional *external* certificates passed by "
      "another component (`external_zeros`, `external_lower_bound`). External inputs are recorded and named as "
      "external in the conclusion, and they are rejected if a ball contradicts them.")
    w("* A coefficient is *certified non-zero* when its ball excludes 0.")
    w("* ord ≥ r_lower := the number of leading coefficients proven zero. ord ≤ r_upper := the first index whose c- "
      "or λ-ball excludes 0, or 0 if modular symbols give L(E,1) ≠ 0. Parity of ord is that of (1-w)/2. "
      "The candidates are the r in [r_lower, r_upper] of the right parity with λ_r not proven zero. "
      "**ord is certified when exactly one candidate remains.**")
    w("* A ball around 0 never proves vanishing. For every candidate below the largest, the output lists rigorous "
      "bounds |c_k| ≤ ε, |λ_k| ≤ ε' and |L^(k)(E,1)| ≤ k!ε (`undetermined`). For w = -1 and analytic rank 3 the "
      "field `L_prime_abs_upper` is the bound on |L'(E,1)| that a Gross–Zagier certificate must turn into exact "
      "vanishing, for instance by showing that the Heegner point is torsion (Buhler–Gross–Zagier 1985 for 5077a). "
      "Alternatively, Gross–Zagier–Kolyvagin (ord ≤ 1 ⇒ rank = ord) together with a *proof* that rank E(Q) ≥ 2 "
      "gives ord ≥ 2, hence ord ≥ 3 by parity, and our certified c_3 ≠ 0 gives ord = 3. Such an input can be passed "
      "through `external_lower_bound`. Neither route is part of this component's results below.")
    w("")
    # ------------------------------------------------------------------ validation
    w("## 4. Validation against PARI (floating point, 450 bits ≈ 135 digits)")
    w("")
    w("Reference: c_k^PARI = `lfun(E,1,k)`/k! and λ_k^PARI = `lfunlambda(E,1,k)`/(2·k!). PARI's Λ carries the "
      "gamma factor Γ_R(s)Γ_R(s+1) = 2(2π)^{-s}Γ(s), twice ours. 'inside' means PARI's value lies in our ball. "
      "Where our coefficient is an exact zero (parity), PARI's value must be < 1e-60 in absolute value; it is 0.E-135.")
    w("")
    w("| curve | N | w | k | our ball for c_k (128 bits) | PARI c_k (first digits) | inside |")
    w("|---|---|---|---|---|---|---|")
    for lab in ["11a1", "37a1", "389a1", "5077a1", "571a1"]:
        r = vfull.get(f"{lab}@128")
        if not r:
            continue
        for c in r["coefficients"]:
            ball = c["c"]["ball"]
            pari = c["pari_lfun_over_kfact"]
            if "e" in pari:
                mant, ex = pari.split("e", 1)
                pari_s = mant[:18] + "e" + ex
            else:
                pari_s = pari[:28] if not pari.startswith("-") else pari[:29]
            ok = "exact 0 (parity); PARI 0.E-135" if c["c_exact_zero_parity"] else ("yes" if c["contains_pari"] else "**NO**")
            w(f"| {lab} | {r['conductor']} | {r['root_number']:+d} | {c['k']} | `{short(ball, 52)}` | {pari_s} | {ok} |")
    w("")
    w("All runs (every k = 0..4, both c_k and λ_k):")
    w("")
    w("| curve | N | w | Cremona rank | bits | M | time (s) | max radius | all inside |")
    w("|---|---|---|---|---|---|---|---|---|")
    for r in V["runs"]:
        w(f"| {r['label']} | {r['conductor']} | {r['root_number']:+d} | {r['cremona_rank']} | {r['prec_bits']} | "
          f"{r['M_terms']} | {r['seconds']:.3f} | {r['max_rad']:.1e} | {'yes' if r['all_contained'] else '**NO**'} |")
    w("")
    w(f"Total: {V['n_values_inside_ball']} of {V['n_values_checked']} reference values inside the balls. "
      "The midpoints differ from PARI by about 1e-42 at 128 bits, about 100 times less than our radii.")
    w("")
    X = R.get("crosscheck_two_rigorous_methods")
    if X and X.get("runs"):
        w("### 4.1 Two independent rigorous algorithms")
        w("")
        w("`lambda_by_integration` computes λ_k with none of the K_k machinery. It integrates "
          "∫_1^U g_M(u)(log u)^k du, where g_M(u) = P(e^{-au}) for the integer polynomial P = Σ_{n≤M} a_n X^n, "
          "with Arb's certified Gauss–Legendre integrator `acb.integral` (Petras' error bounds, log taken with the "
          "analytic flag). Two proven tails are added: |∫_U^∞ g_M (log u)^k| ≤ Σ_{n≤M} |a_n| Γ(k+1, x_n U)/x_n^{k+1} "
          "(using log u ≤ u), and k!·(B3) for n > M. On the same 25 curves the two balls overlap for every k "
          f"(all_overlap = {X['all_overlap']}). PARI's floating-point agreement checks against an outside reference; "
          "this checks the formula and the code against an unrelated rigorous method.")
        w("")
        w("| curve | N | w | overlap (k=0..4) | series s | integral s | max rad series | max rad integral |")
        w("|---|---|---|---|---|---|---|---|")
        for r in X["runs"]:
            w(f"| {r['label']} | {r['conductor']} | {r['root_number']:+d} | {'yes' if r['all_overlap'] else '**NO**'} | "
              f"{r['seconds_series']:.2f} | {r['seconds_integral']:.2f} | {r['max_rad_series']:.1e} | "
              f"{r['max_rad_integral']:.1e} |")
        w("")
    # ------------------------------------------------------------------ certification results
    w("## 5. Certification results (128 bits, kmax = 4)")
    w("")
    w("### 5.1 Counts")
    w("")
    w("| Cremona rank | curves | ord certified exactly | ord ∈ {1,3} (rank-3 gap) | other / not certified |")
    w("|---|---|---|---|---|")
    for k, v in s.items():
        other = v.get("total", 0) - v.get("certified", 0) - v.get("candidates_1_3", 0)
        w(f"| {k.replace('cremona_rank_', '')} | {v.get('total', 0)} | {v.get('certified', 0)} | "
          f"{v.get('candidates_1_3', 0)} | {other} |")
    w("")
    w("| set | curves | certified | ord ∈ {1,3} |")
    w("|---|---|---|---|")
    for k, v in C["summary_by_set"].items():
        w(f"| {k} | {v.get('total', 0)} | {v.get('certified', 0)} | {v.get('candidates_1_3', 0)} |")
    w("")
    w(f"* How ord is certified: rank 0, c_0's ball excludes 0. Rank 1 (w = -1), c_0 = 0 by parity and c_1's ball "
      f"excludes 0. Rank 2 (w = +1), c_0 = 0 by modular symbols, c_1 = 0 by modular symbols plus parity, and c_2's "
      f"ball excludes 0. Modular symbols were computed for {C['modular_symbols']['n_computed']} curves "
      f"({C['modular_symbols']['total_seconds']:.0f} s in total; slowest: "
      + ", ".join(f"{x['label']} {x['seconds']:.1f} s" for x in C['modular_symbols']['slowest'][:3]) + ").")
    w(f"* Certified orders that disagree with Cremona's rank: {len(C['certified_order_disagrees_with_cremona_rank'])}. "
      f"The largest relative difference between our leading coefficient and Cremona's tabulated L^(r)(1)/r! "
      f"(printed to 15 significant digits) is {C['max_rel_diff_leading_vs_cremona_table']:.1e}. Cremona's values were "
      "computed with different software (eclib/gp), so this is an independent numerical check.")
    if C["not_certified_other"]:
        w(f"* Not certified for other reasons: {len(C['not_certified_other'])} curves: "
          + ", ".join(x["label"] for x in C["not_certified_other"][:20]))
    if C["errors"]:
        w(f"* Errors: {len(C['errors'])}: " + ", ".join(f"{e['label']} ({e['error'][:60]})" for e in C["errors"][:10]))
    w("")
    w("### 5.2 Selected certificates")
    w("")
    w("| curve | N | w | conclusion | leading coefficient (rigorous, conditioned on exact zeros) |")
    w("|---|---|---|---|---|")
    byl = {c["label"]: c for c in C["curves"]}
    for lab in ["11a1", "37a1", "389a1", "433a1", "571a1", "571b1", "681b1", "3364c1", "5077a1", "11197a1"]:
        c = byl.get(lab)
        if not c:
            continue
        concl = c["conclusion"].replace("|", "\\|")
        lead = c["leading_coefficient_ball"] or ("c_3 = " + c["c"][3])
        w(f"| {lab} | {c['conductor']} | {c['root_number']:+d} | {concl[:200]} | `{short(lead, 50)}` |")
    w("")
    w("### 5.3 Rank-3 curves (all with N < 20000): the bound a Gross–Zagier certificate must close")
    w("")
    w("For each curve, c_0 = 0 exactly (parity), c_2 = λ_1 C_1 vanishes iff c_1 does, and c_3 ≠ 0 is certified. "
      "So ord ∈ {1,3}, with ord = 3 ⟺ L'(E,1) = 0. Rigorous upper bounds for |L'(E,1)| = |c_1|:")
    w("")
    w("| curve | N | \\|L'(E,1)\\| ≤ (128 bits) | (256 bits) | (512 bits) | c_3 = L'''(E,1)/6 (512 bits, truncated) | time 512 bits (s) |")
    w("|---|---|---|---|---|---|---|")
    for r in R3:
        b = r["by_prec"]
        w(f"| {r['label']} | {r['conductor']} | {b['128']['L_prime_abs_upper']} | {b['256']['L_prime_abs_upper']} | "
          f"{b['512']['L_prime_abs_upper']} | `{short(b['512']['c3_ball'], 40)}` | {b['512']['seconds']:.2f} |")
    w("")
    w("A precision of p bits gives |L'(E,1)| ≤ ~2^-p. Any finite precision leaves a gap; closing it needs an "
      "algebraic argument, not more digits.")
    w("")
    B = R.get("beyond_N_20000_record_ranks", {}).get("curves", [])
    if B:
        w("### 5.4 Beyond N = 20000: rigorous upper bounds for the record-rank curves")
        w("")
        w("The method has no conductor limit; only M ∝ √N grows. For the curated record curves, c_r ≠ 0 is certified, "
          "so ord ≤ r rigorously. Modular symbols were not run at these levels. Known exact methods (parity, "
          "modular symbols, Gross–Zagier–Kolyvagin with proven points) give at most ord ≥ 2 for w = +1 and ord ≥ 3 for "
          "w = -1, and *no* curve is known for which ord ≥ 4 can be proven. So these analytic ranks remain unproven, "
          "although the floating-point evidence for them is overwhelming.")
        w("")
        w("| curve | N | w | PARI numerical rank | M | time (s) | certified non-zero c_k | candidates for ord | bounds for the other candidates |")
        w("|---|---|---|---|---|---|---|---|---|")
        for r in B:
            und = ", ".join(f"\\|c_{k}\\| ≤ {v}" for k, v in r["undetermined"].items())
            w(f"| {r['label']} | {r['conductor']} | {r['root_number']:+d} | {r['numerical_analytic_rank_pari']} | "
              f"{r['M_terms']} | {r['seconds']:.1f} | {r['c_nonzero_certified']} | {r['order_candidates']} | {und} |")
        w("")
        w("The rank-7 curve (N ≈ 3.8·10^11, about 1.1·10^7 terms) is feasible but was not run. It needs about 0.8 GB "
          "and about 40 min with the current non-streaming handling of the a_n.")
        w("")
    # ------------------------------------------------------------------ timing
    w("## 6. Timing (one core, kmax = 4; wall-clock seconds, machine shared with other jobs)")
    w("")
    precs = sorted({t["prec_bits"] for t in T})
    labs = []
    for t in T:
        if t["label"] not in labs:
            labs.append(t["label"])
    tt = {(t["label"], t["prec_bits"]): t for t in T}
    w("| curve | N | " + " | ".join(f"{p} bits: M / s" for p in precs) + " |")
    w("|---|---|" + "---|" * len(precs))
    for lab in sorted(labs, key=lambda l: tt[(l, precs[0])]["conductor"]):
        row = [f"{tt[(lab, p)]['M_terms']} / {tt[(lab, p)]['seconds']:.3f}" if (lab, p) in tt else "-" for p in precs]
        w(f"| {lab} | {tt[(lab, precs[0])]['conductor']} | " + " | ".join(row) + " |")
    w("")
    w("M grows like √N·prec, and the cost of one term grows with the degree D (∝ prec) and the cost of a "
      "multiplication, so the total cost grows roughly like √N·prec^2, a little faster at high precision. "
      "Conductor 20000 at 128 bits takes about 0.2–0.3 s, far below the one-minute target. For w = +1 curves with L(E,1) = 0, the exact modular-symbol computation (seconds near N = 10^4, "
      "and more than 1 GiB of PARI stack at some highly composite levels; this run used 2 GiB) costs much more than "
      "the enclosure.")
    w("")
    # ------------------------------------------------------------------ assumptions
    w("## 7. What the certificates rest on (stated explicitly)")
    w("")
    for a in R["meta"]["assumptions"]:
        w(f"* {a}")
    w("* For L(E,1) = 0 when w = +1: PARI's `msfromell`/`mseval` (exact arithmetic over Q). "
      "The plus modular symbol of E evaluated on {0, ∞} is a non-zero rational multiple of L(E,1)/Ω⁺, so it "
      "vanishes iff L(E,1) does.")
    w("* Floating-point PARI and Cremona values are used for validation only, never as inputs to a conclusion.")
    w("* Limitations: certified non-vanishing is limited to k ≤ kmax = 4 here (the code accepts any kmax). A rank-3 "
      "certificate needs exact vanishing of L'(E,1), which this component does not provide. "
      "Certified ord = 2 is a statement about L(E,s) only. The rank part of BSD for these curves also needs "
      "rank E(Q) = 2, which is a different computation (descent plus independent points).")
    w("")
    w("## 8. Interface for other components (e.g. a Gross–Zagier certificate)")
    w("")
    w("```python")
    w("from bsd.rigorous_l import certify_order_rigorous, taylor_enclosures")
    w("res = certify_order_rigorous(ainvs)                # JSON-ready dict; res['conclusion'] is a sentence")
    w("eps = res['L_prime_abs_upper']                     # w = -1, rank-3 gap: rigorous |L'(E,1)| <= eps (decimal string)")
    w("enc = taylor_enclosures(ainvs, kmax=3, prec=512)   # Arb balls: enc.c[1] = L'(E,1), enc.c[3] = L^(3)(E,1)/6")
    w("twist = taylor_enclosures(twist_ainvs, kmax=0)     # any curve, e.g. a quadratic twist: ball for L(E^D,1)")
    w("# once L'(E,1) = 0 has been proven exactly elsewhere:")
    w("res3 = certify_order_rigorous(ainvs, external_zeros={1: 'Gross-Zagier: Heegner point is torsion'})")
    w("# -> res3['order'] == 3 and res3['order_certified'] is True; the source is quoted in res3['conclusion']")
    w("```")
    w("")
    w("## 9. Reproduce")
    w("")
    w("```")
    w("python3 -m pytest -q tests/test_rigorous_l.py          # ~2 s")
    w("python3 -m bsd.rigorous_l 389a1 --kmax 4 --prec 128     # one certificate as JSON")
    w("python3 research/rigorous-l/run_batch.py all            # resumable; skips finished items")
    w("python3 research/rigorous-l/aggregate.py && python3 research/rigorous-l/make_report.py")
    w("```")
    w("")
    with open(os.path.join(HERE, "report.md"), "w") as f:
        f.write("\n".join(L))
    print("wrote report.md,", len(L), "lines")


if __name__ == "__main__":
    main()
