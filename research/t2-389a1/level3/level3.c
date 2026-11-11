/* Level-3 Riemann sum of the p-adic L-function of 389a1 at a good ordinary prime p, in compiled code.

   Notation as in bsd/kato_bound.py (riemann_coefficients with n = 3): q = p^3, gamma = 1 + p, omega(i) = i^(p^2)
   mod q (Teichmueller), m'(r) the primitive integral plus modular symbol of 389a1 at {oo -> r}. For k in Z/p^2,
       A(k)  = sum_{i=1}^{p-1} m'( omega(i) gamma^k mod q / q ),
       S2(k) = sum_{i=1}^{p-1} m'( omega(i) gamma^k mod p^2 / p^2 )     (depends on k mod p only),
   and the level-3 Riemann sum is sum_k [alpha^-3 A(k) - alpha^-4 S2(k)] (1+T)^k. Writing k = k0 + p k1, the
   p values b = omega(i) gamma^(k0 + p k1) mod q, k1 = 0..p-1, are b0 + k1 i p^2 with b0 = omega(i) gamma^k0 mod q,
   since gamma^p = 1 + p^2 mod p^3. Evenness of the plus symbol (m'(-r) = m'(r), omega(p-i) = -omega(i)) halves
   the sum over i.

   m'(b/q) is evaluated by the continued fraction of b/q: the path oo -> b/q is the sum of the unimodular paths
   between consecutive convergents h_{i-1}/k_{i-1} -> h_i/k_i (determinant (-1)^i), whose values are read from
   PHI[i mod 2][k_{i-1} mod 389][k_i mod 389] (the table, built from PARI's mseval, holds the value of the path
   g.oo -> g.0 for g in SL_2(Z) with bottom row (c, (-1)^s d)). Many fractions are evaluated at once with AVX-512:
   8 lanes per vector, G vectors in flight, each running 8 tasks in lockstep. All arithmetic is exact: the
   integers involved stay below 2^53, and floor(fl(x/y)) = floor(x/y) for x < 2^53.

   For every k0 this program checks the Hecke relation (U_p, a_p = coefficient of 389a1)
       sum_{k1} A(k0 + p k1) = a_p S2(k0) - S1,      S1 = sum_{i=1}^{p-1} m'(i/p),
   which holds exactly because the b above run over all lifts of b0 mod p^2. It prints, per k0,
       k0 S2(k0) hecke_ok C_0 .. C_4 W_0 .. W_4 | k0r S2(k0r) hecke_r R_0 .. R_4 Wr_0 .. Wr_4 | seconds
   with C_j = sum_{k1} binom(k0 + p k1, j) A(k0 + p k1) and W_j = sum_{k1} binom(k0 + p k1, j), both mod p^4, so
   that the T^j coefficient of the level-3 sum is sum_{k0} [alpha^-3 C_j - alpha^-4 S2(k0) W_j] mod p^4.

   Reflection (Fricke symmetry, proof in run_level3.py): m'(b/M) = -w_N m'(b'/M) when N b b' = -1 mod M, hence
   A(c - k) = -w_N A(k) and S2(c - k) = -w_N S2(k) with c = -log_gamma<N> mod p^2. So the class k0r = c - k0 mod p
   is obtained from class k0 without new level-3 evaluations: R_j = sum_{k1} binom(r(k0 + p k1), j) A(k0 + p k1)
   with r(k) = (c - k) mod p^2 (the factor -w_N is applied by the caller), Wr_j its weights, and S2(k0r) is
   evaluated directly (level 2, cheap); hecke_r checks sum_{k1} A = a_p S2(k0r) - S1 for the reflected class.

   Usage: level3 TABLE p a_p c k0_start count           (k0 = k0_start + t mod p, t < count; one line per k0)
          level3 TABLE test                            (stdin: "b q" lines; stdout: m'(b/q), scalar)
          level3 TABLE vtest q count seed              (vector kernel vs scalar on random b, exact)
   Build: gcc -O3 -march=native -o level3 level3.c */
#include <immintrin.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 389
#define G 3
typedef uint64_t u64;
typedef int64_t i64;
typedef unsigned __int128 u128;

static int8_t PHI[2][N][N];
static int32_t PHI32[2 * N * N];

static void load_table(const char *fn) {
    FILE *f = fopen(fn, "rb");
    if (!f) { perror(fn); exit(2); }
    int32_t *buf = malloc(sizeof(int32_t) * 2 * N * N);
    if (fread(buf, sizeof(int32_t), 2 * N * N, f) != 2 * N * N) { fprintf(stderr, "short table\n"); exit(2); }
    fclose(f);
    for (int t = 0; t < 2 * N * N; t++) {
        if (buf[t] < -100 || buf[t] > 100) { fprintf(stderr, "table value out of range\n"); exit(2); }
        PHI32[t] = buf[t];
        PHI[t / (N * N)][(t / N) % N][t % N] = (int8_t)buf[t];
    }
    free(buf);
}

/* scalar reference: m'(b/q) for 0 < b < q, gcd(b, q) = 1 */
static int msym(u64 b, u64 q) {
    int s = PHI[0][0][1];
    u64 x = q, y = b;
    unsigned km2 = 0, km1 = 1, par = 1;
    while (y) {
        u64 a = x / y, r = x - a * y;
        unsigned am = a < N ? (unsigned)a : (unsigned)(a % N);
        unsigned k = (am * km1 + km2) % N;
        s += PHI[par][km1][k];
        km2 = km1; km1 = k; par ^= 1; x = y; y = r;
    }
    return s;
}

/* vector kernel: tasks tb[8m + l] (NaN = padding), results added to acc[td[8m + l]] */
static void run_blocks(const double *tb, const int32_t *td, size_t nblk, u64 q, i64 *acc) {
    const __m512d vq = _mm512_set1_pd((double)q), vN = _mm512_set1_pd((double)N), invN = _mm512_set1_pd(1.0 / N);
    const __m512d one = _mm512_set1_pd(1.0), zero = _mm512_setzero_pd(), vnan = _mm512_set1_pd(NAN);
    const __m512i vNN = _mm512_set1_epi64((i64)N * N), S0 = _mm512_set1_epi64(PHI[0][0][1]);
    __m512d X[G], Y[G], K1[G], K2[G];
    __m512i PO[G], S[G];
    __mmask8 ACT[G], REAL[G];
    size_t blk[G], next = 0, live = 0;
    for (int g = 0; g < G; g++) {
        X[g] = vq; K1[g] = one; K2[g] = zero; PO[g] = vNN; S[g] = S0;
        if (next < nblk) {
            blk[g] = next++; live++;
            Y[g] = _mm512_loadu_pd(tb + 8 * blk[g]);
            ACT[g] = REAL[g] = _mm512_cmp_pd_mask(Y[g], Y[g], _CMP_ORD_Q);
        } else { blk[g] = (size_t)-1; Y[g] = vnan; ACT[g] = REAL[g] = 0; }
    }
    while (live) {
        for (int g = 0; g < G; g++) {
            __m512d A = _mm512_roundscale_pd(_mm512_div_pd(X[g], Y[g]), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
            __m512d R = _mm512_fnmadd_pd(A, Y[g], X[g]);                 /* exact remainder */
            __m512d T = _mm512_fmadd_pd(A, K1[g], K2[g]);                /* exact, < 2^53 */
            __m512d F = _mm512_roundscale_pd(_mm512_mul_pd(T, invN), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
            __m512d Kn = _mm512_fnmadd_pd(F, vN, T);                     /* T mod N, up to one correction */
            Kn = _mm512_mask_add_pd(Kn, _mm512_cmp_pd_mask(Kn, zero, _CMP_LT_OQ), Kn, vN);
            Kn = _mm512_mask_sub_pd(Kn, _mm512_cmp_pd_mask(Kn, vN, _CMP_GE_OQ), Kn, vN);
            __m512i idx = _mm512_add_epi64(_mm512_cvttpd_epi64(_mm512_fmadd_pd(K1[g], vN, Kn)), PO[g]);
            idx = _mm512_maskz_mov_epi64(ACT[g], idx);
            S[g] = _mm512_mask_add_epi64(S[g], ACT[g], S[g], _mm512_cvtepi32_epi64(_mm512_i64gather_epi32(idx, PHI32, 4)));
            __mmask8 done = _mm512_cmp_pd_mask(R, zero, _CMP_EQ_OQ) & ACT[g];
            ACT[g] &= ~done;
            K2[g] = K1[g]; K1[g] = Kn; X[g] = Y[g]; Y[g] = _mm512_mask_mov_pd(R, done, vnan);
            PO[g] = _mm512_sub_epi64(vNN, PO[g]);
            if (__builtin_expect(ACT[g] == 0 && blk[g] != (size_t)-1, 0)) {
                i64 s[8];
                _mm512_storeu_si512(s, S[g]);
                for (int l = 0; l < 8; l++)
                    if (REAL[g] >> l & 1) acc[td[8 * blk[g] + l]] += s[l];
                X[g] = vq; K1[g] = one; K2[g] = zero; PO[g] = vNN; S[g] = S0;
                if (next < nblk) {
                    blk[g] = next++;
                    Y[g] = _mm512_loadu_pd(tb + 8 * blk[g]);
                    ACT[g] = REAL[g] = _mm512_cmp_pd_mask(Y[g], Y[g], _CMP_ORD_Q);
                } else { blk[g] = (size_t)-1; Y[g] = vnan; ACT[g] = REAL[g] = 0; live--; }
            }
        }
    }
}

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1 % m;
    a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}
static u64 binom_mod(u64 k, int j, u64 m) {   /* binom(k, j) mod m for j <= 4, k < 2^30 */
    if ((u64)j > k) return 0;
    u128 num = 1;
    for (int t = 0; t < j; t++) num *= (u128)(k - t);
    static const unsigned fact[5] = {1, 1, 2, 6, 24};
    return (u64)((num / fact[j]) % m);
}

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: see header\n"); return 2; }
    load_table(argv[1]);
    if (!strcmp(argv[2], "test")) {
        unsigned long long b, q;
        while (scanf("%llu %llu", &b, &q) == 2) printf("%d\n", msym(b % q, q));
        return 0;
    }
    if (!strcmp(argv[2], "vtest")) {
        u64 q = strtoull(argv[3], 0, 10), n = strtoull(argv[4], 0, 10), x = strtoull(argv[5], 0, 10) | 1;
        size_t nb = (n + 7) / 8, ndest = 4096;
        double *tb = malloc(nb * 8 * sizeof(double));
        int32_t *td = malloc(nb * 8 * sizeof(int32_t));
        i64 *acc = calloc(ndest, sizeof(i64)), *ref = calloc(ndest, sizeof(i64));
        for (size_t t = 0; t < nb * 8; t++) {
            if (t >= n) { tb[t] = NAN; td[t] = 0; continue; }
            u64 b;
            do { x ^= x << 13; x ^= x >> 7; x ^= x << 17; b = x % q; } while (b == 0 || ({ u64 u = b, v = q; while (v) { u64 w = u % v; u = v; v = w; } u; }) != 1);
            tb[t] = (double)b; td[t] = (int32_t)(t % ndest);
            ref[t % ndest] += msym(b, q);
        }
        double t0 = now();
        run_blocks(tb, td, nb, q, acc);
        double sec = now() - t0;
        size_t bad = 0;
        for (size_t d = 0; d < ndest; d++) bad += acc[d] != ref[d];
        printf("vtest q=%llu n=%llu: %.1f ns/eval, destinations differing: %zu\n", (unsigned long long)q,
               (unsigned long long)n, 1e9 * sec / n, bad);
        return bad != 0;
    }
    u64 p = strtoull(argv[2], 0, 10);
    i64 ap = strtoll(argv[3], 0, 10);
    u64 cref = strtoull(argv[4], 0, 10), k0s = strtoull(argv[5], 0, 10), cnt = strtoull(argv[6], 0, 10);
    u64 p2 = p * p, q = p2 * p, p4 = q * p, h = (p - 1) / 2;
    if (p < 5 || p % 2 == 0 || q >= (1ULL << 52) || k0s >= p || cnt > p || cref >= p2) { fprintf(stderr, "bad arguments\n"); return 2; }
    u64 *om = malloc((h + 1) * sizeof(u64));
    for (u64 i = 1; i <= h; i++) om[i] = powmod(i, p2, q);
    i64 S1 = 0;
    for (u64 i = 1; i <= h; i++) S1 += msym(i, p);
    S1 *= 2;
    size_t ipb = 64;                                           /* values of i per task buffer */
    size_t cap = ((ipb * p + 7) / 8) * 8;
    double *tb = malloc(cap * sizeof(double));
    int32_t *td = malloc(cap * sizeof(int32_t));
    i64 *acc = malloc(p * sizeof(i64));
    for (u64 t = 0; t < cnt; t++) {
        u64 k0 = (k0s + t) % p, g0 = powmod(1 + p, k0, q);
        double t0 = now();
        memset(acc, 0, p * sizeof(i64));
        i64 S2 = 0;
        for (u64 i0 = 1; i0 <= h; i0 += ipb) {
            size_t n = 0;
            for (u64 i = i0; i < i0 + ipb && i <= h; i++) {
                u64 b0 = mulmod(om[i], g0, q), step = mulmod(i, p2, q), b = b0;
                S2 += msym(b0 % p2, p2);
                for (u64 k1 = 0; k1 < p; k1++) {
                    tb[n] = (double)b; td[n] = (int32_t)k1; n++;
                    b += step; if (b >= q) b -= q;
                }
            }
            size_t nb = (n + 7) / 8;
            for (size_t t = n; t < nb * 8; t++) { tb[t] = NAN; td[t] = 0; }
            run_blocks(tb, td, nb, q, acc);
        }
        S2 *= 2;
        i64 tot = 0;
        for (u64 k1 = 0; k1 < p; k1++) { acc[k1] *= 2; tot += acc[k1]; }
        int hecke_ok = tot == ap * S2 - S1;
        printf("%llu %lld %d", (unsigned long long)k0, (long long)S2, hecke_ok);
        u64 W[5];
        for (int j = 0; j <= 4; j++) {
            u128 c = 0, w = 0;
            for (u64 k1 = 0; k1 < p; k1++) {
                i64 a = acc[k1];
                u64 am = a >= 0 ? (u64)a % p4 : (p4 - (u64)(-a) % p4) % p4, bj = binom_mod(k0 + p * k1, j, p4);
                c = (c + (u128)bj * am) % p4;
                w = (w + bj) % p4;
            }
            W[j] = (u64)w;
            printf(" %llu", (unsigned long long)c);
        }
        for (int j = 0; j <= 4; j++) printf(" %llu", (unsigned long long)W[j]);
        /* reflected class k0r = c - k0 mod p */
        u64 k0r = (cref % p + p - k0) % p, gr = powmod(1 + p, k0r, q);
        i64 S2r = 0;
        for (u64 i = 1; i <= h; i++) S2r += msym(mulmod(om[i], gr, q) % p2, p2);
        S2r *= 2;
        printf(" %llu %lld %d", (unsigned long long)k0r, (long long)S2r, (int)(tot == ap * S2r - S1));
        u64 Wr[5];
        for (int j = 0; j <= 4; j++) {
            u128 c = 0, w = 0;
            for (u64 k1 = 0; k1 < p; k1++) {
                i64 a = acc[k1];
                u64 am = a >= 0 ? (u64)a % p4 : (p4 - (u64)(-a) % p4) % p4;
                u64 rk = (cref + p2 - (k0 + p * k1) % p2) % p2, bj = binom_mod(rk, j, p4);
                c = (c + (u128)bj * am) % p4;
                w = (w + bj) % p4;
            }
            Wr[j] = (u64)w;
            printf(" %llu", (unsigned long long)c);
        }
        for (int j = 0; j <= 4; j++) printf(" %llu", (unsigned long long)Wr[j]);
        printf(" %.2f\n", now() - t0);
        fflush(stdout);
    }
    return 0;
}
