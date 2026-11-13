/* Level-2 data of the p-adic L-function of 389a1 at a prime p (ordinary or supersingular), in compiled code.

   For k in Z/p let A2(k) = sum_{i=1}^{p-1} m'( omega(i) gamma^k mod p^2 / p^2 ), omega(i) = i^p mod p^2 (Teichmueller),
   gamma = 1 + p. Then theta_1 = sum_k A2(k) (1+T)^k is the level-p^2 Mazur-Tate element (bsd/kato_bound.py), and the
   ordinary level-2 Riemann sum is sum_k [alpha^-2 A2(k) - alpha^-3 S1] (1+T)^k with S1 = sum_{i=1}^{p-1} m'(i/p).

   The work is organised by rows R_i(k) = m'(omega(i) gamma^k mod p^2 / p^2), k in Z/p. Since gamma^k = 1 + kp mod p^2,
   row i runs over all lifts of i mod p, so each row satisfies the Hecke relation (U_p, p not dividing N)
       sum_k R_i(k) = a_p m'(i/p) - m'(0),
   which is checked exactly for every computed row. Two symmetries reduce the rows to one per orbit of
   {i, -i, tau(i), -tau(i)}, tau(i) = -N^-1 i^-1 mod p:
     R_{-i}(k) = R_i(k)                        (the plus symbol is even, omega(-i) = -omega(i));
     R_{tau(i)}(k) = eps R_i(c - k)            (Fricke lemma, proved in ../level3/run_level3.py: m'(b/M) = -w_N m'(b'/M)
                                                for N b b' = -1 mod M; eps = -w_N, c = -log_gamma<N> mod p).
   For an orbit of size 2 (tau(i) = +-i) the second relation becomes the check R_i(k) = eps R_i(c - k).

   Output (one line):  p a_p c eps S1 m0 rows hecke_fail fricke_fail C_0..C_4 W_0..W_4 evals seconds kernel
   with C_j = sum_k binom(k, j) A2(k) and W_j = sum_k binom(k, j), both mod p^3, and kernel = "lean".

   msym (scalar) and run_blocks (vector) are copied verbatim from ../level3/level3.c, where they are validated.
   The rows are computed by run_blocks_lean (2026-11-13), the same continued-fraction arithmetic with fewer vector
   instructions per step (about 1.3 times faster); its exactness argument is at the function. `vtest` checks it
   against run_blocks and msym on random b/q, and its output line equals that of the run_blocks version at every
   prime compared (../data/beyond_lean_validation.json).
   Usage: level2 TABLE p a_p c eps
          level2 TABLE vtest q count seed     (lean kernel = run_blocks = msym on random b/q, q < 2^40)
   Build: gcc -O3 -march=native -o level2 level2.c   (-DREFERENCE_KERNEL computes the rows with run_blocks) */
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


/* Lean kernel: the arithmetic of run_blocks with fewer vector instructions per continued-fraction step. Valid for
   q < 2^40 (then every quotient is below 2^40 and the expansion of b/q has at most 59 steps).
   - Residues K are kept in [0, N], N standing for 0. F = floor(fl(T * fl(1/N))) has relative error below 2^-51 and
     T = A K1 + K2 < 2^49, so F = floor(T/N), except that F may be one less when N divides T. Hence Kn = T - F N lies
     in [0, N], is computed exactly, and no correction is needed. PHI8 has (N+1)^2 entries per parity, its row and
     column N repeating row and column 0.
   - All lanes of a group step together, so the parity (which of the two tables) is one scalar per group.
   - A finished lane is parked with Y = 1 (X becomes the gcd, 1): every later step has A = 1, R = 0 and an index in
     range, and its sum is no longer updated. Padding lanes hold b = 1 and are not stored.
   - PHI8 holds bytes. Each step adds the 4 bytes read at the index to a 32-bit lane sum, and only the low byte of
     that sum is used. It equals the sum of the table values mod 256. That sum has absolute value at most 2 * 59 < 128
     (|values| <= 2, checked when PHI8 is filled), so the signed low byte is the sum itself. */
#define N1 (N + 1)
#define GL 4
static int8_t PHI8[2 * N1 * N1 + 8];

static void fill_phi8(void) {
    for (int par = 0; par < 2; par++)
        for (int a = 0; a <= N; a++)
            for (int b = 0; b <= N; b++) {
                int v = PHI[par][a % N][b % N];
                if (v < -2 || v > 2) { fprintf(stderr, "table value outside [-2, 2]\n"); exit(2); }
                PHI8[par * N1 * N1 + a * N1 + b] = (int8_t)v;
            }
}

/* tasks tb[0..ntask) (padding up to 8 nblk holds 1.0), results added to acc[td[t]] */
static void run_blocks_lean(const double *tb, const int32_t *td, size_t nblk, size_t ntask, u64 q, i64 *acc) {
    const __m512d vq = _mm512_set1_pd((double)q), vN = _mm512_set1_pd((double)N), vN1 = _mm512_set1_pd((double)N1);
    const __m512d invN = _mm512_set1_pd(1.0 / N), one = _mm512_set1_pd(1.0), zero = _mm512_setzero_pd();
    const i64 s0 = PHI[0][0][1];
    __m512d X[GL], Y[GL], K1[GL], K2[GL];
    __m256i S[GL];
    __mmask8 ACT[GL], REAL[GL];
    int par[GL];
    size_t blk[GL], next = 0, live = 0;
    for (int g = 0; g < GL; g++) {
        X[g] = vq; K1[g] = one; K2[g] = zero; S[g] = _mm256_setzero_si256(); par[g] = 1;
        if (next < nblk) {
            blk[g] = next++; live++;
            Y[g] = _mm512_loadu_pd(tb + 8 * blk[g]);
            size_t left = ntask - 8 * blk[g];
            ACT[g] = REAL[g] = left >= 8 ? 0xFF : (__mmask8)((1u << left) - 1);
        } else { blk[g] = (size_t)-1; X[g] = one; Y[g] = one; ACT[g] = REAL[g] = 0; }
    }
    while (live) {
        for (int g = 0; g < GL; g++) {
            __m512d A = _mm512_roundscale_pd(_mm512_div_pd(X[g], Y[g]), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
            __m512d R = _mm512_fnmadd_pd(A, Y[g], X[g]);                 /* exact remainder */
            __m512d T = _mm512_fmadd_pd(A, K1[g], K2[g]);                /* exact, < 2^49 */
            __m512d F = _mm512_roundscale_pd(_mm512_mul_pd(T, invN), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
            __m512d Kn = _mm512_fnmadd_pd(F, vN, T);                     /* T mod N in [0, N], exact */
            __m256i idx = _mm512_cvttpd_epi32(_mm512_fmadd_pd(K1[g], vN1, Kn));
            const int8_t *base = PHI8 + (par[g] ? N1 * N1 : 0);
            S[g] = _mm256_mask_add_epi32(S[g], ACT[g], S[g], _mm256_i32gather_epi32((const int *)base, idx, 1));
            __mmask8 done = _mm512_cmp_pd_mask(R, zero, _CMP_EQ_OQ) & ACT[g];
            ACT[g] &= ~done;
            K2[g] = K1[g]; K1[g] = Kn; X[g] = Y[g]; Y[g] = _mm512_mask_mov_pd(one, ACT[g], R);   /* parked: Y = 1 */
            par[g] ^= 1;
            if (__builtin_expect(ACT[g] == 0 && blk[g] != (size_t)-1, 0)) {
                int32_t s[8];
                _mm256_storeu_si256((__m256i *)s, S[g]);
                for (int l = 0; l < 8; l++)
                    if (REAL[g] >> l & 1) acc[td[8 * blk[g] + l]] += s0 + (int8_t)(s[l] & 0xFF);
                X[g] = vq; K1[g] = one; K2[g] = zero; S[g] = _mm256_setzero_si256(); par[g] = 1;
                if (next < nblk) {
                    blk[g] = next++;
                    Y[g] = _mm512_loadu_pd(tb + 8 * blk[g]);
                    size_t left = ntask - 8 * blk[g];
                    ACT[g] = REAL[g] = left >= 8 ? 0xFF : (__mmask8)((1u << left) - 1);
                } else { blk[g] = (size_t)-1; X[g] = one; Y[g] = one; ACT[g] = REAL[g] = 0; live--; }
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

static u64 gcd64(u64 u, u64 v) { while (v) { u64 w = u % v; u = v; v = w; } return u; }

int main(int argc, char **argv) {
    if (argc != 6) { fprintf(stderr, "usage: level2 TABLE p a_p c eps | level2 TABLE vtest q count seed\n"); return 2; }
    load_table(argv[1]);
    fill_phi8();
    if (!strcmp(argv[2], "vtest")) {
        u64 q = strtoull(argv[3], 0, 10), n = strtoull(argv[4], 0, 10), x = strtoull(argv[5], 0, 10) | 1;
        if (q < 3 || q >= (1ULL << 40) || n == 0) { fprintf(stderr, "bad arguments\n"); return 2; }
        size_t nb = (n + 7) / 8, ndest = 4096;
        double *tb = malloc(nb * 8 * sizeof(double)), *tl = malloc(nb * 8 * sizeof(double));
        int32_t *td = malloc(nb * 8 * sizeof(int32_t));
        i64 *ref = calloc(ndest, sizeof(i64)), *acc = calloc(ndest, sizeof(i64)), *lean = calloc(ndest, sizeof(i64));
        for (size_t t = 0; t < nb * 8; t++) {
            if (t >= n) { tb[t] = NAN; tl[t] = 1.0; td[t] = 0; continue; }
            u64 b;
            do { x ^= x << 13; x ^= x >> 7; x ^= x << 17; b = x % q; } while (b == 0 || gcd64(b, q) != 1);
            tb[t] = tl[t] = (double)b; td[t] = (int32_t)(t % ndest);
            ref[t % ndest] += msym(b, q);
        }
        double t0 = now();
        run_blocks(tb, td, nb, q, acc);
        double t1 = now();
        run_blocks_lean(tl, td, nb, n, q, lean);
        double t2 = now();
        size_t bad = 0, badl = 0;
        for (size_t d = 0; d < ndest; d++) { bad += acc[d] != ref[d]; badl += lean[d] != ref[d]; }
        printf("vtest q=%llu n=%llu: run_blocks %.1f ns/eval, %zu destinations differ from msym; "
               "run_blocks_lean %.1f ns/eval, %zu differ\n", (unsigned long long)q, (unsigned long long)n,
               1e9 * (t1 - t0) / n, bad, 1e9 * (t2 - t1) / n, badl);
        return bad != 0 || badl != 0;
    }
    u64 p = strtoull(argv[2], 0, 10);
    i64 ap = strtoll(argv[3], 0, 10);
    u64 c = strtoull(argv[4], 0, 10);
    i64 eps = strtoll(argv[5], 0, 10);
    u64 p2 = p * p, p3 = p2 * p;
    if (p < 5 || p % 2 == 0 || p >= (1ULL << 20) || c >= p || (eps != 1 && eps != -1)) { fprintf(stderr, "bad arguments\n"); return 2; }
    double t0 = now();
    i64 m0 = PHI[0][0][1];                       /* m'(0): the path oo -> 0 is the first piece */
    i64 *mp = malloc(p * sizeof(i64));           /* m'(i/p) */
    i64 S1 = 0;
    for (u64 i = 1; i < p; i++) { mp[i] = msym(i, p); S1 += mp[i]; }
    u64 Ninv = powmod(N % p, p - 2, p);
    char *seen = calloc(p, 1);
    u64 *reps = malloc(p * sizeof(u64));
    int *full = malloc(p * sizeof(int));
    size_t nrep = 0;
    for (u64 i = 1; i < p; i++) {
        if (seen[i]) continue;
        u64 ti = (p - mulmod(Ninv, powmod(i, p - 2, p), p)) % p;   /* tau(i) = -N^-1 i^-1 */
        seen[i] = seen[p - i] = seen[ti] = seen[p - ti] = 1;
        reps[nrep] = i; full[nrep] = !(ti == i || ti == p - i); nrep++;
    }
    i64 *A2 = calloc(p, sizeof(i64));
    size_t RB = 256, cap = ((RB * p + 7) / 8) * 8;
    double *tb = malloc(cap * sizeof(double));
    int32_t *td = malloc(cap * sizeof(int32_t));
    i64 *acc = malloc(RB * p * sizeof(i64));
    long hecke_fail = 0, fricke_fail = 0;
    u64 evals = 0;
    for (size_t r0 = 0; r0 < nrep; r0 += RB) {
        size_t nr = nrep - r0 < RB ? nrep - r0 : RB, n = 0;
        memset(acc, 0, nr * p * sizeof(i64));
        for (size_t r = 0; r < nr; r++) {
            u64 i = reps[r0 + r], b = powmod(i, p, p2), step = (i * p) % p2;   /* omega(i) gamma^k = omega(i) + k i p */
            for (u64 k = 0; k < p; k++) {
                tb[n] = (double)b; td[n] = (int32_t)(r * p + k); n++;
                b += step; if (b >= p2) b -= p2;
            }
        }
        evals += n;
        size_t nb = (n + 7) / 8;
#ifdef REFERENCE_KERNEL
        for (size_t t = n; t < nb * 8; t++) { tb[t] = NAN; td[t] = 0; }
        run_blocks(tb, td, nb, p2, acc);
#else
        for (size_t t = n; t < nb * 8; t++) { tb[t] = 1.0; td[t] = 0; }
        run_blocks_lean(tb, td, nb, n, p2, acc);
#endif
        for (size_t r = 0; r < nr; r++) {
            u64 i = reps[r0 + r];
            i64 *R = acc + r * p, tot = 0;
            for (u64 k = 0; k < p; k++) tot += R[k];
            if (tot != ap * mp[i] - m0) hecke_fail++;
            if (full[r0 + r]) {
                for (u64 k = 0; k < p; k++) A2[k] += 2 * R[k] + 2 * eps * R[(c + p - k) % p];
            } else {
                for (u64 k = 0; k < p; k++) {
                    A2[k] += 2 * R[k];
                    if (R[k] != eps * R[(c + p - k) % p]) fricke_fail++;
                }
            }
        }
    }
    printf("%llu %lld %llu %lld %lld %lld %zu %ld %ld", (unsigned long long)p, (long long)ap, (unsigned long long)c,
           (long long)eps, (long long)S1, (long long)m0, nrep, hecke_fail, fricke_fail);
    u64 W[5];
    for (int j = 0; j <= 4; j++) {
        u128 s = 0, w = 0;
        for (u64 k = 0; k < p; k++) {
            i64 a = A2[k];
            u64 am = a >= 0 ? (u64)a % p3 : (p3 - (u64)(-a) % p3) % p3, bj = binom_mod(k, j, p3);
            s = (s + (u128)bj * am) % p3;
            w = (w + bj) % p3;
        }
        W[j] = (u64)w;
        printf(" %llu", (unsigned long long)s);
    }
    for (int j = 0; j <= 4; j++) printf(" %llu", (unsigned long long)W[j]);
#ifdef REFERENCE_KERNEL
    printf(" %llu %.2f run_blocks\n", (unsigned long long)evals, now() - t0);
#else
    printf(" %llu %.2f lean\n", (unsigned long long)evals, now() - t0);
#endif
    return 0;
}
