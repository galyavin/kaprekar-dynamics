/* Exhaustive enumeration of ALL B^5 five-digit strings in base B (reviewer
 * code).  Uses only the definition: x -> desc(x) - asc(x).  For every
 * nonconstant x it computes the first arrival time T at a numerical cycle and
 * the cycle reached.  It then
 *   (a) checks the input-weight formula (18) against the count of strings per
 *       sorted pair (u,v) = (d1-d5, d2-d4);
 *   (b) prints, for every cycle, its period and the histogram of T;
 *   (c) for even B, compares the period-2 cycle's histogram with (7).
 * usage: strings_exhaustive B
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static long B, N;

static long step(long x, long *u, long *v) {
    long d[5];
    for (int k = 0; k < 5; k++) { d[k] = x % B; x /= B; }
    for (int i = 1; i < 5; i++) { long t = d[i]; int j = i - 1; while (j >= 0 && d[j] < t) { d[j + 1] = d[j]; j--; } d[j + 1] = t; }
    long desc = 0, asc = 0;
    for (int k = 0; k < 5; k++) desc = desc * B + d[k];
    for (int k = 4; k >= 0; k--) asc = asc * B + d[k];
    if (u) { *u = d[0] - d[4]; *v = d[1] - d[3]; }
    return desc - asc;
}

int main(int argc, char **argv) {
    B = atol(argv[1]);
    N = B * B * B * B * B;
    int32_t *T = malloc(sizeof(int32_t) * N);
    int32_t *cid = malloc(sizeof(int32_t) * N);
    int32_t *pp = malloc(sizeof(int32_t) * N);
    int64_t *stk = malloc(sizeof(int64_t) * 100000);
    for (long i = 0; i < N; i++) { T[i] = -1; pp[i] = -1; }
    int nc = 0; long crep[256]; int clen[256];
    for (long s = 0; s < N; s++) {
        if (T[s] >= 0) continue;
        long top = 0, x = s;
        while (T[x] < 0 && pp[x] < 0) { pp[x] = top; stk[top++] = x; x = step(x, 0, 0); }
        if (T[x] < 0) {
            long p0 = pp[x], rep = x; int len = 0;
            for (long k = p0; k < top; k++) { long y = stk[k]; T[y] = 0; cid[y] = nc; if (y < rep) rep = y; len++; }
            crep[nc] = rep; clen[nc] = len; nc++;
            top = p0;
        }
        for (long k = top - 1; k >= 0; k--) { long y = stk[k]; long z = step(y, 0, 0); T[y] = T[z] + 1; cid[y] = cid[z]; }
    }
    /* weights */
    long rep1 = 0; for (int k = 0; k < 5; k++) rep1 = rep1 * B + 1;   /* 11111_B */
    long *cnt = calloc((size_t)B * B, sizeof(long));
    long maxT = 0;
    for (long x = 0; x < N; x++) {
        if (x % rep1 == 0) continue;   /* constant strings c*11111_B */
        long u, v, d[5], y = x;
        for (int k = 0; k < 5; k++) { d[k] = y % B; y /= B; }
        for (int i = 1; i < 5; i++) { long t = d[i]; int j = i - 1; while (j >= 0 && d[j] < t) { d[j + 1] = d[j]; j--; } d[j + 1] = t; }
        u = d[0] - d[4]; v = d[1] - d[3];
        cnt[u * B + v]++;
        if (T[x] > maxT) maxT = T[x];
    }
    int wok = 1; long total = 0;
    for (long u = 1; u < B; u++) for (long v = 0; v <= u; v++) {
        long w = v == 0 ? (B - u) * (20 * u - 10) : (u == v ? (B - u) * (30 * u - 10) : (B - u) * (120 * v * (u - v) - 20));
        if (w != cnt[u * B + v]) wok = 0;
        total += cnt[u * B + v];
    }
    printf("B %ld strings %ld nonconstant %ld weights18_ok %d maxT %ld cycles %d\n", B, N, total, wok, maxT, nc - 1);
    for (int c = 0; c < nc; c++) {
        if (crep[c] == 0) continue;   /* the constant-input fixed point 0 */
        long h[64] = {0}; long tot = 0; long tmax = 0;
        for (long x = 0; x < N; x++) if (x % rep1 && cid[x] == c) { if (T[x] < 64) h[T[x]]++; tot++; if (T[x] > tmax) tmax = T[x]; }
        long u, v; step(crep[c], &u, &v);
        long u0 = 0, v0 = 0; { long d[5], y = crep[c]; for (int k = 0; k < 5; k++) { d[k] = y % B; y /= B; }
            for (int i = 1; i < 5; i++) { long t = d[i]; int j = i - 1; while (j >= 0 && d[j] < t) { d[j + 1] = d[j]; j--; } d[j + 1] = t; }
            u0 = d[0] - d[4]; v0 = d[1] - d[3]; }
        printf("  cycle period %d minval %ld pair(min) (%ld,%ld) inputs %ld frac %.6e maxT %ld hist:", clen[c], crep[c], u0, v0, tot, (double)tot / total, tmax);
        for (int t = 0; t <= (tmax < 12 ? tmax : 12); t++) printf(" %ld", h[t]);
        printf("\n");
        if (clen[c] == 2 && B % 2 == 0) {
            long b2 = B * B;
            long f0 = 2, f1 = 35 * b2 - 75 * B - 12, f2 = 5 * b2 - 5 * B;
            int match = (h[0] == f0 && h[1] == f1 && h[2] == f2 && tot == f0 + f1 + f2);
            printf("  T2-formula(7): N0 %ld N1 %ld N2 %ld total %ld | observed N0 %ld N1 %ld N2 %ld total %ld | match %d\n",
                   f0, f1, f2, f0 + f1 + f2, h[0], h[1], h[2], tot, match);
        }
    }
    return 0;
}
