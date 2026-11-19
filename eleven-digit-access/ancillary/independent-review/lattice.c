/* Periodic census of the leading map K on the rational lattice (1/N)Z^5 ∩ D.
 * K(x) sorts (1, x_i, 1-x_i) and takes symmetric differences (Section 2).
 * On Q = N x this is exact integer arithmetic.  All periodic points are found
 * by functional-graph cycle detection.  Prints the periodic points relevant to
 * Lemma 7 / Lemma 3 / Section 4.2:
 *   tag I: all gaps > 0 and e >= 1/30   (must lie in ray / triangle / isolated orbit)
 *   tag O: all gaps > 0, g0 >= 1/6, e < 1/30 (must be center or eq.(two))
 *   tag L (optional 2nd arg): all gaps > 0, g0 < 1/6, e < 1/30 (Lemma 8, coefficient 4/7)
 * plus summary counts.  Usage: lattice N [L]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef int64_t i64;
static uint64_t binom[400][7];
static int N;
static uint64_t rank_of(const int *q) { uint64_t r = 0; for (int i = 0; i < 5; i++) r += binom[q[i] + 4 - i][5 - i]; return r; }
static void Kmap(const int *q, int *o) {
    int s[11]; s[0] = N;
    for (int i = 0; i < 5; i++) { s[1 + i] = q[i]; s[6 + i] = N - q[i]; }
    for (int i = 1; i < 11; i++) { int x = s[i], j = i - 1; while (j >= 0 && s[j] < x) { s[j + 1] = s[j]; j--; } s[j + 1] = x; }
    for (int i = 0; i < 5; i++) o[i] = s[i] - s[10 - i];
}
int main(int argc, char **argv) {
    N = atoi(argv[1]); int printL = argc > 2;
    for (int n = 0; n < 400; n++) { binom[n][0] = 1; for (int k = 1; k < 7; k++) binom[n][k] = n ? binom[n - 1][k - 1] + binom[n - 1][k] : 0; }
    uint64_t M = binom[N + 5][5];          /* values 0..N */
    uint8_t *st = calloc(M, 1);            /* 0 unseen, 1 on path, 2 done-transient, 3 periodic */
    int (*path)[5] = malloc(sizeof(int[5]) * 200000);
    uint64_t *pr = malloc(8 * 200000);
    int q[5];
    long nper = 0, nI = 0, nO = 0, nlow = 0, maxlen = 0;
    for (q[0] = 0; q[0] <= N; q[0]++) for (q[1] = 0; q[1] <= q[0]; q[1]++) for (q[2] = 0; q[2] <= q[1]; q[2]++)
    for (q[3] = 0; q[3] <= q[2]; q[3]++) for (q[4] = 0; q[4] <= q[3]; q[4]++) {
        uint64_t rk = rank_of(q);
        if (st[rk]) continue;
        int len = 0, cur[5]; memcpy(cur, q, sizeof cur); uint64_t c = rk;
        while (st[c] == 0) {
            st[c] = 1; memcpy(path[len], cur, sizeof cur); pr[len] = c; len++;
            int nx[5]; Kmap(cur, nx); memcpy(cur, nx, sizeof cur); c = rank_of(cur);
        }
        if (len > maxlen) maxlen = len;
        int cyc_start = len;
        if (st[c] == 1) { for (int k = 0; k < len; k++) if (pr[k] == c) { cyc_start = k; break; } }
        for (int k = 0; k < len; k++) st[pr[k]] = (k >= cyc_start) ? 3 : 2;
    }
    for (q[0] = 0; q[0] <= N; q[0]++) for (q[1] = 0; q[1] <= q[0]; q[1]++) for (q[2] = 0; q[2] <= q[1]; q[2]++)
    for (q[3] = 0; q[3] <= q[2]; q[3]++) for (q[4] = 0; q[4] <= q[3]; q[4]++) {
        if (st[rank_of(q)] != 3) continue;
        nper++;
        int G[6] = {N - q[0], q[0] - q[1], q[1] - q[2], q[2] - q[3], q[3] - q[4], q[4]};
        i64 s2 = 0; int pos = 1;
        for (int i = 0; i < 6; i++) { i64 t = 6 * G[i] - N; s2 += t * t; if (G[i] <= 0) pos = 0; }
        int ehigh = 5 * s2 >= (i64)6 * N * N;      /* e >= 1/30 */
        if (!ehigh) nlow++;
        if (pos && ehigh) { nI++; printf("I %d %d %d %d %d %d\n", N, q[0], q[1], q[2], q[3], q[4]); }
        if (pos && !ehigh && 6 * G[0] >= N) { nO++; printf("O %d %d %d %d %d %d\n", N, q[0], q[1], q[2], q[3], q[4]); }
        if (printL && pos && !ehigh && 6 * G[0] < N) printf("L %d %d %d %d %d %d\n", N, q[0], q[1], q[2], q[3], q[4]);
    }
    fprintf(stderr, "N=%d points=%llu periodic=%ld low-energy periodic=%ld tagI=%ld tagO=%ld max_path=%ld\n", N,
            (unsigned long long)M, nper, nlow, nI, nO, maxlen);
    return 0;
}
