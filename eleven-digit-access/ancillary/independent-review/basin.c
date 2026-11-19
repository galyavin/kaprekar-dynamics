/* Complete inverse basin of the B=5m exceptional four-cycle, by brute force:
 * every nonzero difference state B>q1>=...>=q5>=0 is classified by following
 * its forward F_B orbit (positional subtraction, kap.h) with memoisation.
 * No depth bound and no inverse rule are assumed.
 * Usage: basin m      (B = 5m)
 * Output: layer sizes by distance, min-gap multiplicities, sum(minG+1),
 *         exact sorted-word count, cycle check.
 */
#include <stdio.h>
#include "kap.h"

static i64 B;
static uint64_t binom[300][7];

static uint64_t rank_of(const i64 *q) {      /* combinatorial number system */
    uint64_t r = 0;
    for (int i = 0; i < 5; i++) r += binom[q[i] + 4 - i][5 - i];
    return r;
}

int main(int argc, char **argv) {
    i64 m = atoll(argv[1]); B = 5 * m;
    for (int n = 0; n < 300; n++) { binom[n][0] = 1; for (int k = 1; k < 7; k++) binom[n][k] = n ? binom[n - 1][k - 1] + binom[n - 1][k] : 0; }
    uint64_t N = binom[B + 4][5];
    uint8_t *st = calloc(N, 1);           /* 0 unknown, 1 on stack, 2 outside basin, 3+d distance d */
    if (!st) { fprintf(stderr, "oom\n"); return 1; }
    /* the claimed cycle */
    i64 cyc[4][5] = {
        {4 * m + 1, 3 * m, 3 * m - 1, 2 * m, m + 2},
        {4 * m, 3 * m, 2 * m - 2, m, m - 1},
        {4 * m + 1, 3 * m, 3 * m, 2 * m + 1, m + 2},
        {4 * m, 3 * m, 2 * m - 2, m + 1, m - 1}};
    int cyc_ok = 1;
    for (int k = 0; k < 4; k++) {
        i64 nx[5]; F_state(cyc[k], nx, B);
        if (memcmp(nx, cyc[(k + 1) % 4], sizeof nx)) cyc_ok = 0;
        st[rank_of(cyc[k])] = 3;
    }
    printf("B=%lld m=%lld states=%llu cycle_closes=%d\n", (long long)B, (long long)m, (unsigned long long)N, cyc_ok);
    i64 (*path)[5] = malloc(sizeof(i64[5]) * 1000000);
    uint64_t *prk = malloc(sizeof(uint64_t) * 1000000);
    i64 q[5];
    /* enumerate all states */
    for (q[0] = 0; q[0] < B; q[0]++)
    for (q[1] = 0; q[1] <= q[0]; q[1]++)
    for (q[2] = 0; q[2] <= q[1]; q[2]++)
    for (q[3] = 0; q[3] <= q[2]; q[3]++)
    for (q[4] = 0; q[4] <= q[3]; q[4]++) {
        uint64_t rk = rank_of(q);
        if (st[rk]) continue;
        if (q[0] == 0) { st[rk] = 2; continue; }     /* zero state: constant words, excluded */
        int len = 0; i64 cur[5]; memcpy(cur, q, sizeof cur); uint64_t crk = rk;
        int res;
        for (;;) {
            if (st[crk] == 1) { res = 2; break; }            /* closed a new cycle: not the target */
            if (st[crk] >= 2) { res = st[crk]; break; }
            st[crk] = 1; memcpy(path[len], cur, sizeof cur); prk[len] = crk; len++;
            if (len >= 1000000) { fprintf(stderr, "path overflow\n"); return 1; }
            i64 nx[5]; F_state(cur, nx, B);
            memcpy(cur, nx, sizeof cur); crk = rank_of(cur);
        }
        for (int k = len - 1; k >= 0; k--) {
            if (res == 2) st[prk[k]] = 2;
            else { res = res + 1; if (res > 255) { fprintf(stderr, "depth overflow\n"); return 1; } st[prk[k]] = res; }
        }
    }
    /* tally */
    long layer[260] = {0}, mg[200] = {0}; long tot = 0; long sum_mg1 = 0;
    __int128 words = 0;
    for (q[0] = 0; q[0] < B; q[0]++)
    for (q[1] = 0; q[1] <= q[0]; q[1]++)
    for (q[2] = 0; q[2] <= q[1]; q[2]++)
    for (q[3] = 0; q[3] <= q[2]; q[3]++)
    for (q[4] = 0; q[4] <= q[3]; q[4]++) {
        uint8_t s = st[rank_of(q)];
        if (s >= 3) {
            layer[s - 3]++; tot++;
            i64 G[6]; gaps(q, B, G);
            i64 mn = min_gap(q, B); mg[mn]++; sum_mg1 += mn + 1;
            __int128 c = G[0]; for (int i = 1; i < 6; i++) c *= (G[i] + 1);
            words += c;
            if (layer[s - 3] <= 3 && s - 3 >= 7) printf("  deep state d=%d q=(%lld,%lld,%lld,%lld,%lld) minG=%lld\n", s - 3,
                (long long)q[0], (long long)q[1], (long long)q[2], (long long)q[3], (long long)q[4], (long long)mn);
        }
    }
    printf("basin size=%ld layers:", tot);
    for (int d = 0; d < 260; d++) if (layer[d]) printf(" d%d:%ld", d, layer[d]);
    printf("\nmin-gap multiplicities:");
    for (int k = 0; k < 200; k++) if (mg[k]) printf(" %d:%ld", k, mg[k]);
    printf("\nsum(minG+1)=%ld  sorted-word count=%.6e  bound 391*B^5=%.6e  ratio=%.4f\n", sum_mg1, (double)words,
           391.0 * (double)B * B * B * B * B, (double)words / (391.0 * (double)B * B * B * B * B));
    return 0;
}
