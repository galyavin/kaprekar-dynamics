/* Random and adversarial tests of Lemma 12 (centered lower hits) with own Psi (eq.(centered)).
 * X in Z^5 with X_i = i*beta (mod 6); R(X) = -min ghat(X).  Claim: if R(X) > 120 then some
 * 1<=j<=32 has R(Psi^j X) <= R(X) - 6.  Records the largest first-hit delay found, by uniform
 * sampling at many scales and by greedy local search maximising the delay.
 * Usage: lowerhit nrandom nclimb seed
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "kap.h"
static void psi(const i64 *X, i64 *Y) {
    i64 c[11]; c[0] = -6;
    for (int i = 1; i <= 5; i++) {
        i64 a = X[i - 1] - 6 * (i == 5), b = -X[6 - i - 1] - 6 * (i < 5);
        c[2 * i - 1] = a > b ? a : b; c[2 * i] = a > b ? b : a;
    }
    for (int k = 1; k <= 5; k++) Y[k - 1] = c[k - 1] - c[11 - k];
}
static i64 Rof(const i64 *X) {
    i64 g[6] = {-X[0], X[0] - X[1], X[1] - X[2], X[2] - X[3], X[3] - X[4], X[4]};
    i64 m = g[0]; for (int i = 1; i < 6; i++) if (g[i] < m) m = g[i];
    return -m;
}
/* first-hit delay, 0 if R<=120, 99 if none within 32 */
static int delay(const i64 *X0) {
    i64 R0 = Rof(X0); if (R0 <= 120) return 0;
    i64 X[5], Y[5]; memcpy(X, X0, sizeof X);
    for (int j = 1; j <= 32; j++) { psi(X, Y); memcpy(X, Y, sizeof X); if (Rof(X) <= R0 - 6) return j; }
    return 99;
}
static void randX(rng_t *r, i64 *X, i64 scale) {
    int beta = rng_below(r, 6);
    for (int i = 0; i < 5; i++) X[i] = 6 * (rng_below(r, 2 * scale + 1) - scale) + (i + 1) * beta;
}
int main(int argc, char **argv) {
    long nr = atol(argv[1]), nc = atol(argv[2]); uint64_t seed = strtoull(argv[3], 0, 10);
    rng_t r; rng_seed(&r, seed);
    long hist[100] = {0}, fails = 0, tested = 0; int best = 0; i64 bestX[5] = {0};
    i64 scales[] = {25, 40, 60, 100, 300, 1000, 10000, 1000000};
    for (long t = 0; t < nr; t++) {
        i64 X[5]; randX(&r, X, scales[t % 8]);
        int d = delay(X); if (!d) continue;
        tested++; hist[d]++; if (d == 99) fails++;
        if (d > best && d != 99) { best = d; memcpy(bestX, X, sizeof X); }
    }
    /* greedy local search for long delays */
    int cbest = 0; i64 cX[5] = {0}; long cfails = 0;
    for (long c = 0; c < nc; c++) {
        i64 X[5]; randX(&r, X, scales[c % 6]);
        int d = delay(X); if (!d) continue;
        for (int it = 0; it < 3000; it++) {
            i64 Z[5]; memcpy(Z, X, sizeof Z);
            int k = rng_below(&r, 5); i64 step = 6 * (rng_below(&r, 7) - 3);
            if (rng_below(&r, 4) == 0) { int beta = ((Z[0] % 6) + 6) % 6; (void)beta; }
            Z[k] += step;
            int dz = delay(Z);
            if (dz == 99) { cfails++; d = dz; memcpy(X, Z, sizeof X); break; }
            if (dz >= d && dz > 0) { d = dz; memcpy(X, Z, sizeof X); }
        }
        if (d > cbest) { cbest = d; memcpy(cX, X, sizeof X); }
    }
    printf("random: tested(R>120)=%ld failures(no hit in 32)=%ld max delay=%d at X=(%lld,%lld,%lld,%lld,%lld) R=%lld\n", tested, fails, best,
           (long long)bestX[0], (long long)bestX[1], (long long)bestX[2], (long long)bestX[3], (long long)bestX[4], (long long)Rof(bestX));
    printf("delay histogram:"); for (int d = 1; d < 100; d++) if (hist[d]) printf(" %d:%ld", d, hist[d]); printf("\n");
    printf("local search: climbs=%ld failures=%ld max delay=%d at X=(%lld,%lld,%lld,%lld,%lld) R=%lld\n", nc, cfails, cbest,
           (long long)cX[0], (long long)cX[1], (long long)cX[2], (long long)cX[3], (long long)cX[4], (long long)Rof(cX));
    return 0;
}
