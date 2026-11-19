/* Brute-force phase work w(p) from its definition (Section 4.1):
 *   w(p) = lim_R limsup_B sup_{q physical, |q - Bp|_inf <= R} (E_B(F_B q) - E_B(q))/B,
 * with E_B(q) = sum_i (G_i - B/6)^2, so Delta E_B = sum G'^2 - sum G^2 (integers).
 * For p = P/N, B p mod 1 depends only on B mod N, so the limsup over all bases is
 * approximated by the max over one large base in every residue class mod lcm(N,2).
 * F_B is positional subtraction (kap.h).
 * Usage: work N R B0scale P1 P2 P3 P4 P5
 * Prints per phase: approx w, the maximising residue/offset, then W and W/(h*|y|_inf).
 */
#include <stdio.h>
#include <math.h>
#include "kap.h"

static int N;
static void Kmap(const i64 *q, i64 *o) {
    i64 s[11]; s[0] = N;
    for (int i = 0; i < 5; i++) { s[1 + i] = q[i]; s[6 + i] = N - q[i]; }
    sort_desc(s, 11);
    for (int i = 0; i < 5; i++) o[i] = s[i] - s[10 - i];
}
static i64 floordiv(i64 a, i64 b) { i64 q = a / b; if ((a % b != 0) && ((a < 0) != (b < 0))) q--; return q; }

int main(int argc, char **argv) {
    N = atoi(argv[1]); int R = atoi(argv[2]); i64 scale = atoll(argv[3]);
    i64 P[5]; for (int i = 0; i < 5; i++) P[i] = atoll(argv[4 + i]);
    /* leading orbit */
    i64 orb[800][5]; int h = 0;
    memcpy(orb[0], P, sizeof P);
    for (h = 1; h < 800; h++) {
        Kmap(orb[h - 1], orb[h]);
        if (!memcmp(orb[h], orb[0], sizeof P)) break;
    }
    if (h >= 800) { printf("not periodic\n"); return 1; }
    int M = (N % 2) ? 2 * N : N;   /* residues: lcm(N,2) */
    double W = 0;
    printf("N=%d h=%d R=%d residues=%d\n", N, h, R, M);
    for (int ph = 0; ph < h; ph++) {
        i64 *p = orb[ph];
        double best = -1e300; int bres = -1; i64 bq[5] = {0}; i64 bB = 0;
        for (int r = 0; r < M; r++) {
            i64 B = scale * M + r;
            /* integer box: q_i in [ceil(B p_i - R), floor(B p_i + R)] */
            i64 lo[5], hi[5];
            for (int i = 0; i < 5; i++) {
                __int128 c = (__int128)B * p[i];             /* N * (B p_i) */
                lo[i] = -floordiv((i64)(-(c - (__int128)R * N)), N);   /* ceil((c - RN)/N) */
                hi[i] = floordiv((i64)(c + (__int128)R * N), N);
            }
            i64 q[5];
            for (q[0] = lo[0]; q[0] <= hi[0]; q[0]++)
            for (q[1] = lo[1]; q[1] <= hi[1]; q[1]++) { if (q[1] > q[0]) break;
            for (q[2] = lo[2]; q[2] <= hi[2]; q[2]++) { if (q[2] > q[1]) break;
            for (q[3] = lo[3]; q[3] <= hi[3]; q[3]++) { if (q[3] > q[2]) break;
            for (q[4] = lo[4]; q[4] <= hi[4]; q[4]++) { if (q[4] > q[3]) break;
                if (q[0] >= B || q[4] < 0) continue;
                i64 qn[5], G[6], Gn[6];
                F_state(q, qn, B);
                gaps(q, B, G); gaps(qn, B, Gn);
                __int128 d = 0;
                for (int i = 0; i < 6; i++) d += (__int128)Gn[i] * Gn[i] - (__int128)G[i] * G[i];
                double v = (double)d / (double)B;
                if (v > best) { best = v; bres = r; memcpy(bq, q, sizeof q); bB = B; }
            }}}}
        }
        W += best;
        printf(" phase %d p=(%lld,%lld,%lld,%lld,%lld)/%d  w~%.6f  at B mod %d=%d offset=(", ph,
               (long long)p[0], (long long)p[1], (long long)p[2], (long long)p[3], (long long)p[4], N, best, M, bres);
        for (int i = 0; i < 5; i++) printf("%s%.3f", i ? "," : "", (double)bq[i] - (double)bB * p[i] / N);
        printf(")\n");
    }
    i64 G0[6] = {N - P[0], P[0] - P[1], P[1] - P[2], P[2] - P[3], P[3] - P[4], P[4]};
    double yinf = 0, e = 0;
    for (int i = 0; i < 6; i++) { double y = (double)G0[i] / N - 1.0 / 6; if (fabs(y) > yinf) yinf = fabs(y); e += y * y; }
    printf("W~%.6f  h|y|inf=%.6f  W/(h|y|inf)=%.6f  e=%.6f (s*=0.033333) g0=%.6f\n", W, h * yinf, W / (h * yinf), e, (double)G0[0] / N);
    return 0;
}
