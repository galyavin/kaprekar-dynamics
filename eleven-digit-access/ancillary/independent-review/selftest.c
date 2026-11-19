/* Exact integer tests of the elementary identities of Sections 1-2:
 *  (a) eq.(digits) multiset vs positional subtraction (q5>0)
 *  (b) eq.(track): |F_B(q) - B K(q/B)|_inf <= 1, gap error <= 2   (q5>0)
 *  (c) eq.(outer): G0(F q) = min(G0+1-[q1=q2], q5)                  (q5>0)
 *  (d) Lemma 2 majorization g(Kx) < g(x) at x=q/B, exact integers
 *  (e) actual energy increment <= -delta + 4/B (2-Lipschitz on simplex)
 *  (f) the inverse-rule digit description for q5=0 (h positive diffs)
 * Usage: selftest B ntrials seed
 */
#include <stdio.h>
#include "kap.h"

static void lead_K(const i64 *q, i64 B, i64 *k) { /* B*K(q/B), exact */
    i64 s[ND];
    s[0] = B;
    for (int i = 0; i < 5; i++) { s[1 + i] = q[i]; s[6 + i] = B - q[i]; }
    sort_desc(s, ND);
    for (int i = 0; i < 5; i++) k[i] = s[i] - s[ND - 1 - i];
}
static int cmp_desc(const void *a, const void *b) {
    i64 x = *(const i64 *)a, y = *(const i64 *)b; return (x < y) - (x > y);
}
static __int128 energyB(const i64 *G, i64 B) { /* 36 * E_B = sum (6G - B)^2 */
    __int128 e = 0;
    for (int i = 0; i < 6; i++) { __int128 t = 6 * (__int128)G[i] - B; e += t * t; }
    return e;
}

int main(int argc, char **argv) {
    i64 B = atoll(argv[1]); long n = atol(argv[2]); uint64_t seed = strtoull(argv[3], 0, 10);
    rng_t r; rng_seed(&r, seed);
    long bad[8] = {0}, cnt5 = 0, cnt0 = 0;
    i64 maxerr = 0, maxgerr = 0;
    double worst_e = -1e9;
    for (long t = 0; t < n; t++) {
        i64 q[5];
        /* random sorted state; mix uniform and near-degenerate */
        int mode = rng_below(&r, 4);
        for (int i = 0; i < 5; i++) {
            if (mode == 0) q[i] = rng_below(&r, B);
            else if (mode == 1) q[i] = rng_below(&r, 6) * (B / 6) + rng_below(&r, 5) - 2;
            else if (mode == 2) q[i] = (i64)((5 - i) * B / 6) + rng_below(&r, 41) - 20;
            else q[i] = rng_below(&r, 3) == 0 ? 0 : rng_below(&r, B);
            if (q[i] < 0) q[i] = 0;
            if (q[i] > B - 1) q[i] = B - 1;
        }
        qsort(q, 5, sizeof(i64), cmp_desc);
        if (q[0] == 0) continue;
        i64 w[ND] = {0}, o[ND], qn[5], G[6], Gn[6];
        for (int i = 0; i < 5; i++) w[i] = q[i];
        kap_word(w, o, B);
        word_diffs(o, qn);
        gaps(q, B, G); gaps(qn, B, Gn);
        if (q[4] > 0) {
            cnt5++;
            /* (a) digits multiset */
            i64 d[ND], e[ND];
            d[0] = B - 1;
            for (int i = 1; i <= 5; i++) { d[i] = q[i - 1] - (i == 5); d[5 + i] = B - q[i - 1] - (i > 1); }
            memcpy(e, o, sizeof e);
            sort_desc(d, ND); sort_desc(e, ND);
            if (memcmp(d, e, sizeof d)) bad[0]++;
            /* (b) tracking */
            i64 k[5], Gk[6];
            lead_K(q, B, k); gaps(k, B, Gk);
            for (int i = 0; i < 5; i++) { i64 er = llabs(qn[i] - k[i]); if (er > maxerr) maxerr = er; if (er > 1) bad[1]++; }
            for (int i = 0; i < 6; i++) { i64 er = llabs(Gn[i] - Gk[i]); if (er > maxgerr) maxgerr = er; if (er > 2) bad[1]++; }
            /* (c) outer gap */
            i64 a = G[0] + 1 - (q[0] == q[1]);
            i64 pred = a < q[4] ? a : q[4];
            if (Gn[0] != pred) bad[2]++;
            /* (d) majorization of gaps of K vs gaps of x (integers) */
            i64 s1[6], s2[6];
            memcpy(s1, Gk, sizeof s1); memcpy(s2, G, sizeof s2);
            sort_desc(s1, 6); sort_desc(s2, 6);
            i64 c1 = 0, c2 = 0;
            for (int i = 0; i < 6; i++) { c1 += s1[i]; c2 += s2[i]; if (c1 > c2) { bad[3]++; break; } }
            /* (e) energy: E(Fq)-E(q) <= -(E(q)-E(Kq)) + 4B  i.e. E(Fq)-E(Kq) <= 4B */
            __int128 eF = energyB(Gn, B), eK = energyB(Gk, B);
            double ex = (double)(eF - eK) / 36.0 / (double)B;   /* (E_B(Fq)-E_B(Kq))/B */
            if (ex > worst_e) worst_e = ex;
            if (eF - eK > (__int128)36 * 4 * B + 36 * 24) bad[4]++;   /* 4/B plus O(1/B^2) slack */
        } else {
            cnt0++;
            /* (f) zero-innermost digits from the inverse rule */
            int h = 0; while (h < 5 && q[h] > 0) h++;
            i64 d[ND]; int m = 0;
            if (h == 1) { d[m++] = q[0] - 1; d[m++] = B - q[0]; }
            else {
                d[m++] = q[0]; d[m++] = B - q[0];
                for (int i = 2; i < h; i++) { d[m++] = q[i - 1]; d[m++] = B - q[i - 1] - 1; }
                d[m++] = q[h - 1] - 1; d[m++] = B - q[h - 1] - 1;
            }
            while (m < ND) d[m++] = B - 1;
            i64 e[ND]; memcpy(e, o, sizeof e);
            sort_desc(d, ND); sort_desc(e, ND);
            if (memcmp(d, e, sizeof d)) bad[5]++;
        }
    }
    printf("B=%lld trials=%ld q5>0:%ld q5=0:%ld\n", (long long)B, n, cnt5, cnt0);
    printf("digits_mismatch=%ld track_viol=%ld outer_viol=%ld major_viol=%ld energy_viol=%ld inverse_digits_mismatch=%ld\n",
           bad[0], bad[1], bad[2], bad[3], bad[4], bad[5]);
    printf("max |F-BK|=%lld  max gap err=%lld  max (E(Fq)-E(Kq))/B=%.4f\n", (long long)maxerr, (long long)maxgerr, worst_e);
    return 0;
}
