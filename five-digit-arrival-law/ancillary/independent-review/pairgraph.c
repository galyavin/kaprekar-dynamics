/* Independent exhaustive pair-graph analysis for five-digit Kaprekar dynamics.
 * Independent review code, written from the paper's text only.
 *
 * The pair map is computed from the DEFINITION: N = desc - asc of the sorted
 * digits (d1..d5) equals (B^4-1)u + (B^3-B)v; we form that integer, extract
 * its five base-B digits, sort them, and read off (u', v') = (s1-s5, s2-s4).
 * It does not use the borrowed-digit formula (17) or the paper's branches.
 *
 * Input weights: number of ordered digit strings with sorted pair (u,v).  By
 * default the closed form (18) is used; `strings_exhaustive.c` checks (18)
 * against complete enumeration of all B^5 strings for small B.
 *
 * usage: pairgraph B [checks]
 *   prints exact distribution of T_B (ordered-input weights), cycle labels,
 *   maximum, witness (37), and finite-base tests of (3),(9),(22)-(23),(36),
 *   the six-step entry lemma and the rare-basin formula (7).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef unsigned __int128 u128;
static long B;
static long npairs;

static inline long idx(long u, long v) { return (u - 1) * (u + 2) / 2 + v; }
static long *U_of, *V_of;

static void pair_step(long u, long v, long *uo, long *vo) {
    u128 Bb = (u128)B;
    u128 N = (Bb * Bb * Bb * Bb - 1) * (u128)u + (Bb * Bb * Bb - Bb) * (u128)v;
    long d[5];
    for (int k = 0; k < 5; k++) { d[k] = (long)(N % Bb); N /= Bb; }
    if (N != 0) { fprintf(stderr, "overflow digits\n"); exit(2); }
    /* insertion sort descending */
    for (int i = 1; i < 5; i++) {
        long x = d[i]; int j = i - 1;
        while (j >= 0 && d[j] < x) { d[j + 1] = d[j]; j--; }
        d[j + 1] = x;
    }
    *uo = d[0] - d[4];
    *vo = d[1] - d[3];
}

static u128 weight(long u, long v) {
    /* (18) */
    if (v == 0) return (u128)(B - u) * (u128)(20 * u - 10);
    if (u == v) return (u128)(B - u) * (u128)(30 * u - 10);
    return (u128)(B - u) * (u128)(120 * (u128)v * (u128)(u - v) - 20);
}

static int inH(long u, long v) {
    long a = u + v - B, b = 2 * u - v - B;
    return 2 * v <= B - 1 && 2 * (u - v) <= B - 3 && a >= -2 && b >= -2;
}
static int inJ(long u, long v) {
    return 2 * u >= B - 1 && 2 * v <= B + 1 && 2 * (u - v) <= B + 2;
}
static int inT2(long u, long v) {
    if (B % 2) return 0;
    long m = B / 2;
    return (u == m && v == m - 1) || (u == m + 1 && v == 0);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: pairgraph B\n"); return 1; }
    B = atol(argv[1]);
    npairs = (B - 1) * (B + 2) / 2;
    uint32_t *nxt = malloc(sizeof(uint32_t) * npairs);
    int32_t *h = malloc(sizeof(int32_t) * npairs);     /* pair-cycle arrival */
    int32_t *cyc = malloc(sizeof(int32_t) * npairs);   /* terminal cycle id */
    int32_t *pathpos = malloc(sizeof(int32_t) * npairs);
    uint32_t *stack = malloc(sizeof(uint32_t) * (npairs < 4000000 ? npairs : 4000000));
    U_of = malloc(sizeof(long) * npairs);
    V_of = malloc(sizeof(long) * npairs);
    if (!nxt || !h || !cyc || !pathpos || !stack || !U_of || !V_of) { fprintf(stderr, "oom\n"); return 2; }
    for (long u = 1; u < B; u++)
        for (long v = 0; v <= u; v++) {
            long i = idx(u, v), u2, v2;
            U_of[i] = u; V_of[i] = v;
            pair_step(u, v, &u2, &v2);
            if (u2 < 1 || u2 >= B || v2 < 0 || v2 > u2) { fprintf(stderr, "bad image\n"); return 3; }
            nxt[i] = (uint32_t)idx(u2, v2);
        }
    /* arrival times via path walking */
    for (long i = 0; i < npairs; i++) { h[i] = -1; pathpos[i] = -1; }
    int ncyc = 0;
    long cyc_rep[64]; int cyc_len[64];
    for (long s = 0; s < npairs; s++) {
        if (h[s] >= 0) continue;
        long top = 0, x = s;
        while (h[x] < 0 && pathpos[x] < 0) { pathpos[x] = (int32_t)top; stack[top++] = (uint32_t)x; x = nxt[x]; }
        if (h[x] < 0) { /* new cycle found starting at pathpos[x] */
            long p0 = pathpos[x];
            if (ncyc >= 64) { fprintf(stderr, "too many cycles\n"); return 4; }
            long rep = x; int len = 0;
            for (long k = p0; k < top; k++) { long y = stack[k]; h[y] = 0; cyc[y] = ncyc; if ((long)y < rep) rep = y; len++; }
            cyc_rep[ncyc] = rep; cyc_len[ncyc] = len; ncyc++;
            top = p0;
        }
        for (long k = top - 1; k >= 0; k--) { long y = stack[k]; h[y] = h[nxt[y]] + 1; cyc[y] = cyc[nxt[y]]; }
        for (long k = 0; k < top; k++) pathpos[stack[k]] = -1;
        /* clear also cycle nodes' pathpos */
        if (pathpos[x] >= 0) {}
    }
    for (long i = 0; i < npairs; i++) pathpos[i] = -1;

    /* label cycles: by (a,b) membership in the paper's lists, or T2 */
    char label[64][8];
    for (int c = 0; c < ncyc; c++) {
        long r = cyc_rep[c]; long u = U_of[r], v = V_of[r];
        long a = u + v - B, b = 2 * u - v - B;
        strcpy(label[c], "OTHER");
        /* walk the cycle and test membership */
        long y = r; int isT1 = 1, isT3 = 1, isT4 = 1, isT5 = 1, isT2 = 1;
        static const int T1s[1][2] = {{0,0}};
        static const int T3s[3][2] = {{-2,0},{1,3},{3,-2}};
        static const int T4s[4][2] = {{-1,-1},{-1,2},{2,2},{2,-1}};
        static const int T5s[4][2] = {{-1,0},{0,2},{2,0},{0,-1}};
        for (int k = 0; k < cyc_len[c]; k++) {
            long uu = U_of[y], vv = V_of[y]; long aa = uu + vv - B, bb = 2 * uu - vv - B;
            int f;
            f = 0; for (int t = 0; t < 1; t++) f |= (aa == T1s[t][0] && bb == T1s[t][1]); isT1 &= f;
            f = 0; for (int t = 0; t < 3; t++) f |= (aa == T3s[t][0] && bb == T3s[t][1]); isT3 &= f;
            f = 0; for (int t = 0; t < 4; t++) f |= (aa == T4s[t][0] && bb == T4s[t][1]); isT4 &= f;
            f = 0; for (int t = 0; t < 4; t++) f |= (aa == T5s[t][0] && bb == T5s[t][1]); isT5 &= f;
            isT2 &= inT2(uu, vv);
            y = nxt[y];
        }
        if (isT1 && cyc_len[c] == 1) strcpy(label[c], "T1");
        if (isT3 && cyc_len[c] == 3) strcpy(label[c], "T3");
        if (isT4 && cyc_len[c] == 4) strcpy(label[c], "T4");
        if (isT5 && cyc_len[c] == 4) strcpy(label[c], "T5");
        if (isT2 && cyc_len[c] == 2) strcpy(label[c], "T2");
        (void)a; (void)b;
    }

    /* exact distribution of T over ordered inputs */
    long maxh = 0; for (long i = 0; i < npairs; i++) if (h[i] > maxh) maxh = h[i];
    long maxT = maxh + 1;
    u128 *hist = calloc(maxT + 2, sizeof(u128));
    u128 *histc = calloc((size_t)(maxT + 2) * ncyc, sizeof(u128));
    u128 tot = 0; long ncycpairs = 0;
    for (long i = 0; i < npairs; i++) {
        u128 w = weight(U_of[i], V_of[i]);
        tot += w;
        hist[h[i] + 1] += w;
        histc[(size_t)(h[i] + 1) * ncyc + cyc[i]] += w;
        if (h[i] == 0) ncycpairs++;
    }
    /* numerical cycle inputs: K(q) for cyclic q has pair f(q); time 0 */
    for (long i = 0; i < npairs; i++) if (h[i] == 0) {
        long j = nxt[i];
        hist[1] -= 1; hist[0] += 1;
        histc[(size_t)1 * ncyc + cyc[j]] -= 1; histc[(size_t)0 * ncyc + cyc[j]] += 1;
    }
    u128 expect_tot = (u128)B * B * B * B * B - (u128)B;
    printf("B %ld pairs %ld total_weight_ok %d cycles %d cyclic_pairs %ld\n", B, npairs, tot == expect_tot, ncyc, ncycpairs);
    for (int c = 0; c < ncyc; c++) {
        long r = cyc_rep[c]; u128 m = 0;
        for (long t = 0; t <= maxT; t++) m += histc[(size_t)t * ncyc + c];
        printf("cycle %d label %s len %d rep (%ld,%ld) (a,b)=(%ld,%ld) mass %.12f count %.0f\n", c, label[c], cyc_len[c],
               U_of[r], V_of[r], U_of[r] + V_of[r] - B, 2 * U_of[r] - V_of[r] - B, (double)m / (double)tot, (double)m);
    }
    /* CDF at requested t */
    double ts[] = {1.0/25, 0.1, 0.2, 0.3, 1.0/3};
    for (int k = 0; k < 5; k++) {
        double t = ts[k]; u128 acc = 0;
        for (long T = 0; T <= maxT; T++) if ((double)T <= t * (double)B + 1e-9) acc += hist[T];
        printf("cdf t=%.6f P(T<=tB)=%.10f", t, (double)acc / (double)tot);
        for (int c = 0; c < ncyc; c++) {
            u128 a2 = 0, m = 0;
            for (long T = 0; T <= maxT; T++) { m += histc[(size_t)T * ncyc + c]; if ((double)T <= t * (double)B + 1e-9) a2 += histc[(size_t)T * ncyc + c]; }
            if ((double)m / (double)tot > 1e-6) printf(" | %s joint %.10f cond %.10f", label[c], (double)a2 / (double)tot, (double)a2 / (double)m);
        }
        printf("\n");
    }
    long double m1 = 0, m2 = 0;
    for (long T = 0; T <= maxT; T++) { long double x = (long double)T / B; m1 += x * (long double)hist[T]; m2 += x * x * (long double)hist[T]; }
    printf("moments E[T/B]=%.10Lf E[(T/B)^2]=%.10Lf\n", m1 / (long double)tot, m2 / (long double)tot);
    /* maximum */
    long arg = 0; for (long i = 0; i < npairs; i++) if (h[i] > h[arg]) arg = i;
    double sq = sqrt((double)B), ce = ceil(sq);
    printf("maxT %ld at pair (%ld,%ld)  maxT-B/3-sqrtB/3 = %.4f  bound3 B/3+ceil(sqrtB)+20 = %.2f ok3 %d  lower9 %.3f upper9 %.3f ok9 %d\n",
           maxT, U_of[arg], V_of[arg], maxT - B / 3.0 - sq / 3.0, B / 3.0 + ce + 20, maxT <= B / 3.0 + ce + 20,
           B / 3.0 + sq / 3.0 - 5, B / 3.0 + sq / 3.0 + 60, (maxT >= B / 3.0 + sq / 3.0 - 5) && (maxT <= B / 3.0 + sq / 3.0 + 60));
    /* witness (37) */
    {
        long vs = (B - 1) / 2, n = 1;
        while ((n + 1) * n / 2 <= B - 3 - vs) n++;
        long u0 = B - 1, v0 = vs + n * (n - 1) / 2;
        long T = h[idx(u0, v0)] + 1;
        printf("witness37 n %ld pair (%ld,%ld) T %ld  T-B/3-sqrtB/3 = %.4f  lower_claim(B/3+sqrtB/3-4) ok %d  maxT-T %ld\n",
               n, u0, v0, T, T - B / 3.0 - sq / 3.0, T >= B / 3.0 + sq / 3.0 - 4, maxT - T);
    }
    if (argc >= 3) {
        /* finite-base tests of the reductions */
        int ok23 = 1; long lo3 = 1000, hi3 = -1000;   /* 3h - 2M */
        long maxsix = 0; long maxsteps = 0; int never = 0;
        long ex36 = -1000000; int ok36 = 1;
        uint8_t *dist = malloc(npairs);
        for (long i = 0; i < npairs; i++) dist[i] = 255;
        for (long s = 0; s < npairs; s++) {
            long x = s, k = 0;
            while (!(inH(U_of[x], V_of[x]) || inT2(U_of[x], V_of[x])) && k < 100000) { x = nxt[x]; k++; }
            if (k >= 100000) never = 1;
            if (k > maxsteps) maxsteps = k;
            if (inJ(U_of[s], V_of[s]) && k > maxsix) maxsix = k;
            if (inH(U_of[s], V_of[s])) {
                long u = U_of[s], v = V_of[s], a = u + v - B, b = 2 * u - v - B, M = a > b ? a : b;
                long d = 3 * h[s] - 2 * M;
                if (d < lo3) lo3 = d; if (d > hi3) hi3 = d;
                if (d < -9 || d > 18) ok23 = 0;
            }
            if (inJ(U_of[s], V_of[s]) && 4 * U_of[s] >= 3 * B) {
                long y = s, j = 0;
                while (!inH(U_of[y], V_of[y]) && j <= 6) { y = nxt[y]; j++; }
                if (j > 6) ok36 = 0;
                else {
                    long u = U_of[y], v = V_of[y], a = u + v - B, b = 2 * u - v - B, M = a > b ? a : b;
                    long e = 2 * M - (2 * U_of[s] - B);
                    if (e > ex36) ex36 = e;
                    if (e > 40) ok36 = 0;
                }
            }
        }
        printf("check23 3h-2M in [%ld,%ld] (claim [-9,18]) ok %d\n", lo3, hi3, ok23);
        printf("checksix max steps J_B->H_B u T2 = %ld (claim <=6) ok %d\n", maxsix, maxsix <= 6);
        printf("checkglobal max steps ->H_B u T2 = %ld  ceil(sqrtB)+13 = %.0f ok %d never %d\n", maxsteps, ce + 13, maxsteps <= ce + 13, never);
        printf("check36 max of 2max(a_j,b_j)-(2u0-B) = %ld (claim <=40) ok %d\n", ex36, ok36);
        /* rare basin (7) */
        if (B % 2 == 0) {
            int t2 = -1; for (int c = 0; c < ncyc; c++) if (!strcmp(label[c], "T2")) t2 = c;
            long basin = 0; for (long i = 0; i < npairs; i++) if (t2 >= 0 && cyc[i] == t2) basin++;
            double N0 = 0, N1 = 0, N2 = 0, Nother = 0;
            if (t2 >= 0) for (long T = 0; T <= maxT; T++) {
                double val = (double)histc[(size_t)T * ncyc + t2];
                if (T == 0) N0 = val; else if (T == 1) N1 = val; else if (T == 2) N2 = val; else Nother += val;
            }
            double b2 = (double)B * B;
            printf("rare T2 basin_pairs %ld N0 %.0f N1 %.0f N2 %.0f Nother %.0f | formula N0 2 N1 %.0f N2 %.0f | match %d\n",
                   basin, N0, N1, N2, Nother, 35 * b2 - 75 * B - 12, 5 * b2 - 5 * B,
                   basin == 3 && N0 == 2 && N1 == 35 * b2 - 75 * B - 12 && N2 == 5 * b2 - 5 * B && Nother == 0);
        }
    }
    return 0;
}
