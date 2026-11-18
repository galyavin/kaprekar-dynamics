/* Exhaustive difference-level checks at fixed bases, written from the paper's definitions.
 *
 * For each base B in [B1,B2] (step S):
 *   - F_B by the hand-derived borrow digits + insertion sort (checked against __int128
 *     digit extraction of L_B(q) on every state when B <= CHECKMAX);
 *   - Lemma A.5: one-step V bound on every positive state; depth-nine strict decrease on
 *     every positive state with q..F^9 q positive and V(q) > 200;
 *   - Lemma A.4 (i) (B>=18), (ii), (iii) on every relevant state;
 *   - first cyclic index h_B on every state; max h, max T_B = h+1; cyclic states and the
 *     max |4u-3B|,|4v-2B|,|4w-B| over cyclic states; same over first cyclic destinations;
 *   - "permanently positive by step six": the last time an orbit has w=0 is <= 5.
 *   - weighted law (if WEIGHTS): exact word counts N(q), mean/var of T_B/B, CDF on a grid.
 * Usage: exhaust B1 B2 [step] [weights 0/1] [lemmas 0/1]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

typedef long long ll;
static int B;

static inline ll IDX(int u, int v, int w) { return (ll)u * (u + 1) * (u + 2) / 6 + (ll)v * (v + 1) / 2 + w; }

static inline void Fstep(int u, int v, int w, int *ou, int *ov, int *ow) {
    int d[7];
    if (w > 0) { d[0]=u; d[1]=v; d[2]=w-1; d[3]=B-1; d[4]=B-w-1; d[5]=B-v-1; d[6]=B-u; }
    else if (v > 0) { d[0]=u; d[1]=v-1; d[2]=B-1; d[3]=B-1; d[4]=B-1; d[5]=B-v-1; d[6]=B-u; }
    else { d[0]=u-1; d[1]=B-1; d[2]=B-1; d[3]=B-1; d[4]=B-1; d[5]=B-1; d[6]=B-u; }
    for (int i = 1; i < 7; i++) { int x = d[i], j = i - 1; while (j >= 0 && d[j] < x) { d[j+1] = d[j]; j--; } d[j+1] = x; }
    *ou = d[0] - d[6]; *ov = d[1] - d[5]; *ow = d[2] - d[4];
}

static void Fint128(int u, int v, int w, int *ou, int *ov, int *ow) {
    __int128 Bp[7]; Bp[0] = 1; for (int i = 1; i < 7; i++) Bp[i] = Bp[i-1] * B;
    __int128 L = (Bp[6] - 1) * u + (Bp[5] - B) * (__int128)v + (Bp[4] - Bp[2]) * (__int128)w;
    int d[7];
    for (int i = 0; i < 7; i++) { d[i] = (int)(L % B); L /= B; }
    for (int i = 1; i < 7; i++) { int x = d[i], j = i - 1; while (j >= 0 && d[j] < x) { d[j+1] = d[j]; j--; } d[j+1] = x; }
    *ou = d[0] - d[6]; *ov = d[1] - d[5]; *ow = d[2] - d[4];
}

static inline int Vf(int u, int v, int w) {
    int g[4] = {B - u, u - v, v - w, w}, s = 0;
    for (int i = 0; i < 4; i++) s += abs(4 * g[i] - B);
    return s;
}

/* exact number of 7-digit words (leading zeros allowed) with difference vector (u,v,w) */
static double fact[8];
static double wordcount(int u, int v, int w) {
    int G[3] = {u - v, v - w, w};
    /* gaps g1..g6 : pairs (g1,g6) sum G0, (g2,g5) sum G1, (g3,g4) sum G2.
       options per pair: list of (zero-pattern, multiplicity) */
    int pat[3][3][2]; double mult[3][3]; int nopt[3];
    for (int p = 0; p < 3; p++) {
        int g = G[p]; nopt[p] = 0;
        if (g == 0) { pat[p][0][0] = 0; pat[p][0][1] = 0; mult[p][0] = 1; nopt[p] = 1; }
        else {
            pat[p][0][0] = 0; pat[p][0][1] = 1; mult[p][0] = 1;
            pat[p][1][0] = 1; pat[p][1][1] = 0; mult[p][1] = 1;
            nopt[p] = 2;
            if (g >= 2) { pat[p][2][0] = 1; pat[p][2][1] = 1; mult[p][2] = g - 1; nopt[p] = 3; }
        }
    }
    double tot = 0;
    for (int a = 0; a < nopt[0]; a++) for (int b = 0; b < nopt[1]; b++) for (int c = 0; c < nopt[2]; c++) {
        int nz[6]; /* nonzero flag of g1..g6 */
        nz[0] = pat[0][a][0]; nz[5] = pat[0][a][1];
        nz[1] = pat[1][b][0]; nz[4] = pat[1][b][1];
        nz[2] = pat[2][c][0]; nz[3] = pat[2][c][1];
        double m = fact[7]; int run = 1;
        for (int i = 0; i < 6; i++) { if (nz[i]) { m /= fact[run]; run = 1; } else run++; }
        m /= fact[run];
        tot += m * mult[0][a] * mult[1][b] * mult[2][c];
    }
    return tot * (B - u);
}

int main(int argc, char **argv) {
    int B1 = atoi(argv[1]), B2 = atoi(argv[2]);
    int ST = argc > 3 ? atoi(argv[3]) : 1;
    int WEIGHTS = argc > 4 ? atoi(argv[4]) : 0;
    int LEMMAS = argc > 5 ? atoi(argv[5]) : 1;
    int CHECKMAX = 200;
    fact[0] = 1; for (int i = 1; i < 8; i++) fact[i] = fact[i-1] * i;
    for (B = B1; B <= B2; B += ST) {
        ll N = (ll)B * (B + 1) * (B + 2) / 6;
        int32_t *h = malloc(sizeof(int32_t) * N);
        for (ll i = 0; i < N; i++) h[i] = -1;
        int8_t *Z = malloc(N); /* last time w==0 along forward orbit, -1 if none, 99 if a cycle has w==0 */
        ll fail_int = 0, fail_one = 0, fail_nine = 0, fail_ai = 0, fail_aii = 0, fail_aiii = 0, nine_cases = 0;
        int max_one = -1000000, max_nine = -1000000;
        /* lemma checks */
        for (int u = 1; u < B; u++) for (int v = 0; v <= u; v++) for (int w = 0; w <= v; w++) {
            int a, b, c;
            Fstep(u, v, w, &a, &b, &c);
            if (B <= CHECKMAX) { int x, y, z; Fint128(u, v, w, &x, &y, &z); if (x != a || y != b || z != c) fail_int++; }
            if (!LEMMAS) continue;
            if (w > 0) {
                int dv = Vf(a, b, c) - Vf(u, v, w);
                if (dv > max_one) max_one = dv;
                if (dv > 8) fail_one++;
                int V0 = Vf(u, v, w);
                if (V0 > 200) {
                    int qu = u, qv = v, qw = w, ok = 1;
                    for (int j = 0; j < 9; j++) { int x, y, z; Fstep(qu, qv, qw, &x, &y, &z); qu = x; qv = y; qw = z; if (qw <= 0) { ok = 0; break; } }
                    if (ok) { nine_cases++; int d9 = Vf(qu, qv, qw) - V0; if (d9 > max_nine) max_nine = d9; if (d9 >= 0) { fail_nine++; if (fail_nine < 5) printf("  NINE FAIL B=%d q=(%d,%d,%d)\n", B, u, v, w); } }
                }
                if (B >= 18 && c > 0) {
                    int x, y, z, x2, y2, z2; Fstep(a, b, c, &x, &y, &z);
                    if (z > 0) { Fstep(x, y, z, &x2, &y2, &z2); if (z2 <= 0) fail_ai++; }
                }
            } else if (v > 0 && B >= 18) {
                int bb = B / 2, beta = B % 2;
                int listed = (u == bb && v == bb) || (beta && ((u == bb + 1 && v == bb) || (u == bb + 1 && v == bb + 1)));
                if (c <= 0) fail_aii++;
                else {
                    int x, y, z; Fstep(a, b, c, &x, &y, &z);
                    if ((z == 0) != listed) fail_aii++;
                    if (z == 0 && !(x == bb + 1 && y == 1)) fail_aii++;
                    if (z > 0) { int x2, y2, z2; Fstep(x, y, z, &x2, &y2, &z2); if (z2 <= 0) fail_aii++; }
                }
            } else if (v == 0) {
                int eu = (B - u > u - 1) ? B - u : u - 1, ev = (B - u > u - 1) ? u - 1 : B - u;
                if (!(a == eu && b == ev && c == 0)) fail_aiii++;
            }
        }
        /* first cyclic index for all states */
        int cap = 1 << 22; int *stk = malloc(sizeof(int) * 3 * cap);
        for (int u0 = 1; u0 < B; u0++) for (int v0 = 0; v0 <= u0; v0++) for (int w0 = 0; w0 <= v0; w0++) {
            ll i0 = IDX(u0, v0, w0);
            if (h[i0] >= 0) continue;
            int n = 0, u = u0, v = v0, w = w0;
            while (1) {
                ll i = IDX(u, v, w);
                if (h[i] >= 0) { /* known */
                    int base = h[i]; int zb = Z[i];
                    for (int k = n - 1; k >= 0; k--) { base++; ll ik = IDX(stk[3*k], stk[3*k+1], stk[3*k+2]); h[ik] = base;
                        if (zb == 99) zb = 99; else if (stk[3*k+2] == 0) zb = (zb < 0) ? 0 : (zb + 1 > 98 ? 98 : zb + 1); else zb = (zb < 0) ? -1 : (zb + 1 > 98 ? 98 : zb + 1);
                        Z[ik] = zb; }
                    break;
                }
                if (h[i] == -2) { /* cycle found: states from position of (u,v,w) on stack */
                    int pos = n - 1;
                    while (!(stk[3*pos] == u && stk[3*pos+1] == v && stk[3*pos+2] == w)) pos--;
                    int cz = 0;
                    for (int k = pos; k < n; k++) if (stk[3*k+2] == 0) cz = 1;
                    for (int k = pos; k < n; k++) { ll ik = IDX(stk[3*k], stk[3*k+1], stk[3*k+2]); h[ik] = 0; Z[ik] = cz ? 99 : -1; }
                    int base = 0; int zb = cz ? 99 : -1;
                    for (int k = pos - 1; k >= 0; k--) { base++; ll ik = IDX(stk[3*k], stk[3*k+1], stk[3*k+2]); h[ik] = base;
                        if (zb == 99) zb = 99; else if (stk[3*k+2] == 0) zb = (zb < 0) ? 0 : (zb + 1 > 98 ? 98 : zb + 1); else zb = (zb < 0) ? -1 : (zb + 1 > 98 ? 98 : zb + 1);
                        Z[ik] = zb; }
                    break;
                }
                h[i] = -2;
                if (n >= cap) { fprintf(stderr, "stack overflow\n"); exit(1); }
                stk[3*n] = u; stk[3*n+1] = v; stk[3*n+2] = w; n++;
                int a, b, c; Fstep(u, v, w, &a, &b, &c); u = a; v = b; w = c;
            }
        }
        /* statistics */
        int maxh = 0; ll ncyc = 0; int maxXc = 0; int mu = 0, mv = 0, mw = 0;
        int last_zero_max = 0; /* max over states of last time with w==0 along orbit up to h+cycle */
        for (int u = 1; u < B; u++) for (int v = 0; v <= u; v++) for (int w = 0; w <= v; w++) {
            int hh = h[IDX(u, v, w)];
            if (hh > maxh) { maxh = hh; mu = u; mv = v; mw = w; }
            if (hh == 0) {
                ncyc++;
                int X1 = abs(4 * u - 3 * B), X2 = abs(4 * v - 2 * B), X3 = abs(4 * w - B);
                int m = X1 > X2 ? X1 : X2; m = m > X3 ? m : X3;
                if (m > maxXc) maxXc = m;
                
            }
        }
        for (ll i = 0; i < N; i++) { /* index 0..B-1 with u=0 are invalid states: skip by checking via loop below */ }
        for (int u = 1; u < B; u++) for (int v = 0; v <= u; v++) for (int w = 0; w <= v; w++) {
            int zz = Z[IDX(u, v, w)]; if (zz > last_zero_max) last_zero_max = zz; }
        printf("B=%d states=%lld int128mismatch=%lld | A5 one-step fails=%lld (max dV=%d) | A5 nine-step cases=%lld fails=%lld (max dV9=%d) | A4 i/ii/iii fails=%lld/%lld/%lld | last zero time=%d\n",
               B, N - 1, fail_int, fail_one, max_one, nine_cases, fail_nine, max_nine, fail_ai, fail_aii, fail_aiii, last_zero_max);
        printf("B=%d max h=%d at (%d,%d,%d)  max T=%d  T/B=%.5f  cyclic states=%lld  max|X| on cycles=%d (bound 136) => max dist*B=%.2f (bound 34)\n",
               B, maxh, mu, mv, mw, maxh + 1, (maxh + 1.0) / B, ncyc, maxXc, maxXc / 4.0);
        if (WEIGHTS) {
            /* exact weighted law of T_B/B over all B^7-B words */
            double tot = 0, s1 = 0, s2 = 0; int NG = 26; double cdf[26]; memset(cdf, 0, sizeof cdf);
            double grid[26]; for (int k = 0; k < NG; k++) grid[k] = 0.01 * k;
            double mexp1 = 0, mexp2 = 0; /* E exp(theta T/B) for theta = 4, -4 */
            for (int u = 1; u < B; u++) for (int v = 0; v <= u; v++) for (int w = 0; w <= v; w++) {
                double n = wordcount(u, v, w);
                int hh = h[IDX(u, v, w)];
                double T = hh + 1;
                double t = T / B;
                if (hh == 0) { /* one cyclic word with T=0, the rest T=1 */
                    tot += 1;
                    for (int k = 0; k < NG; k++) if (0 <= grid[k]) cdf[k] += 1;
                    mexp1 += 1; mexp2 += 1;
                    n -= 1;
                }
                tot += n; s1 += n * t; s2 += n * t * t;
                for (int k = 0; k < NG; k++) if (t <= grid[k]) cdf[k] += n;
                mexp1 += n * exp(4 * t); mexp2 += n * exp(-4 * t);
            }
            double expected = pow((double)B, 7) - B;
            printf("B=%d weights total=%.6e expected=%.6e rel=%.2e  mean(T/B)=%.6f var=%.6f  E e^{4T/B}=%.6f E e^{-4T/B}=%.6f\n",
                   B, tot, expected, tot / expected - 1, s1 / tot, s2 / tot - (s1 / tot) * (s1 / tot), mexp1 / tot, mexp2 / tot);
            printf("B=%d CDF", B);
            for (int k = 0; k < NG; k++) printf(" %.5f", cdf[k] / tot);
            printf("\n");
        }
        fflush(stdout);
        free(h); free(stk); free(Z);
    }
    return 0;
}
