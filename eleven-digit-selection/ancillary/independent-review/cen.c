/* Centered paired-level map with route compression.
 * Written from proof.md's description of the centered map only.
 *
 * Input lines: B G0 G1 G2 G3 G4 [ODD]   (ODD=1: base is B and the state is
 *   computed from the given gaps directly; the caller supplies G+e4 at B+1.)
 * Output: "T <T> MU <mu> LAM <lam> NTOK <n>" then n token lines
 *   "<p> <count> (<delta>) <code_0>,...,<code_{p-1}>": greedy compression of
 *   the 5-bit within-level order codes for steps 0..mu-1 into maximal runs
 *   (>= 4 full copies) of words that translate the actual orbit by a
 *   constant nonzero delta every copy; other steps are single-code tokens
 *   with delta printed as 0. Then "END".
 *   T = mu + 1 is the word clock implied by the difference clock mu.
 *
 * Pairs (left,right) at levels 1..5:
 *   (X1,-X5-6) (X2,-X4-6) (X3,-X3-6) (X4,-X2-6) (X5-6,-X1);
 * bit k = [left >= right], level 1 is the least significant bit.
 * Sorted list e = (-6, hi1, lo1, ..., hi5, lo5); X'_i = e_i - e_{12-i}.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;

static unsigned char step(const ll *X, ll *Y) {
    ll L[5] = {X[0], X[1], X[2], X[3], X[4] - 6};
    ll R[5] = {-X[4] - 6, -X[3] - 6, -X[2] - 6, -X[1] - 6, -X[0]};
    ll hi[5], lo[5];
    unsigned char code = 0;
    for (int k = 0; k < 5; k++) {
        if (L[k] >= R[k]) { code |= (unsigned char)(1u << k); hi[k] = L[k]; lo[k] = R[k]; }
        else { hi[k] = R[k]; lo[k] = L[k]; }
    }
    /* e1=-6, e2=hi1, e3=lo1, e4=hi2, e5=lo2, e6=hi3, e7=lo3, e8=hi4, e9=lo4, e10=hi5, e11=lo5 */
    Y[0] = -6 - lo[4];
    Y[1] = hi[0] - hi[4];
    Y[2] = lo[0] - lo[3];
    Y[3] = hi[1] - hi[3];
    Y[4] = lo[1] - lo[2];
    return code;
}

static ll resid(const ll *X) {
    ll g[6] = {-X[0], X[0] - X[1], X[1] - X[2], X[2] - X[3], X[3] - X[4], X[4]};
    ll m = g[0];
    for (int i = 1; i < 6; i++) if (g[i] < m) m = g[i];
    return -m;
}

/* Check that the word of length p starting at state X is a translation on
 * the actual orbit: X_{p}-X_0 == X_{2p}-X_{p} != 0 (then, the codes being
 * periodic, F(X)=MX+b with M fixing delta, so every copy translates). */
static int orbit_translation(const ll *X, ll p, ll *d) {
    ll A[5], Bv[5], Cv[5], tmp[5];
    memcpy(A, X, sizeof A);
    memcpy(Bv, X, sizeof Bv);
    for (ll t = 0; t < p; t++) { step(Bv, tmp); memcpy(Bv, tmp, sizeof Bv); }
    memcpy(Cv, Bv, sizeof Cv);
    for (ll t = 0; t < p; t++) { step(Cv, tmp); memcpy(Cv, tmp, sizeof Cv); }
    int nz = 0;
    for (int r = 0; r < 5; r++) {
        if (Bv[r] - A[r] != Cv[r] - Bv[r]) return 0;
        d[r] = Bv[r] - A[r];
        if (d[r]) nz = 1;
    }
    return nz;
}

static int eq(const ll *a, const ll *b) { return memcmp(a, b, 5 * sizeof(ll)) == 0; }

int main(void) {
    ll B, G[5];
    size_t cap = 0;
    unsigned char *codes = NULL;
    while (scanf("%lld %lld %lld %lld %lld %lld", &B, &G[0], &G[1], &G[2], &G[3], &G[4]) == 6) {
        ll X[5], Y[5];
        ll s = 0;
        for (int i = 0; i < 5; i++) { s += G[i]; X[i] = (ll)(i + 1) * B - 6 * s; }
        /* X_i = iB - 6 sum_{j<i} G_j  for i = 1..5 */
        size_t need = (size_t)(B / 20) + 100000;
        if (need > cap) { free(codes); codes = malloc(need); cap = need; if (!codes) { fprintf(stderr, "oom\n"); return 1; } }
        ll t = 0;
        ll W[5]; memcpy(W, X, sizeof W);
        while (resid(W) > 504) {
            if ((size_t)t >= cap) { fprintf(stderr, "cap\n"); return 1; }
            codes[t++] = step(W, Y); memcpy(W, Y, sizeof W);
        }
        ll t0 = t;
        /* Brent from W */
        ll tort[5], hare[5], tmp[5];
        ll power = 1, lam = 1;
        memcpy(tort, W, sizeof tort); step(W, hare);
        while (!eq(tort, hare)) {
            if (power == lam) { memcpy(tort, hare, sizeof tort); power *= 2; lam = 0; }
            step(hare, tmp); memcpy(hare, tmp, sizeof hare); lam++;
        }
        memcpy(tort, W, sizeof tort); memcpy(hare, W, sizeof hare);
        for (ll i = 0; i < lam; i++) { step(hare, tmp); memcpy(hare, tmp, sizeof hare); }
        ll mu = 0;
        while (!eq(tort, hare)) {
            if ((size_t)(t0 + mu) >= cap) { fprintf(stderr, "cap2\n"); return 1; }
            codes[t0 + mu] = step(tort, tmp); memcpy(tort, tmp, sizeof tort);
            step(hare, tmp); memcpy(hare, tmp, sizeof hare); mu++;
        }
        if (mu == 0) { printf("FAIL start-cyclic\n"); fflush(stdout); continue; }
        ll MU = t0 + mu;
        /* greedy compression of codes[0..MU) */
        /* first pass count tokens, second pass print; simpler: buffer output */
        size_t outcap = 1 << 20, outlen = 0; char *out = malloc(outcap); ll ntok = 0;
        ll pos = 0;
        ll Xc[5]; memcpy(Xc, X, sizeof Xc);
        while (pos < MU) {
            ll bestp = 1, bestc = 1, bestcov = 1, bd[5] = {0, 0, 0, 0, 0};
            for (ll p = 1; p <= 64 && pos + p <= MU; p++) {
                ll m = 0;
                while (pos + p + m < MU && codes[pos + p + m] == codes[pos + m]) m++;
                ll c = 1 + m / p;
                ll d[5];
                if (c >= 4 && c * p > bestcov && orbit_translation(Xc, p, d)) {
                    bestp = p; bestc = c; bestcov = c * p; memcpy(bd, d, sizeof bd);
                }
            }
            if (outlen + 1024 > outcap) { outcap *= 2; out = realloc(out, outcap); }
            outlen += (size_t)sprintf(out + outlen, "%lld %lld (%lld,%lld,%lld,%lld,%lld) ", bestp, bestc, bd[0], bd[1], bd[2], bd[3], bd[4]);
            for (ll j = 0; j < bestp; j++)
                outlen += (size_t)sprintf(out + outlen, j ? ",%d" : "%d", codes[pos + j]);
            outlen += (size_t)sprintf(out + outlen, "\n");
            ntok++;
            for (ll j = 0; j < bestp * bestc; j++) {
                unsigned char cc = step(Xc, Y);
                if (cc != codes[pos + j]) { fprintf(stderr, "replay mismatch\n"); return 1; }
                memcpy(Xc, Y, sizeof Xc);
            }
            pos += bestp * bestc;
        }
        printf("T %lld MU %lld LAM %lld NTOK %lld\n", MU + 1, MU, lam, ntok);
        fwrite(out, 1, outlen, stdout);
        printf("END\n");
        fflush(stdout);
        free(out);
    }
    return 0;
}
