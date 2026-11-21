/* Index-2 primes: strictness and W_p != {0} straight from P10's row definition.
 *
 * Index 2 means <-1,2> has index 2 in (Z/p)^x, i.e. <-1,2> = squares; we test it as
 * p = 1 mod 8 and ord_p(2) = (p-1)/2 (asserted equivalent in the report).
 * Then J = {0} U {slots j in 1..k with (j/p) = +1}; there is one orbit outside J (the non-residue
 * slots), so W_p is spanned by the indicator profile of the non-residue slots (minus its mean),
 * and with v = 1_{NR slots}:
 *     D_i = d_i(v) = #NR in slots [m-i, m-1]  -  #NR in slots [0, i-1],   2 <= i <= p0 = floor(m/2).
 * (raw suffix/prefix definition; no use of eq. (prefix) or of the F-identity).
 *   W_p != 0  <=>  all D_i >= 0  or  all D_i <= 0
 *   strict    <=>  all D_i >  0  or  all D_i <  0
 * Usage: idx2scan N mode   (mode 0: strict only, early exit; mode 1: also W != 0, full scan)
 * Output: per-bin counts and the list of strict primes; the largest row index at which a
 * non-strict prime was first decided.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int jacobi(uint64_t a, uint64_t n) { /* n odd positive */
    int t = 1;
    a %= n;
    while (a) {
        while ((a & 1) == 0) {
            a >>= 1;
            uint64_t r = n & 7;
            if (r == 3 || r == 5) t = -t;
        }
        uint64_t tmp = a; a = n; n = tmp;
        if ((a & 3) == 3 && (n & 3) == 3) t = -t;
        a %= n;
    }
    return n == 1 ? t : 0;
}
static uint64_t pw(uint64_t b, uint64_t e, uint64_t m) {
    unsigned __int128 r = 1, x = b % m;
    while (e) { if (e & 1) r = r * x % m; x = x * x % m; e >>= 1; }
    return (uint64_t)r;
}
static uint32_t *plist; static int np;

int main(int argc, char **argv) {
    uint64_t N = strtoull(argv[1], 0, 10);
    int mode = atoi(argv[2]);
    uint8_t *comp = calloc(N + 1, 1);
    for (uint64_t i = 2; i * i <= N; i++) if (!comp[i]) for (uint64_t j = i * i; j <= N; j += i) comp[j] = 1;
    uint64_t sq = 1; while (sq * sq <= N) sq++;
    plist = malloc(sizeof(uint32_t) * (sq + 10)); np = 0;
    for (uint64_t i = 2; i <= sq; i++) if (!comp[i]) plist[np++] = (uint32_t)i;
    /* bins */
    double edges[] = {0, 6000, 30000, 150000, 300000, 600000, 1000000, 2000000, 3000000, 10000000, 30000000, 100000000, 1e18};
    int nb = sizeof(edges) / sizeof(edges[0]) - 1;
    long cnt[16] = {0}, cW[16] = {0}, cS[16] = {0};
    long total = 0, maxdecide = 0; uint64_t maxdecide_p = 0;
    long below[6] = {0}; uint64_t lim[6] = {6000, 20000, 30000, 200000, 1000000, 30000000};
    long Wbelow[6] = {0};
    for (uint64_t p = 17; p <= N; p += 8) {
        if (comp[p]) continue;
        uint64_t h = (p - 1) / 2;
        if (pw(2, h, p) != 1) { fprintf(stderr, "2 not QR at %llu\n", (unsigned long long)p); return 1; }
        /* ord_p(2) == h ? */
        uint64_t x = h; int ok = 1;
        for (int t = 0; t < np && (uint64_t)plist[t] * plist[t] <= x; t++) {
            uint64_t q = plist[t];
            if (x % q == 0) { if (pw(2, h / q, p) == 1) { ok = 0; break; } while (x % q == 0) x /= q; }
        }
        if (ok && x > 1 && pw(2, h / x, p) == 1) ok = 0;
        if (!ok) continue;
        total++;
        int b = 0; while (!(p >= edges[b] && p < edges[b + 1])) b++;
        cnt[b]++;
        for (int u = 0; u < 6; u++) if (p < lim[u]) below[u]++;
        uint64_t k = h, m = k + 1, p0 = m / 2;
        long D = 0; int allpos = 1, allneg = 1, sp = 1, sn = 1; long decided = 0;
        for (uint64_t i = 1; i <= p0; i++) {
            uint64_t sslot = m - i, pslot = i - 1;
            if (sslot >= 1 && jacobi(sslot, p) == -1) D++;
            if (pslot >= 1 && jacobi(pslot, p) == -1) D--;
            if (i >= 2) {
                if (D < 0) allpos = 0, sp = 0;
                if (D > 0) allneg = 0, sn = 0;
                if (D == 0) sp = 0, sn = 0;
                if (!sp && !sn && decided == 0) decided = (long)i;
                if (mode == 0 && !sp && !sn) break;
                if (mode == 1 && !allpos && !allneg) break;
            }
        }
        if (decided > maxdecide) { maxdecide = decided; maxdecide_p = p; }
        int W = (mode == 1) ? (allpos || allneg) : -1;
        if (W == 1) { cW[b]++; for (int u = 0; u < 6; u++) if (p < lim[u]) Wbelow[u]++; }
        if (sp || sn) { cS[b]++; printf("STRICT %llu (orientation %s)\n", (unsigned long long)p, sp ? "D>0" : "D<0"); }
        if (mode == 1 && W == 1 && allpos) printf("W!=0 with D>=0 orientation at %llu\n", (unsigned long long)p);
    }
    printf("N=%llu mode=%d  index-2 primes 17<=p<=N: %ld\n", (unsigned long long)N, mode, total);
    for (int u = 0; u < 6; u++) if (lim[u] <= N + 1) printf("  below %llu: %ld index-2 primes%s", (unsigned long long)lim[u], below[u], mode ? "" : "\n"),
        mode ? printf(", W!=0: %ld\n", Wbelow[u]) : 0;
    for (int b = 0; b < nb; b++) if (cnt[b]) {
        printf("  bin [%.0f, %.0f): idx2=%ld strict=%ld", edges[b], edges[b + 1] > 1e17 ? (double)N : edges[b + 1], cnt[b], cS[b]);
        if (mode == 1) printf(" W!=0=%ld share=%.4f", cW[b], (double)cW[b] / cnt[b]);
        printf("\n");
    }
    printf("  largest row index deciding non-strictness: %ld (at p=%llu)\n", maxdecide, (unsigned long long)maxdecide_p);
    return 0;
}
