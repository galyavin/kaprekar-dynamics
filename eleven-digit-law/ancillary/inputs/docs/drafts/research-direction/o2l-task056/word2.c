/* Faster whole-word replay with the same output as word.c ("T period" per line).
   One pass to N0 = 0.135 B storing a checkpoint every 512 steps; the period is found by
   stepping from the state at N0; the first cyclic time by binary search over checkpoints
   and a final two-pointer scan inside one block. Falls back to Brent when the state at N0
   is not cyclic within 100000 steps. Ground truth is word.c; validate on samples. */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int64_t d[11]; } wd;
static int64_t BASE;

static inline wd step(wd w) {
    int64_t s[11];
    memcpy(s, w.d, sizeof s);
    for (int i = 1; i < 11; i++) { int64_t v = s[i]; int j = i - 1;
        while (j >= 0 && s[j] < v) { s[j+1] = s[j]; j--; } s[j+1] = v; }
    wd r; int64_t borrow = 0;
    for (int i = 10; i >= 0; i--) {
        int64_t v = s[i] - s[10 - i] - borrow;
        if (v < 0) { v += BASE; borrow = 1; } else borrow = 0;
        r.d[i] = v;
    }
    return r;
}
static inline int eq(wd a, wd b) { return memcmp(a.d, b.d, sizeof a.d) == 0; }

static void brent(wd w0, long long *mu_out, long long *lam_out) {
    long long power = 1, lam = 1;
    wd tort = w0, hare = step(w0);
    while (!eq(tort, hare)) {
        if (power == lam) { tort = hare; power *= 2; lam = 0; }
        hare = step(hare); lam++;
    }
    wd x = w0, y = w0;
    for (long long i = 0; i < lam; i++) y = step(y);
    long long mu = 0;
    while (!eq(x, y)) { x = step(x); y = step(y); mu++; }
    *mu_out = mu; *lam_out = lam;
}

static int cyclic_at(wd w, long long lam) {   /* is w periodic with period lam? */
    wd y = w;
    for (long long i = 0; i < lam; i++) y = step(y);
    return eq(w, y);
}

int main(void) {
    long long B, a, b, c, d, e;
    const int BLK = 512;
    while (scanf("%lld %lld %lld %lld %lld %lld", &B, &a, &b, &c, &d, &e) == 6) {
        BASE = B;
        wd w0; memset(&w0, 0, sizeof w0);
        w0.d[0] = a; w0.d[1] = b; w0.d[2] = c; w0.d[3] = d; w0.d[4] = e;
        long long N0 = (long long)(0.135 * (double)B);
        long long ncp = N0 / BLK + 1;
        wd *cp = malloc(sizeof(wd) * (size_t)ncp);
        wd w = w0;
        for (long long t = 0; t < N0; t++) { if (t % BLK == 0) cp[t / BLK] = w; w = step(w); }
        /* period from the state at N0 */
        wd y = step(w); long long lam = 1;
        while (!eq(y, w) && lam < 100000) { y = step(y); lam++; }
        long long mu;
        if (!eq(y, w)) { brent(w0, &mu, &lam); printf("%lld %lld\n", mu, lam); free(cp); continue; }
        /* binary search for the first cyclic checkpoint index */
        long long lo = 0, hi = (N0 - 1) / BLK;        /* checkpoints 0..hi exist */
        if (!cyclic_at(cp[hi], lam)) { lo = hi; }      /* first cyclic time lies after cp[hi] */
        else if (cyclic_at(cp[0], lam)) { lo = -1; }
        else { while (hi - lo > 1) { long long mid = (lo + hi) / 2;
                   if (cyclic_at(cp[mid], lam)) hi = mid; else lo = mid; } }
        if (lo < 0) mu = 0;
        else {
            wd x = cp[lo], z = cp[lo];
            for (long long i = 0; i < lam; i++) z = step(z);
            long long t = lo * BLK;
            while (!eq(x, z)) { x = step(x); z = step(z); t++; }
            mu = t;
        }
        printf("%lld %lld\n", mu, lam);
        free(cp);
    }
    return 0;
}
