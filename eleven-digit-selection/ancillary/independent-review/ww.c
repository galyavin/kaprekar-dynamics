/* Independent whole-word eleven-digit Kaprekar clock.
 *
 * Written from the digit-subtraction definitions.
 *
 * Input lines:  B G0 G1 G2 G3 G4      (G5 = B - sum, must be >= 0, G0 >= 1)
 * Output lines: B G0 G1 G2 G3 G4 T MUD LAM T0
 *   T   = first t >= 0 such that the word w_t lies on a cycle of the
 *         whole-word map w -> desc(w) - asc(w) (11 base-B digits, leading
 *         zeros kept), starting from the canonical word with these gaps;
 *   MUD = first t such that the five-difference state of w_t is cyclic;
 *   LAM = least period of the terminal word cycle;
 *   T0  = offset where cycle search started (diagnostic only).
 *
 * Each update sorts the 11 digits and performs schoolbook base-B
 * subtraction with borrows.  No centered model or pairing lemma is used
 * for correctness: Brent's algorithm is run from an offset t0, and the
 * result is accepted only if the word at t0 is not itself cyclic (so that
 * T = t0 + mu_rel exactly); otherwise the search is redone from t0 = 0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;
#define N 11

static ll BASE;

static void sort_desc(const ll *w, ll *s) {
    for (int i = 0; i < N; i++) s[i] = w[i];
    for (int i = 1; i < N; i++) {
        ll v = s[i]; int j = i - 1;
        while (j >= 0 && s[j] < v) { s[j + 1] = s[j]; j--; }
        s[j + 1] = v;
    }
}

/* w[p] = digit at position p (p = 0 least significant). */
static void kstep(const ll *w, ll *out) {
    ll s[N];
    sort_desc(w, s);              /* s[0] largest digit */
    /* D has digit s[10-p] at position p; A has digit s[p] at position p. */
    ll borrow = 0;
    for (int p = 0; p < N; p++) {
        ll d = s[N - 1 - p] - s[p] - borrow;
        if (d < 0) { d += BASE; borrow = 1; } else borrow = 0;
        out[p] = d;
    }
    if (borrow) { fprintf(stderr, "negative difference\n"); exit(2); }
}

static void diffs(const ll *w, ll *q) {
    ll s[N];
    sort_desc(w, s);
    for (int i = 0; i < 5; i++) q[i] = s[i] - s[N - 1 - i];
}

static int weq(const ll *a, const ll *b) { return memcmp(a, b, sizeof(ll) * N) == 0; }

static ll residual(const ll *w) {
    ll q[5]; diffs(w, q);
    ll g[6];
    g[0] = BASE - q[0];
    for (int i = 0; i < 4; i++) g[i + 1] = q[i] - q[i + 1];
    g[5] = q[4];
    ll R = 0;
    for (int i = 0; i < 6; i++) { ll h = BASE - 6 * g[i]; if (h > R) R = h; }
    return R;
}

/* Brent from x0; returns lam, mu (relative). */
static void brent(const ll *x0, ll *lam_out, ll *mu_out, ll *prev_out, int *has_prev) {
    ll tort[N], hare[N], tmp[N];
    ll power = 1, lam = 1;
    memcpy(tort, x0, sizeof tort);
    kstep(x0, hare);
    while (!weq(tort, hare)) {
        if (power == lam) { memcpy(tort, hare, sizeof tort); power *= 2; lam = 0; }
        kstep(hare, tmp); memcpy(hare, tmp, sizeof hare);
        lam++;
    }
    memcpy(tort, x0, sizeof tort);
    memcpy(hare, x0, sizeof hare);
    for (ll i = 0; i < lam; i++) { kstep(hare, tmp); memcpy(hare, tmp, sizeof hare); }
    ll mu = 0;
    *has_prev = 0;
    while (!weq(tort, hare)) {
        memcpy(prev_out, tort, sizeof(ll) * N); *has_prev = 1;
        kstep(tort, tmp); memcpy(tort, tmp, sizeof tort);
        kstep(hare, tmp); memcpy(hare, tmp, sizeof hare);
        mu++;
    }
    *lam_out = lam; *mu_out = mu;
}

int main(void) {
    ll B, G[6];
    while (scanf("%lld %lld %lld %lld %lld %lld", &B, &G[0], &G[1], &G[2], &G[3], &G[4]) == 6) {
        BASE = B;
        G[5] = B - (G[0] + G[1] + G[2] + G[3] + G[4]);
        if (G[0] < 1 || G[5] < 0) { printf("%lld %lld %lld %lld %lld %lld BAD\n", B, G[0], G[1], G[2], G[3], G[4]); continue; }
        if (G[1] < 0 || G[2] < 0 || G[3] < 0 || G[4] < 0) { printf("%lld BAD\n", B); continue; }
        ll q[5];
        for (int i = 0; i < 5; i++) { q[i] = 0; for (int j = i + 1; j < 6; j++) q[i] += G[j]; }
        /* canonical word: digits q1..q5 and six zeros (sorted pairs give q) */
        ll w0[N];
        for (int p = 0; p < N; p++) w0[p] = 0;
        for (int i = 0; i < 5; i++) w0[N - 1 - i] = q[i];
        {   /* sanity: diffs(w0) == q */
            ll qq[5]; diffs(w0, qq);
            for (int i = 0; i < 5; i++) if (qq[i] != q[i]) { fprintf(stderr, "canon\n"); exit(3); }
        }
        ll w[N], tmp[N];
        memcpy(w, w0, sizeof w);
        ll t0 = 0;
        ll cap = 64 * B + 1000000;
        while (residual(w) > 504 && t0 < cap) { kstep(w, tmp); memcpy(w, tmp, sizeof w); t0++; }
        ll lam, mu, prev[N]; int hp;
        brent(w, &lam, &mu, prev, &hp);
        ll T;
        if (mu == 0) {
            /* fallback: full search from w0 (lemma-free) */
            t0 = 0;
            brent(w0, &lam, &mu, prev, &hp);
        }
        T = t0 + mu;
        /* difference clock: MUD in {T-1, T}; check whether d_{T-1} is cyclic */
        ll MUD = T;
        if (T >= 1 && hp) {
            ll a[N], b[N], qa[5], qb[5];
            memcpy(a, prev, sizeof a);     /* w_{T-1} */
            memcpy(b, prev, sizeof b);
            for (ll i = 0; i < lam; i++) { kstep(b, tmp); memcpy(b, tmp, sizeof b); }
            diffs(a, qa); diffs(b, qb);
            if (memcmp(qa, qb, sizeof qa) == 0) MUD = T - 1;
        } else if (T >= 1 && !hp) {
            /* mu==0 relative but t0>0 cannot happen here; full search gave mu==0 -> T==0 */
            MUD = T;
        }
        printf("%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld\n", B, G[0], G[1], G[2], G[3], G[4], T, MUD, lam, t0);
        fflush(stdout);
    }
    return 0;
}
