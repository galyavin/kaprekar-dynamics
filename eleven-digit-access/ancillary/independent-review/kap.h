/* Independent encoding of eleven-digit Kaprekar subtraction.
 * Written from the paper's definitions only (Section 1): sort digits,
 * subtract increasing arrangement from decreasing one positionally with
 * borrow, keep leading zeros.  No use of eq. (digits) here. */
#ifndef KAP_H
#define KAP_H
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#define ND 11
typedef int64_t i64;

static inline void sort_desc(i64 *a, int n) {
    for (int i = 1; i < n; i++) {
        i64 x = a[i]; int j = i - 1;
        while (j >= 0 && a[j] < x) { a[j + 1] = a[j]; j--; }
        a[j + 1] = x;
    }
}

/* one Kaprekar step on a word w[0..10] (w[0] most significant) */
static inline void kap_word(const i64 *w, i64 *out, i64 B) {
    i64 D[ND];
    memcpy(D, w, sizeof D);
    sort_desc(D, ND);
    i64 borrow = 0;
    for (int k = ND - 1; k >= 0; k--) {
        i64 t = D[k] - D[ND - 1 - k] - borrow;   /* A[k] = D[10-k] */
        if (t < 0) { t += B; borrow = 1; } else borrow = 0;
        out[k] = t;
    }
    /* borrow must be 0 at the end since D >= A as numbers */
    if (borrow) abort();
}

/* differences q[0..4] (q_1..q_5) of a word */
static inline void word_diffs(const i64 *w, i64 *q) {
    i64 D[ND];
    memcpy(D, w, sizeof D);
    sort_desc(D, ND);
    for (int i = 0; i < 5; i++) q[i] = D[i] - D[ND - 1 - i];
}

/* F_B on difference states, via the canonical word (q1..q5,0,...,0) */
static inline void F_state(const i64 *q, i64 *qn, i64 B) {
    i64 w[ND] = {0}, o[ND];
    for (int i = 0; i < 5; i++) w[i] = q[i];
    kap_word(w, o, B);
    word_diffs(o, qn);
}

static inline void gaps(const i64 *q, i64 B, i64 *G) {
    G[0] = B - q[0];
    for (int i = 1; i < 5; i++) G[i] = q[i - 1] - q[i];
    G[5] = q[4];
}
static inline i64 min_gap(const i64 *q, i64 B) {
    i64 G[6]; gaps(q, B, G);
    i64 m = G[0];
    for (int i = 1; i < 6; i++) if (G[i] < m) m = G[i];
    return m;
}
static inline i64 max_gap(const i64 *q, i64 B) {
    i64 G[6]; gaps(q, B, G);
    i64 m = G[0];
    for (int i = 1; i < 6; i++) if (G[i] > m) m = G[i];
    return m;
}
static inline i64 MB(i64 B) { return (B + 7) / 7 + 64; } /* ceil((B+1)/7)+64 */
static inline int in_S(const i64 *q, i64 B) { return min_gap(q, B) >= MB(B); }
/* centered residual R = -min(6G_i - B) = B - 6 min G */
static inline i64 resid(const i64 *q, i64 B) { return B - 6 * min_gap(q, B); }

/* Paired rule P_B from its definition (Section 4.2): positive innermost source,
 * level i = (L_i, R_{6-i}), L_i = q_i - [i=5], R_i = B - q_i - [i>1];
 * B-1 first, sort within each pair, keep level order, symmetric differences.
 * Here the digits L,R are taken from eq. (digits) -- used only to test P_B vs F_B. */
static inline void P_state(const i64 *q, i64 *qn, i64 B) {
    i64 L[6], R[6], s[ND];
    for (int i = 1; i <= 5; i++) {
        L[i] = q[i - 1] - (i == 5);
        R[i] = B - q[i - 1] - (i > 1);
    }
    s[0] = B - 1;
    for (int i = 1; i <= 5; i++) {
        i64 a = L[i], b = R[6 - i];
        s[2 * i - 1] = a > b ? a : b;
        s[2 * i] = a > b ? b : a;
    }
    for (int i = 0; i < 5; i++) qn[i] = s[i] - s[ND - 1 - i];
}

/* xoshiro256** */
typedef struct { uint64_t s[4]; } rng_t;
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rng_next(rng_t *r) {
    uint64_t *s = r->s;
    const uint64_t result = rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl(s[3], 45);
    return result;
}
static inline void rng_seed(rng_t *r, uint64_t seed) {
    uint64_t z = seed;
    for (int i = 0; i < 4; i++) {           /* splitmix64 */
        z += 0x9e3779b97f4a7c15ULL;
        uint64_t x = z;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        r->s[i] = x ^ (x >> 31);
    }
}
static inline i64 rng_below(rng_t *r, i64 n) {   /* unbiased */
    uint64_t lim = UINT64_MAX - (UINT64_MAX % (uint64_t)n);
    uint64_t x;
    do { x = rng_next(r); } while (x >= lim);
    return (i64)(x % (uint64_t)n);
}
#endif
