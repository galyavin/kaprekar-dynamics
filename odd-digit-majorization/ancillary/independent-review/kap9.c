/* Independent nine-digit Kaprekar simulator (direct column subtraction).
 *
 * The transition is the definition: realise q by the sorted string
 * (q1,q2,q3,q4,q4,0,0,0,0), subtract the ascending arrangement from the
 * descending one with column borrows, sort the nine result digits and take
 * symmetric differences.  No borrow formula and no guarded branches.
 *
 * Modes:
 *   descent B            exhaustive check of eq. (descent) at base B
 *   maxT B               exact max_N T_B and distribution by memoised orbit graph
 *   orbits               stdin lines "B q1 q2 q3 q4" -> "mu lam z1 z2 z3 z4 pos paired"
 *   uniform B n seed     uniform random strings N: prints "T z1 z2 z3 z4 x1 x2 x3 x4"
 *   descrand n seed      random large-B stratified test of eq. (descent)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef long long ll;
typedef struct { ll q[4]; } st;

static inline void sort9desc(ll *d) {
    for (int i = 1; i < 9; i++) {
        ll x = d[i]; int j = i - 1;
        while (j >= 0 && d[j] < x) { d[j + 1] = d[j]; j--; }
        d[j + 1] = x;
    }
}

/* K_B on a digit string (most significant first), column subtraction */
static inline void kap_digits(const ll *digits, ll B, ll *out) {
    ll desc[9];
    memcpy(desc, digits, sizeof desc);
    sort9desc(desc);
    ll borrow = 0;
    for (int j = 8; j >= 0; j--) {
        ll t = desc[j] - desc[8 - j] - borrow;
        if (t < 0) { t += B; borrow = 1; } else borrow = 0;
        out[j] = t;
    }
    if (borrow) { fprintf(stderr, "negative difference?\n"); exit(1); }
}

static inline void qvec(const ll *digits, ll *q) {
    ll d[9];
    memcpy(d, digits, sizeof d);
    sort9desc(d);
    for (int i = 0; i < 4; i++) q[i] = d[i] - d[8 - i];
}

static inline void realize(const ll *q, ll *dig) {
    dig[0] = q[0]; dig[1] = q[1]; dig[2] = q[2]; dig[3] = q[3]; dig[4] = q[3];
    dig[5] = dig[6] = dig[7] = dig[8] = 0;
}

static inline st F(st s, ll B) {
    ll dig[9], out[9];
    st r;
    realize(s.q, dig);
    kap_digits(dig, B, out);
    qvec(out, r.q);
    return r;
}

static inline int eq(st a, st b) { return a.q[0] == b.q[0] && a.q[1] == b.q[1] && a.q[2] == b.q[2] && a.q[3] == b.q[3]; }

static inline ll Gfun(st s, ll B) {
    ll g[5] = {B - s.q[0], s.q[0] - s.q[1], s.q[1] - s.q[2], s.q[2] - s.q[3], s.q[3]};
    ll G = 0;
    for (int i = 0; i < 5; i++) for (int j = i + 1; j < 5; j++) G += g[i] > g[j] ? g[i] - g[j] : g[j] - g[i];
    return G;
}

static inline ll Pfun(st s, ll B) {
    ll g[5] = {B - s.q[0], s.q[0] - s.q[1], s.q[1] - s.q[2], s.q[2] - s.q[3], s.q[3]};
    ll m = g[0], M = g[0];
    for (int i = 1; i < 5; i++) { if (g[i] < m) m = g[i]; if (g[i] > M) M = g[i]; }
    ll mn = M < 2 * m ? M : 2 * m;
    return B - 4 * m - mn;
}

static inline int paired(st s, ll B) {
    if (s.q[3] <= 0) return 1; /* only defined for positive states */
    ll a[4] = {s.q[0], s.q[1], s.q[2], s.q[3] - 1};
    ll b[4] = {B - 1 - s.q[3], B - 1 - s.q[2], B - 1 - s.q[1], B - s.q[0]};
    for (int i = 0; i < 3; i++) {
        ll mn = a[i] < b[i] ? a[i] : b[i];
        ll mx = a[i + 1] > b[i + 1] ? a[i + 1] : b[i + 1];
        if (mn < mx) return 0;
    }
    return 1;
}

/* Brent: mu = first cyclic index, lam = period. Also flags along states 0..mu+lam-1. */
static void orbit(st x0, ll B, ll *mu_o, ll *lam_o, st *succ, int *allpos, int *allpaired) {
    ll power = 1, lam = 1;
    st tort = x0, hare = F(x0, B);
    while (!eq(tort, hare)) {
        if (power == lam) { tort = hare; power *= 2; lam = 0; }
        hare = F(hare, B);
        lam++;
    }
    tort = x0; hare = x0;
    for (ll i = 0; i < lam; i++) hare = F(hare, B);
    ll mu = 0;
    int ap = 1, apr = 1;
    while (!eq(tort, hare)) {
        /* flags exclude the input state (index 0): pairing is tested from the first image on */
        if (mu > 0) { if (tort.q[3] <= 0) ap = 0; else if (!paired(tort, B)) apr = 0; }
        tort = F(tort, B); hare = F(hare, B); mu++;
    }
    /* tort = first cyclic state; walk the cycle for flags */
    st c = tort;
    for (ll i = 0; i < lam; i++) {
        if (c.q[3] <= 0) ap = 0; else if (!paired(c, B)) apr = 0;
        c = F(c, B);
    }
    *mu_o = mu; *lam_o = lam; *succ = F(tort, B); *allpos = ap; *allpaired = apr;
}

/* ---------------------------------------------------------------- RNG */
static uint64_t rs = 88172645463325252ULL;
static inline uint64_t xr(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline ll rnd(ll n) { return (ll)(xr() % (uint64_t)n); }

/* ---------------------------------------------------------------- modes */
static int mode_descent(ll B) {
    ll checked = 0, viol = 0, allpos = 0, nodesc_lowP = 0;
    ll firststep[10] = {0};
    ll minP_nodesc = 1LL << 60;
    st s;
    for (s.q[0] = 1; s.q[0] < B; s.q[0]++)
    for (s.q[1] = 1; s.q[1] <= s.q[0]; s.q[1]++)
    for (s.q[2] = 1; s.q[2] <= s.q[1]; s.q[2]++)
    for (s.q[3] = 1; s.q[3] <= s.q[2]; s.q[3]++) {
        ll G0 = Gfun(s, B), P = Pfun(s, B);
        st t = s; int pos = 1, first = 0;
        for (int i = 1; i <= 9; i++) {
            t = F(t, B);
            if (t.q[3] <= 0) { pos = 0; break; }
            if (!first && Gfun(t, B) < G0) first = i;
        }
        if (!pos) continue;
        allpos++;
        if (P > 100) {
            checked++;
            if (!first) {
                viol++;
                if (viol <= 10) printf("VIOLATION B=%lld q=(%lld,%lld,%lld,%lld) G=%lld P=%lld\n", B, s.q[0], s.q[1], s.q[2], s.q[3], G0, P);
            } else firststep[first]++;
        } else if (!first) {
            nodesc_lowP++;
            if (P < minP_nodesc) minP_nodesc = P;
        }
        if (P <= 100 && !first) { /* track max P among no-descent states */ }
    }
    printf("B=%lld all-ten-positive=%lld checked(P>100)=%lld violations=%lld  firststep:", B, allpos, checked, viol);
    for (int i = 1; i <= 9; i++) printf(" %lld", firststep[i]);
    printf("  no-descent states with P<=100: %lld\n", nodesc_lowP);
    return 0;
}

/* max over no-descent states of P, regardless of threshold: margin probe */
static int mode_descmargin(ll B) {
    ll maxP = -(1LL << 60); st arg = {{0, 0, 0, 0}};
    ll cnt = 0;
    st s;
    for (s.q[0] = 1; s.q[0] < B; s.q[0]++)
    for (s.q[1] = 1; s.q[1] <= s.q[0]; s.q[1]++)
    for (s.q[2] = 1; s.q[2] <= s.q[1]; s.q[2]++)
    for (s.q[3] = 1; s.q[3] <= s.q[2]; s.q[3]++) {
        ll G0 = Gfun(s, B);
        st t = s; int pos = 1, first = 0;
        for (int i = 1; i <= 9; i++) {
            t = F(t, B);
            if (t.q[3] <= 0) { pos = 0; break; }
            if (Gfun(t, B) < G0) { first = i; break; }
        }
        if (!pos || first) continue;
        cnt++;
        ll P = Pfun(s, B);
        if (P > maxP) { maxP = P; arg = s; }
    }
    printf("B=%lld no-descent(all ten positive) states=%lld  max P among them=%lld at q=(%lld,%lld,%lld,%lld)\n",
           B, cnt, maxP, arg.q[0], arg.q[1], arg.q[2], arg.q[3]);
    return 0;
}

static inline ll C2(ll n) { return n < 2 ? 0 : n * (n - 1) / 2; }
static inline ll C3(ll n) { return n < 3 ? 0 : n * (n - 1) * (n - 2) / 6; }
static inline ll C4(ll n) { return n < 4 ? 0 : n * (n - 1) * (n - 2) * (n - 3) / 24; }
static inline ll idx(st s) { return C4(s.q[0] + 3) + C3(s.q[1] + 2) + C2(s.q[2] + 1) + s.q[3]; }

static int mode_maxT(ll B) {
    ll N = C4(B + 3);
    int32_t *h = malloc(sizeof(int32_t) * N);
    if (!h) { fprintf(stderr, "oom\n"); return 1; }
    for (ll i = 0; i < N; i++) h[i] = -1;
    ll cap = 1 << 20;
    st *path = malloc(sizeof(st) * cap);
    ll *pidx = malloc(sizeof(ll) * cap);
    int32_t ONPATH = -2;
    st s;
    for (s.q[0] = 0; s.q[0] < B; s.q[0]++)
    for (s.q[1] = 0; s.q[1] <= s.q[0]; s.q[1]++)
    for (s.q[2] = 0; s.q[2] <= s.q[1]; s.q[2]++)
    for (s.q[3] = 0; s.q[3] <= s.q[2]; s.q[3]++) {
        ll i0 = idx(s);
        if (h[i0] >= 0) continue;
        ll len = 0; st t = s; ll it = i0;
        while (h[it] == -1) {
            if (len >= cap) { cap *= 2; path = realloc(path, sizeof(st) * cap); pidx = realloc(pidx, sizeof(ll) * cap); }
            path[len] = t; pidx[len] = it; len++;
            h[it] = ONPATH;
            t = F(t, B); it = idx(t);
        }
        ll base;
        if (h[it] == ONPATH) { /* new cycle: from position of it */
            ll p = len - 1;
            while (pidx[p] != it) p--;
            for (ll j = p; j < len; j++) h[pidx[j]] = 0;
            len = p; base = 0;
        } else base = h[it];
        for (ll j = len - 1; j >= 0; j--) { base++; h[pidx[j]] = (int32_t)base; }
    }
    /* T_B(N) = 1 + h(q(N)) for transient q(N); q=0 excluded (constant strings) */
    ll maxh = 0; st arg = {{0, 0, 0, 0}};
    ll ncyc = 0;
    for (s.q[0] = 1; s.q[0] < B; s.q[0]++)
    for (s.q[1] = 0; s.q[1] <= s.q[0]; s.q[1]++)
    for (s.q[2] = 0; s.q[2] <= s.q[1]; s.q[2]++)
    for (s.q[3] = 0; s.q[3] <= s.q[2]; s.q[3]++) {
        ll v = h[idx(s)];
        if (v == 0) ncyc++;
        if (v > maxh) { maxh = v; arg = s; }
    }
    printf("B=%lld states=%lld cyclic=%lld max_h=%lld maxT=%lld maxT/B=%.5f argmax q=(%lld,%lld,%lld,%lld)\n",
           B, N, ncyc, maxh, maxh + 1, (double)(maxh + 1) / B, arg.q[0], arg.q[1], arg.q[2], arg.q[3]);
    free(h); free(path); free(pidx);
    return 0;
}

static int mode_orbits(void) {
    ll B, q[4];
    while (scanf("%lld %lld %lld %lld %lld", &B, &q[0], &q[1], &q[2], &q[3]) == 5) {
        st s = {{q[0], q[1], q[2], q[3]}}, z;
        ll mu, lam; int ap, apr;
        orbit(s, B, &mu, &lam, &z, &ap, &apr);
        printf("%lld %lld %lld %lld %lld %lld %d %d\n", mu, lam, z.q[0], z.q[1], z.q[2], z.q[3], ap, apr);
    }
    return 0;
}

static int mode_uniform(ll B, ll n, uint64_t seed) {
    rs ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) xr();
    for (ll k = 0; k < n; k++) {
        ll dig[9];
        int cons = 1;
        do {
            for (int i = 0; i < 9; i++) dig[i] = rnd(B);
            cons = 1; for (int i = 1; i < 9; i++) if (dig[i] != dig[0]) cons = 0;
        } while (cons);
        st s; qvec(dig, s.q);
        st z; ll mu, lam; int ap, apr;
        orbit(s, B, &mu, &lam, &z, &ap, &apr);
        ll T;
        st Z;
        if (mu > 0) { T = mu + 1; Z = z; }
        else {
            /* q(N) cyclic: N cyclic iff N == number of the cycle predecessor */
            st c = s; for (ll i = 0; i < lam - 1; i++) c = F(c, B);
            ll dg[9], out[9]; realize(c.q, dg); kap_digits(dg, B, out);
            int same = 1; for (int i = 0; i < 9; i++) if (out[i] != dig[i]) same = 0;
            if (same) { T = 0; Z = s; } else { T = 1; Z = F(s, B); }
        }
        printf("%lld %lld %lld %lld %lld %lld %lld %lld %lld\n", T, Z.q[0], Z.q[1], Z.q[2], Z.q[3], s.q[0], s.q[1], s.q[2], s.q[3]);
    }
    return 0;
}

/* random large-B test of the descent assertion, stratified */
static int mode_descrand(ll n, uint64_t seed) {
    rs ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) xr();
    ll checked = 0, viol = 0, skipped_pos = 0, skipped_P = 0;
    ll firststep[10] = {0};
    for (ll k = 0; k < n; k++) {
        ll B;
        int bt = (int)rnd(4);
        if (bt == 0) B = 101 + rnd(900);
        else if (bt == 1) B = 1000 + rnd(100000);
        else if (bt == 2) B = 100000 + rnd(10000000);
        else B = 1000000000LL + rnd(1000000000LL);
        ll g[5];
        int stype = (int)rnd(5);
        /* gap vectors summing to B with q4>=1 (g[4]>=1) */
        if (stype == 0) { /* uniform on the simplex (sorted uniforms) */
            ll c[4]; for (int i = 0; i < 4; i++) c[i] = rnd(B + 1);
            for (int i = 1; i < 4; i++) for (int j = i; j > 0 && c[j] < c[j - 1]; j--) { ll t = c[j]; c[j] = c[j - 1]; c[j - 1] = t; }
            g[0] = c[0]; g[1] = c[1] - c[0]; g[2] = c[2] - c[1]; g[3] = c[3] - c[2]; g[4] = B - c[3];
        } else if (stype == 1) { /* near-equal gaps: P just above 100 */
            ll m = B / 5 - 40 - rnd(40);
            for (int i = 0; i < 5; i++) g[i] = m;
            ll rest = B - 5 * m;
            while (rest > 0) { int i = (int)rnd(5); ll a = 1 + rnd(rest < 7 ? rest : 7); g[i] += a; rest -= a; }
        } else if (stype == 2) { /* small minimum gap */
            ll m = rnd(6);
            ll c[3]; for (int i = 0; i < 3; i++) c[i] = rnd(B - 2 * m + 1);
            for (int i = 1; i < 3; i++) for (int j = i; j > 0 && c[j] < c[j - 1]; j--) { ll t = c[j]; c[j] = c[j - 1]; c[j - 1] = t; }
            ll rem[4] = {c[0], c[1] - c[0], c[2] - c[1], B - 2 * m - c[2]};
            int pos = (int)rnd(5), pos2 = (int)rnd(5); while (pos2 == pos) pos2 = (int)rnd(5);
            int r = 0; for (int i = 0; i < 5; i++) { if (i == pos) g[i] = m; else if (i == pos2) g[i] = m + rem[r++] * 0; else g[i] = 0; }
            r = 0; for (int i = 0; i < 5; i++) if (i != pos) { g[i] = (i == pos2 ? m : 0) + rem[r++]; }
        } else if (stype == 3) { /* ties among gaps */
            ll a = rnd(B / 3 + 1), b = rnd(B / 3 + 1);
            ll cpat[5] = {a, a, b, b, 0}; ll sum = 2 * a + 2 * b; cpat[4] = B - sum;
            for (int i = 4; i > 0; i--) { int j = (int)rnd(i + 1); ll t = cpat[i]; cpat[i] = cpat[j]; cpat[j] = t; }
            for (int i = 0; i < 5; i++) g[i] = cpat[i];
        } else { /* near the core point (0.86,0.64,0.40,0.18) and Q(4/5) */
            double y[4];
            if (rnd(2)) { y[0] = 0.86; y[1] = 0.64; y[2] = 0.40; y[3] = 0.18; }
            else { y[0] = 0.8; y[1] = 0.6; y[2] = 0.4; y[3] = 0.2; }
            ll q[4]; ll sp = B / 50 + 1;
            for (int i = 0; i < 4; i++) q[i] = (ll)(y[i] * B) + rnd(2 * sp + 1) - sp;
            for (int i = 1; i < 4; i++) for (int j = i; j > 0 && q[j] > q[j - 1]; j--) { ll t = q[j]; q[j] = q[j - 1]; q[j - 1] = t; }
            if (q[0] > B - 1) q[0] = B - 1;
            g[0] = B - q[0]; g[1] = q[0] - q[1]; g[2] = q[1] - q[2]; g[3] = q[2] - q[3]; g[4] = q[3];
        }
        int okg = 1; for (int i = 0; i < 5; i++) if (g[i] < 0) okg = 0;
        if (!okg || g[4] < 1 || g[0] < 1) { k--; continue; }
        st s = {{B - g[0], B - g[0] - g[1], B - g[0] - g[1] - g[2], g[4]}};
        if (Pfun(s, B) <= 100) { skipped_P++; continue; }
        ll G0 = Gfun(s, B);
        st t = s; int pos = 1, first = 0;
        for (int i = 1; i <= 9; i++) {
            t = F(t, B);
            if (t.q[3] <= 0) { pos = 0; break; }
            if (!first && Gfun(t, B) < G0) first = i;
        }
        if (!pos) { skipped_pos++; continue; }
        checked++;
        if (!first) { viol++; if (viol < 10) printf("VIOLATION B=%lld q=(%lld,%lld,%lld,%lld)\n", B, s.q[0], s.q[1], s.q[2], s.q[3]); }
        else firststep[first]++;
    }
    printf("descrand: checked=%lld violations=%lld skipped(P<=100)=%lld skipped(nonpositive in 10)=%lld firststep:", checked, viol, skipped_P, skipped_pos);
    for (int i = 1; i <= 9; i++) printf(" %lld", firststep[i]);
    printf("\n");
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage\n"); return 1; }
    if (!strcmp(argv[1], "descent")) return mode_descent(atoll(argv[2]));
    if (!strcmp(argv[1], "descmargin")) return mode_descmargin(atoll(argv[2]));
    if (!strcmp(argv[1], "maxT")) return mode_maxT(atoll(argv[2]));
    if (!strcmp(argv[1], "orbits")) return mode_orbits();
    if (!strcmp(argv[1], "uniform")) return mode_uniform(atoll(argv[2]), atoll(argv[3]), strtoull(argv[4], 0, 10));
    if (!strcmp(argv[1], "descrand")) return mode_descrand(atoll(argv[2]), strtoull(argv[3], 0, 10));
    fprintf(stderr, "unknown mode\n");
    return 1;
}
