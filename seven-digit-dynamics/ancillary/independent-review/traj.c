/* Fresh-base tests by direct word-level digit subtraction (no difference-vector shortcut).
 *
 * kappa_B: sort the 7 digits, subtract ascending from descending by column subtraction with
 * borrow, keep 7 digits.  T_B(word) = index of the first word on a cycle (Brent's algorithm).
 * Real-map prediction (Theorem 1): iterate f(xi) = C(sort|2xi-1|) to the first visit to C(K),
 * pick the canonical rotation (x+z>2y, x>y+z), d = x-3z, clocks from the table.
 *
 * modes:
 *   clocks B NXI SEED     : NXI points xi ~ rho; for each, 512 words whose difference vectors fill
 *                           q0+[0,8)^3 (all residues mod 8); print T/B atoms and fast fraction.
 *   pop B NS SEED         : NS uniform random nonconstant words; print T_B values (one per line).
 *   limit NS SEED         : NS points xi ~ rho; print canonical data for the limit law.
 *   diffT B u v w         : T of a word with difference vector (u,v,w), plus h (=T-1 normally).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

static uint64_t rs;
static uint64_t rnd64(void) { uint64_t z = (rs += 0x9E3779B97F4A7C15ULL); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }
static double unif(void) { return (rnd64() >> 11) * (1.0 / 9007199254740992.0); }
static double gamma2(void) { return -log(1.0 - unif()) - log(1.0 - unif()); }

static int B;
typedef struct { int d[7]; } word;

static void kappa(const word *in, word *out) {
    int s[7];
    memcpy(s, in->d, sizeof s);
    for (int i = 1; i < 7; i++) { int x = s[i], j = i - 1; while (j >= 0 && s[j] < x) { s[j+1] = s[j]; j--; } s[j+1] = x; }
    /* hi = s[0..6] (high digit first), lo = reversed */
    int borrow = 0;
    for (int i = 6; i >= 0; i--) {
        int t = s[i] - s[6 - i] - borrow;
        if (t < 0) { t += B; borrow = 1; } else borrow = 0;
        out->d[i] = t;
    }
    if (borrow) { fprintf(stderr, "negative result?\n"); exit(1); }
}
static int weq(const word *a, const word *b) { return memcmp(a->d, b->d, sizeof a->d) == 0; }

/* Brent: returns mu (first cyclic index), sets *lam */
static long brent(const word *x0, long *lamout) {
    long power = 1, lam = 1;
    word tort = *x0, hare; kappa(x0, &hare);
    while (!weq(&tort, &hare)) {
        if (power == lam) { tort = hare; power *= 2; lam = 0; }
        word t; kappa(&hare, &t); hare = t; lam++;
    }
    word a = *x0, b = *x0;
    for (long i = 0; i < lam; i++) { word t; kappa(&b, &t); b = t; }
    long mu = 0;
    while (!weq(&a, &b)) { word t; kappa(&a, &t); a = t; kappa(&b, &t); b = t; mu++; }
    *lamout = lam;
    return mu;
}

static void word_from_diff(int u, int v, int w, word *W) {
    int d7 = (B - 1 - u) / 2, d1 = d7 + u;
    int d6 = d7 + (u - v) / 2, d2 = d6 + v;
    int d5 = d6 + (v - w) / 2, d3 = d5 + w;
    int d4 = d5 + w / 2;
    int s[7] = {d1, d2, d3, d4, d5, d6, d7};
    /* scramble order deterministically */
    int perm[7] = {3, 6, 0, 5, 1, 4, 2};
    for (int i = 0; i < 7; i++) W->d[i] = s[perm[i]];
}

/* ---------- real map ---------- */
static int inCK(long double x, long double y, long double z, long double eps) {
    long double a = 2 * x - 1, b = 2 * y - a, c = 2 * z - b;
    long double g[9] = {1 - a, a - b, b - c, c, 2 * a - 1, a + b - 1, a - b - c, b + 2 * c - 2 * a + 1, 2 - a - 2 * b - c};
    for (int i = 0; i < 9; i++) if (g[i] < eps) return 0;
    return 1;
}
static void freal(long double *x, long double *y, long double *z) {
    long double r[3] = {fabsl(2 * *x - 1), fabsl(2 * *y - 1), fabsl(2 * *z - 1)};
    for (int i = 1; i < 3; i++) { long double t = r[i]; int j = i - 1; while (j >= 0 && r[j] < t) { r[j+1] = r[j]; j--; } r[j+1] = t; }
    *x = (1 + r[0]) / 2; *y = (r[0] + r[1]) / 2; *z = (r[1] + r[2]) / 2;
}
/* returns steps to first visit of C(K) (strict, margin eps), -1 if none within maxit;
   fills canonical (cx,cy,cz), rotation index, t_f, t_s, side */
static int predict(long double x, long double y, long double z, long double eps, int maxit,
                   long double *cx, long double *cy, long double *cz, int *rot, long double *tf, long double *ts, int *side) {
    int n = 0;
    while (!inCK(x, y, z, eps)) { if (n >= maxit) return -1; freal(&x, &y, &z); n++; }
    int found = -1; long double px = x, py = y, pz = z;
    for (int k = 0; k < 3; k++) {
        if (px + pz > 2 * py && px > py + pz) { if (found >= 0) return -2; found = k; *cx = px; *cy = py; *cz = pz; }
        freal(&px, &py, &pz);
    }
    if (found < 0) return -3;
    *rot = found;
    long double X = *cx, Z = *cz, d = X - 3 * Z;
    if (d > 0) { *side = 1; *tf = (9 * X - 3 * Z - 6) / 8; *ts = (5 * X - 3 * Z - 3) / 4; }
    else { *side = -1; *tf = 3 * (X + Z - 1) / 4; *ts = (2 * X + 6 * Z - 3) / 4; }
    return n;
}
static void sample_rho(double *x, double *y, double *z) {
    double g[4]; double s = 0;
    for (int i = 0; i < 4; i++) { g[i] = gamma2(); s += g[i]; }
    for (int i = 0; i < 4; i++) g[i] /= s;
    /* gaps (1-x, x-y, y-z, z) */
    *z = g[3]; *y = g[3] + g[2]; *x = g[3] + g[2] + g[1];
}

int main(int argc, char **argv) {
    const char *mode = argv[1];
    if (!strcmp(mode, "clocks")) {
        B = atoi(argv[2]); int NXI = atoi(argv[3]); rs = strtoull(argv[4], 0, 10);
        int done = 0;
        while (done < NXI) {
            double x, y, z; sample_rho(&x, &y, &z);
            long double cx, cy, cz, tf, ts; int rot, side;
            int n = predict(x, y, z, 1e-6L, 200, &cx, &cy, &cz, &rot, &tf, &ts, &side);
            if (n < 0) continue;
            long double dd = cx - 3 * cz;
            if (fabsl(dd) < 0.03L) continue; /* need clock gap >> O(1)/B errors */
            int u0 = (int)(x * B) & ~7, v0 = (int)(y * B) & ~7, w0 = (int)(z * B) & ~7;
            if (!(u0 - 8 > v0 && v0 - 8 > w0 && w0 > 8 && u0 + 8 < B)) continue;
            int nfast = 0, nslow = 0, nother = 0; long double maxdev = 0; long Tmin = 1L << 60, Tmax = 0;
            for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) for (int k = 0; k < 8; k++) {
                word W; word_from_diff(u0 + i, v0 + j, w0 + k, &W);
                long lam; long T = brent(&W, &lam);
                if (T < Tmin) Tmin = T; if (T > Tmax) Tmax = T;
                long double df = T - B * tf, ds = T - B * ts;
                if (fabsl(df) < fabsl(ds)) { nfast++; if (fabsl(df) > maxdev) maxdev = fabsl(df); }
                else { nslow++; if (fabsl(ds) > maxdev) maxdev = fabsl(ds); }
                if (fabsl(df) > 400 && fabsl(ds) > 400) nother++;
            }
            printf("xi %.6f %.6f %.6f n %d side %+d d %.4Lf tf %.5Lf ts %.5Lf gapB %.0Lf | fast %d slow %d far %d frac %.5f maxdev %.1Lf Trange %ld %ld\n",
                   x, y, z, n, side, dd, tf, ts, (ts - tf) * B, nfast, nslow, nother, nfast / 512.0, maxdev, Tmin, Tmax);
            fflush(stdout);
            done++;
        }
    } else if (!strcmp(mode, "pop")) {
        B = atoi(argv[2]); long NS = atol(argv[3]); rs = strtoull(argv[4], 0, 10);
        for (long s = 0; s < NS; s++) {
            word W; int allsame;
            do { allsame = 1; for (int i = 0; i < 7; i++) { W.d[i] = (int)(rnd64() % (uint64_t)B); if (W.d[i] != W.d[0]) allsame = 0; } } while (allsame);
            long lam; long T = brent(&W, &lam);
            printf("%ld\n", T);
        }
    } else if (!strcmp(mode, "limit")) {
        long NS = atol(argv[2]); rs = strtoull(argv[3], 0, 10);
        long fail = 0;
        for (long s = 0; s < NS; s++) {
            double x, y, z; sample_rho(&x, &y, &z);
            long double cx, cy, cz, tf, ts; int rot, side;
            int n = predict(x, y, z, 0.0L, 10000, &cx, &cy, &cz, &rot, &tf, &ts, &side);
            if (n < 0) { fail++; continue; }
            printf("%d %d %.9Lf %.9Lf\n", n, side, tf, ts);
        }
        fprintf(stderr, "limit: %ld failures of %ld\n", fail, NS);
    } else if (!strcmp(mode, "limitstats")) {
        /* limit law nu_beta: mixture p_beta delta_tf + (1-p_beta) delta_ts, xi ~ rho */
        long NS = atol(argv[2]); rs = strtoull(argv[3], 0, 10);
        long fail = 0, npos = 0, nneg = 0, maxn = 0; long hist_n[64] = {0};
        double m1[2] = {0, 0}, m2[2] = {0, 0}, e4[2] = {0, 0}, em4[2] = {0, 0}, cdf[2][26];
        double mtf = 0, mts = 0, tmin = 1, tmax = 0, pfast[2] = {0, 0};
        memset(cdf, 0, sizeof cdf);
        for (long s = 0; s < NS; s++) {
            double x, y, z; sample_rho(&x, &y, &z);
            long double cx, cy, cz, tf, ts; int rot, side;
            int n = predict(x, y, z, 0.0L, 10000, &cx, &cy, &cz, &rot, &tf, &ts, &side);
            if (n < 0) { fail++; fprintf(stderr, "entry failure code %d at xi=(%.17g, %.17g, %.17g)\n", n, x, y, z); continue; }
            if (n > maxn) maxn = n; hist_n[n < 63 ? n : 63]++;
            if (side > 0) npos++; else nneg++;
            if (!(0 < tf && tf < ts && ts < 0.25)) fprintf(stderr, "clock order violated: %Lf %Lf\n", tf, ts);
            if (tf < tmin) tmin = tf; if (ts > tmax) tmax = ts;
            mtf += tf; mts += ts;
            for (int beta = 0; beta < 2; beta++) {
                double p = beta == 0 ? (side > 0 ? 0.5 : 0.75) : 0.375;
                pfast[beta] += p;
                m1[beta] += p * tf + (1 - p) * ts;
                m2[beta] += p * tf * tf + (1 - p) * ts * ts;
                e4[beta] += p * exp(4 * tf) + (1 - p) * exp(4 * ts);
                em4[beta] += p * exp(-4 * tf) + (1 - p) * exp(-4 * ts);
                for (int k = 0; k < 26; k++) { double g = 0.01 * k; cdf[beta][k] += p * (tf <= g) + (1 - p) * (ts <= g); }
            }
        }
        long N = NS - fail;
        printf("limit: samples=%ld entry failures=%ld side+ %.5f side- %.5f max entry steps %ld  t_f range min %.6f  t_s max %.6f\n", NS, fail, (double)npos / N, (double)nneg / N, maxn, tmin, tmax);
        printf("entry-steps histogram:"); for (int i = 0; i < 20; i++) printf(" %ld", hist_n[i]); printf("\n");
        for (int beta = 0; beta < 2; beta++) {
            double mu = m1[beta] / N;
            printf("beta=%d P(fast)=%.6f mean=%.6f var=%.6f E e^{4t}=%.6f E e^{-4t}=%.6f\n", beta, pfast[beta] / N, mu, m2[beta] / N - mu * mu, e4[beta] / N, em4[beta] / N);
            printf("beta=%d CDF", beta); for (int k = 0; k < 26; k++) printf(" %.5f", cdf[beta][k] / N); printf("\n");
        }
        printf("mean(odd)-mean(even)=%.6f\n", (m1[1] - m1[0]) / N);
    } else if (!strcmp(mode, "diffT")) {
        B = atoi(argv[2]);
        word W; word_from_diff(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), &W);
        long lam; long T = brent(&W, &lam);
        printf("%ld %ld\n", T, lam);
    }
    return 0;
}
