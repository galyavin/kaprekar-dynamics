/* Fresh-base test of Theorem 5 (conditional overshoot law) -- reviewer code.
 * Samples difference pairs (u,v) in the band E_B with probability proportional
 * to the ordered-input weight (18), computes the pair arrival time h by Floyd's
 * algorithm on the pair map obtained from numerical subtraction, and reports
 * statistics of Y_B = (T_B - B/3)/sqrt(B), T_B = h+1.  It also checks the
 * paper's uniform formula (44), |e_B| <= 8, and (41), on every sample.
 * usage: overshoot_mc B alpha beta gamma zeta nsamples seed
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef unsigned __int128 u128;
static long B;
static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t next64(void) { const uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17; s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45); return r; }
static uint64_t splitmix(uint64_t *x) { uint64_t z = (*x += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }
static double unif01(void) { return (next64() >> 11) * (1.0 / 9007199254740992.0); }

static inline uint64_t enc(long u, long v) { return ((uint64_t)u << 32) | (uint64_t)v; }
static uint64_t f(uint64_t q) {
    long u = (long)(q >> 32), v = (long)(q & 0xffffffffu);
    u128 Bb = (u128)B;
    u128 N = (Bb * Bb * Bb * Bb - 1) * (u128)u + (Bb * Bb * Bb - Bb) * (u128)v;
    long d[5];
    for (int k = 0; k < 5; k++) { d[k] = (long)(N % Bb); N /= Bb; }
    for (int i = 1; i < 5; i++) { long x = d[i]; int j = i - 1; while (j >= 0 && d[j] < x) { d[j + 1] = d[j]; j--; } d[j + 1] = x; }
    return enc(d[0] - d[4], d[1] - d[3]);
}
static double w18(long u, long v) { return (double)(B - u) * (120.0 * (double)v * (double)(u - v) - 20.0); }

int main(int argc, char **argv) {
    B = atol(argv[1]);
    double al = atof(argv[2]), be = atof(argv[3]), ga = atof(argv[4]), ze = atof(argv[5]);
    long ns = atol(argv[6]); uint64_t seed = strtoull(argv[7], 0, 10);
    for (int i = 0; i < 4; i++) s[i] = splitmix(&seed);
    double sq = sqrt((double)B);
    long rlo = (long)ceil(al * sq - 1e-12), rhi = (long)floor(be * sq + 1e-12);
    long vlo = (long)ceil(ga * B - 1e-12), vhi = (long)floor(ze * B + 1e-12);
    /* exact band probability */
    long double mass = 0; double wmax = 0;
    for (long r = rlo; r <= rhi; r++) for (long v = vlo; v <= vhi; v++) {
        long u = B - r; if (!(0 < v && v < u)) continue;
        double w = w18(u, v); mass += w; if (w > wmax) wmax = w;
    }
    long double tot = (long double)B * B * B * B * B - B;
    long H = (B + 1) / 2;
    double ys[7] = {-1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25};
    long cnt[7] = {0};
    double sy = 0, sy2 = 0; double maxe = 0; long bad41 = 0, badbranch = 0;
    for (long it = 0; it < ns; it++) {
        long r, v, u;
        for (;;) {
            r = rlo + (long)(unif01() * (rhi - rlo + 1)); if (r > rhi) r = rhi;
            v = vlo + (long)(unif01() * (vhi - vlo + 1)); if (v > vhi) v = vhi;
            u = B - r;
            if (unif01() * wmax <= w18(u, v)) break;
        }
        uint64_t x0 = enc(u, v);
        /* Floyd: h = tail length of the pair orbit */
        uint64_t t = f(x0), hh = f(f(x0));
        while (t != hh) { t = f(t); hh = f(f(hh)); }
        long mu = 0; t = x0; while (t != hh) { t = f(t); hh = f(hh); mu++; }
        long T = mu + 1;
        double Y = (T - B / 3.0) / sq;
        sy += Y; sy2 += Y * Y;
        for (int k = 0; k < 7; k++) if (Y <= ys[k]) cnt[k]++;
        /* high segment: simulate exit, compare with (40)-(41) and (44) */
        long n = 0, uu = u, vv = v;
        while (vv > H) {
            uint64_t q = f(enc(uu, vv)); long u2 = (long)(q >> 32), v2 = (long)(q & 0xffffffffu);
            if (u2 != uu - 1 || v2 != uu + vv - B + 1) badbranch++;
            uu = u2; vv = v2; n++;
        }
        long R = r + n, nu = H - vv;
        if (nu < 0 || nu > R - 3) bad41++;
        double mn = nu < R - nu ? nu : R - nu;
        double e = T - n - (2.0 / 3.0) * (B / 2.0 - R - mn);
        if (fabs(e) > maxe) maxe = fabs(e);
    }
    double my = sy / ns, sdy = sqrt(sy2 / ns - my * my);
    printf("B %ld band lambda[%g,%g] z[%g,%g] samples %ld | B*P(E_B) = %.6Lf\n", B, al, be, ga, ze, ns, mass / tot * B);
    printf("E[Y|E] = %.5f  (se %.5f)  sd %.5f\n", my, sdy / sqrt((double)ns), sdy);
    for (int k = 0; k < 7; k++) { double p = (double)cnt[k] / ns; printf("P(Y<=%.2f|E) = %.5f (se %.5f)\n", ys[k], p, sqrt(p * (1 - p) / ns)); }
    printf("max |e_B| in (44) = %.4f (claim <= 8); violations of (41): %ld; high-branch mismatches: %ld\n", maxe, bad41, badbranch);
    return 0;
}
