/* Monte Carlo of T_B at fresh bases, using ONLY the definition of the
 * Kaprekar step on numbers (no pair coordinates, no reduction from the paper).
 * For each sampled nonconstant string x, Floyd's cycle-finding algorithm gives
 * the exact tail length mu = first arrival time at the terminal cycle, and the
 * cycle (period, minimal element).
 * usage: mc_fresh B nsamples seed
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

static uint64_t B;
static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t next64(void) {      /* xoshiro256** */
    const uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
    return r;
}
static uint64_t splitmix(uint64_t *x) { uint64_t z = (*x += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }
static uint64_t unif(uint64_t n) {  /* unbiased in [0,n) by rejection */
    uint64_t lim = UINT64_MAX - UINT64_MAX % n;
    uint64_t r; do r = next64(); while (r >= lim);
    return r % n;
}
static uint64_t g(uint64_t x) {
    uint64_t d[5];
    for (int k = 0; k < 5; k++) { d[k] = x % B; x /= B; }
    for (int i = 1; i < 5; i++) { uint64_t t = d[i]; int j = i - 1; while (j >= 0 && d[j] < t) { d[j + 1] = d[j]; j--; } d[j + 1] = t; }
    uint64_t desc = 0, asc = 0;
    for (int k = 0; k < 5; k++) desc = desc * B + d[k];
    for (int k = 4; k >= 0; k--) asc = asc * B + d[k];
    return desc - asc;
}

#define MAXC 16
int main(int argc, char **argv) {
    B = strtoull(argv[1], 0, 10);
    long ns = atol(argv[2]);
    uint64_t seed = strtoull(argv[3], 0, 10);
    for (int i = 0; i < 4; i++) s[i] = splitmix(&seed);
    uint64_t cmin[MAXC]; int cper[MAXC]; long ccount[MAXC]; int nc = 0;
    double ts[4] = {1.0 / 25, 0.1, 0.2, 0.3};
    long hit[4] = {0}, hitc[4][MAXC] = {{0}};
    long maxT = 0; uint64_t argmax = 0; double sum = 0, sum2 = 0;
    uint64_t rep1 = 0; for (int k = 0; k < 5; k++) rep1 = rep1 * B + 1;
    for (long it = 0; it < ns; it++) {
        uint64_t x0;
        do { x0 = 0; for (int k = 0; k < 5; k++) x0 = x0 * B + unif(B); } while (x0 % rep1 == 0);
        /* Floyd */
        uint64_t t = g(x0), h = g(g(x0));
        while (t != h) { t = g(t); h = g(g(h)); }
        long mu = 0; t = x0;
        while (t != h) { t = g(t); h = g(h); mu++; }
        int lam = 1; uint64_t mn = t; h = g(t);
        while (h != t) { if (h < mn) mn = h; h = g(h); lam++; }
        int c;
        for (c = 0; c < nc; c++) if (cmin[c] == mn) break;
        if (c == nc) { cmin[nc] = mn; cper[nc] = lam; ccount[nc] = 0; nc++; }
        ccount[c]++;
        double y = (double)mu / (double)B;
        sum += y; sum2 += y * y;
        if (mu > maxT) { maxT = mu; argmax = x0; }
        for (int k = 0; k < 4; k++) if ((double)mu <= ts[k] * (double)B + 1e-9) { hit[k]++; hitc[k][c]++; }
    }
    printf("B %llu samples %ld seed %s\n", (unsigned long long)B, ns, argv[3]);
    for (int c = 0; c < nc; c++)
        printf("cycle period %d minval %llu count %ld frac %.6f se %.6f\n", cper[c], (unsigned long long)cmin[c], ccount[c],
               (double)ccount[c] / ns, sqrt((double)ccount[c] / ns * (1 - (double)ccount[c] / ns) / ns));
    for (int k = 0; k < 4; k++) {
        double p = (double)hit[k] / ns;
        printf("cdf t=%.4f p=%.6f se=%.6f", ts[k], p, sqrt(p * (1 - p) / ns));
        for (int c = 0; c < nc; c++) if (ccount[c] > 0) printf(" | per%d cond %.6f", cper[c], (double)hitc[k][c] / ccount[c]);
        printf("\n");
    }
    printf("mean T/B %.6f se %.6f\n", sum / ns, sqrt((sum2 / ns - (sum / ns) * (sum / ns)) / ns));
    printf("max T observed %ld at x0=%llu\n", maxT, (unsigned long long)argmax);
    return 0;
}
