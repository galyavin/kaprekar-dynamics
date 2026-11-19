/* Fresh-base Monte Carlo of the headline laws (Theorem 1), exact integer arithmetic.
 * Uniform nonconstant ordered 11-digit words in base B; Kaprekar subtraction on words.
 * Per sample:
 *   sigma  = first j with differences in S_B (min gap >= ceil((B+1)/7)+64)
 *   mu,lam = first cyclic difference index and period (Brent on q_j, q_{j+1}=F(q_j))
 *   T      = first cyclic WORD index (mu or mu+1, decided by comparing words)
 *   dist   = B*|Z_B - (5,4,3,2,1)/6|_inf at the first cyclic word
 * and, along j in [sigma, mu+lam): pairing F_B==P_B and q5>0, R_{j+1}<=R_j+6,
 * lower hit within 32 when R_j>120, envelope R_j<=R_sigma+384 before first R<=120,
 * R_j<=504 after it, and R<=504 on the cycle (Section 6.2 / Lemma 12 consequences).
 * Usage: sim B nsamples seed [capfactor]
 */
#include <stdio.h>
#include "kap.h"

static i64 B;
static int eq5(const i64 *a, const i64 *b) { return !memcmp(a, b, 5 * sizeof(i64)); }
static void canon_word(const i64 *q, i64 *w) { memset(w, 0, ND * sizeof(i64)); for (int i = 0; i < 5; i++) w[i] = q[i]; }

int main(int argc, char **argv) {
    B = atoll(argv[1]); long ns = atol(argv[2]); uint64_t seed = strtoull(argv[3], 0, 10);
    i64 capf = argc > 4 ? atoll(argv[4]) : 5000;
    i64 cap = capf * B;
    rng_t rg; rng_seed(&rg, seed * 1000003ULL + (uint64_t)B);
    static const i64 v[5] = {5, 4, 3, 2, 1};
    i64 *Rbuf = NULL; size_t Rcap = 0;
    for (long smp = 0; smp < ns; smp++) {
        i64 w0[ND];
        for (;;) {
            int nonconst = 0;
            for (int k = 0; k < ND; k++) { w0[k] = rng_below(&rg, B); if (k && w0[k] != w0[0]) nonconst = 1; }
            if (nonconst) break;
        }
        i64 q0[5]; word_diffs(w0, q0);
        /* consistency: word iteration vs F_state for a few steps */
        int consist = 1;
        { i64 w[ND], wn[ND], q[5], qf[5], qn[5]; memcpy(w, w0, sizeof w); memcpy(q, q0, sizeof q);
          for (int k = 0; k < 50; k++) { kap_word(w, wn, B); word_diffs(wn, qn); F_state(q, qf, B);
            if (!eq5(qn, qf)) consist = 0; memcpy(w, wn, sizeof w); memcpy(q, qn, sizeof q); } }
        /* Brent */
        i64 tort[5], hare[5]; memcpy(tort, q0, sizeof tort);
        i64 power = 1, lam = 1, hidx = 1, sigma = -1;
        if (in_S(q0, B)) sigma = 0;
        F_state(tort, hare, B);
        int capped = 0;
        while (!eq5(tort, hare)) {
            if (sigma < 0 && in_S(hare, B)) sigma = hidx;
            if (power == lam) { memcpy(tort, hare, sizeof tort); power *= 2; lam = 0; }
            i64 t[5]; F_state(hare, t, B); memcpy(hare, t, sizeof hare); lam++; hidx++;
            if (hidx > cap) { capped = 1; break; }
        }
        if (capped) { printf("B=%lld s=%ld CAPPED sigma=%lld\n", (long long)B, smp, (long long)sigma); continue; }
        /* mu */
        i64 a[5], b[5], pa[5], pb[5]; memcpy(a, q0, sizeof a); memcpy(b, q0, sizeof b);
        for (i64 k = 0; k < lam; k++) { memcpy(pb, b, sizeof pb); i64 t[5]; F_state(b, t, B); memcpy(b, t, sizeof b); }
        /* pb = q_{lam-1} */
        i64 mu = 0;
        memcpy(pa, a, sizeof pa);
        while (!eq5(a, b)) {
            memcpy(pa, a, sizeof pa); memcpy(pb, b, sizeof pb);
            i64 t[5]; F_state(a, t, B); memcpy(a, t, sizeof a); F_state(b, t, B); memcpy(b, t, sizeof b); mu++;
        }
        /* now a = q_mu (cyclic); if mu>0: pa = q_{mu-1}, pb = q_{mu+lam-1} */
        i64 T;
        {
            i64 wm[ND], wml[ND], c[ND];
            if (mu == 0) { memcpy(wm, w0, sizeof wm); canon_word(pb, c); kap_word(c, wml, B); }
            else { canon_word(pa, c); kap_word(c, wm, B); canon_word(pb, c); kap_word(c, wml, B); }
            T = memcmp(wm, wml, sizeof wm) ? mu + 1 : mu;
        }
        i64 qT[5]; memcpy(qT, a, sizeof qT);
        if (T == mu + 1) { i64 t[5]; F_state(qT, t, B); memcpy(qT, t, sizeof qT); }
        i64 dist6 = 0;
        for (int i = 0; i < 5; i++) { i64 d = llabs(6 * qT[i] - B * v[i]); if (d > dist6) dist6 = d; }
        /* post-sigma checks */
        long pair_bad = 0, jump_bad = 0, hit_bad = 0, env_bad = 0, cyc_bad = 0;
        i64 Rsig = -1, maxR_after = -1, maxjump = -1000000; i64 hitmax = 0;
        if (sigma >= 0) {
            i64 q[5]; memcpy(q, q0, sizeof q);
            for (i64 k = 0; k < sigma; k++) { i64 t[5]; F_state(q, t, B); memcpy(q, t, sizeof q); }
            i64 end = mu + lam;  /* through one full cycle */
            if (end < sigma + lam) end = sigma + lam;
            size_t L = (size_t)(end - sigma + 40);
            if (L > Rcap) { Rcap = L; Rbuf = realloc(Rbuf, Rcap * sizeof(i64)); }
            Rsig = resid(q, B);
            int reached = 0;
            for (i64 j = sigma; j < end + 33; j++) {
                i64 R = resid(q, B);
                Rbuf[j - sigma] = R;
                if (j < end) {
                    if (!reached && R > Rsig + 384) env_bad++;
                    if (reached && R > 504) env_bad++;
                    if (R <= 120) reached = 1;
                    if (reached && R > maxR_after) maxR_after = R;
                    if (j >= mu && R > 504) cyc_bad++;
                    i64 f[5], p[5]; F_state(q, f, B); P_state(q, p, B);
                    if (q[4] <= 0 || !eq5(f, p)) pair_bad++;
                }
                i64 t[5]; F_state(q, t, B); memcpy(q, t, sizeof q);
            }
            for (i64 j = sigma; j < end; j++) {
                i64 R = Rbuf[j - sigma];
                i64 d = Rbuf[j - sigma + 1] - R; if (d > maxjump) maxjump = d; if (d > 6) jump_bad++;
                if (R > 120) {
                    int ok = 0;
                    for (int k = 1; k <= 32; k++) if (Rbuf[j - sigma + k] <= R - 6) { ok = 1; if (k > hitmax) hitmax = k; break; }
                    if (!ok) hit_bad++;
                }
            }
        }
        printf("B=%lld s=%ld sigma=%lld T=%lld mu=%lld lam=%lld dist=%.4f Rsig=%lld maxRafter=%lld maxjump=%lld "
               "maxhit=%lld bad[pair=%ld jump=%ld hit=%ld env=%ld cyc=%ld consist=%d]\n",
               (long long)B, smp, (long long)sigma, (long long)T, (long long)mu, (long long)lam, dist6 / 6.0,
               (long long)Rsig, (long long)maxR_after, (long long)maxjump, (long long)hitmax,
               pair_bad, jump_bad, hit_bad, env_bad, cyc_bad, consist);
        fflush(stdout);
    }
    return 0;
}
