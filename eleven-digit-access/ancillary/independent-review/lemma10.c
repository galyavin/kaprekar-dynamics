/* Exhaustive box tests of Lemma 10 (integer recurrence eq.(paired-return)).
 * Own encoding of the displayed formulas (Appendix D).
 *  mode A: second assertion. For the 50 triples (p,b0,r), all (l,w) with -27<=l<=LMAX, -WMAX<=w<=2:
 *          no (l,w) keeps all 18 guard pairs with all 18 cumulative changes < 2.
 *  mode B: composite claim used by Lemma 5, on a 5-dim box with p<=2: within 19 returns a prefix
 *          gains >=1 in q4-q5, or a guard fails; plus the first assertion's box for states whose
 *          first three returns are guarded with nonpositive cumulative changes.
 *  mode C: composite claim on random states with large coordinates (application regime:
 *          l>=0 up to 1e6, w<=0 down to -1e6, |p|,|s|,|r| up to 1e4, p<=2).
 * Usage: lemma10 A LMAX WMAX | lemma10 B L0 L1 W0 W1 P0 S R | lemma10 C N SEED
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef long long ll;
static ll pos(ll x) { return x > 0 ? x : 0; }
static ll mx(ll a, ll b) { return a > b ? a : b; }
static ll mn(ll a, ll b) { return a < b ? a : b; }
typedef struct { ll l, w, p, s, r; } st;
static int guards(st x) {
    ll Kf = pos(x.r + 1);
    return (x.w - x.l + x.p <= -1) && (x.l - 2 * x.w - 2 * Kf >= -1);
}
static st ret(st x, ll *dmin) {
    ll Cf = pos(-x.l + 2 * x.s + 1), Af = mx(-1, -x.r);
    ll J = pos(-x.r - 1) - pos(x.r - 1), Kf = pos(x.r + 1);
    ll Hf = pos(x.w - J - x.p - 1), t = x.s - x.p - Cf;
    ll Df = mx(Af - 1, t);
    ll Ef = mx(Af, t - 1) - mx(-t - 1, -Af);
    ll Ff = mn(Af, t - 1) + J + 1;
    st y = { x.l - 6 * Df, x.w + J - x.p - 4 * Df, Ef - 2 * Df, Ff - 3 * Df + Hf, -x.p - Kf + Df - Hf };
    *dmin = x.r + J + Kf + Hf;
    return y;
}
int main(int argc, char **argv) {
    if (argv[1][0] == 'A') {
        ll LMAX = atoll(argv[2]), WMAX = atoll(argv[3]);
        long ntrip = 0, cex = 0; ll tested = 0;
        for (ll p = -2; p <= 2; p++) for (ll b0 = -4; b0 <= -1; b0++) for (ll r = -b0 - 4; r <= 0; r++) {
            ntrip++;
            for (ll l = -27; l <= LMAX; l++) for (ll w = -WMAX; w <= 2; w++) {
                st x = { l, w, p, p + b0, r };
                ll cum = 0; int ok = 0;
                for (int k = 0; k < 18; k++) {
                    if (!guards(x)) { ok = 1; break; }
                    ll d; x = ret(x, &d); cum += d;
                    if (cum >= 2) { ok = 1; break; }
                }
                tested++;
                if (!ok) { cex++; if (cex <= 5) printf("CEX p=%lld b0=%lld r=%lld l=%lld w=%lld\n", p, b0, r, l, w); }
            }
        }
        printf("mode A: triples=%ld (l,w) pairs=%lld counterexamples=%ld\n", ntrip, tested, cex);
    } else if (argv[1][0] == 'C') {
        long n = atol(argv[2]); unsigned long long z = strtoull(argv[3], 0, 10);
        long guarded = 0, cex = 0; int maxret = 0;
        for (long i = 0; i < n; i++) {
            ll v[5];
            for (int k = 0; k < 5; k++) { z = z * 6364136223846793005ULL + 1442695040888963407ULL; v[k] = (ll)(z >> 33); }
            ll sc = (i % 3 == 0) ? 100 : (i % 3 == 1) ? 10000 : 1000000;
            ll ps_ = sc < 10000 ? sc : 10000;
            st x = { v[0] % (sc + 1), -(v[1] % (sc + 1)), 2 - v[2] % (ps_ + 1), v[3] % (2 * ps_ + 1) - ps_, v[4] % (2 * ps_ + 1) - ps_ };
            if (!guards(x)) continue;
            guarded++;
            ll cum = 0; int ok = 0;
            for (int k = 0; k < 19; k++) {
                if (!guards(x)) { ok = 1; break; }
                ll d; x = ret(x, &d); cum += d;
                if (cum >= 1) { ok = 1; if (k + 1 > maxret) maxret = k + 1; break; }
            }
            if (!ok) cex++;
        }
        printf("mode C: random states=%ld guarded=%ld counterexamples=%ld max returns to first gain=%d\n", n, guarded, cex, maxret);
    } else {
        ll L0 = atoll(argv[2]), L1 = atoll(argv[3]), W0 = atoll(argv[4]), W1 = atoll(argv[5]);
        ll P0 = atoll(argv[6]), S = atoll(argv[7]), R = atoll(argv[8]);
        ll tested = 0, cex = 0, boxviol = 0, firstbox = 0, maxret = 0;
        for (ll p = P0; p <= 2; p++) for (ll s = -S; s <= S; s++) for (ll r = -R; r <= R; r++)
        for (ll l = L0; l <= L1; l++) for (ll w = W0; w <= W1; w++) {
            st x0 = { l, w, p, s, r }, x = x0, x1 = x0;
            ll cum = 0; int ok = 0, k; int cumnp3 = 1, guard3 = 1;
            for (k = 0; k < 19; k++) {
                if (!guards(x)) { ok = 1; if (k < 3) guard3 = 0; break; }
                ll d; x = ret(x, &d); cum += d;
                if (k == 0) x1 = x;
                if (k < 3 && cum > 0) cumnp3 = 0;
                if (cum >= 1) { ok = 1; if (k + 1 > maxret) maxret = k + 1; break; }
            }
            tested++;
            if (!ok) { cex++; if (cex <= 5) printf("CEX l=%lld w=%lld p=%lld s=%lld r=%lld\n", l, w, p, s, r); }
            /* first assertion: if first three returns guarded with cumulative <= 0, x1 in box */
            if (guard3 && cumnp3 && k >= 3) {
                firstbox++;
                ll b0 = x1.s - x1.p;
                if (!(x1.l >= -27 && x1.w <= 2 && x1.p >= -2 && x1.p <= 2 && b0 >= -4 && b0 <= -1 && x1.r >= -b0 - 4 && x1.r <= 0)) {
                    boxviol++;
                    if (boxviol <= 5) printf("BOXVIOL x0=(%lld,%lld,%lld,%lld,%lld) x1=(%lld,%lld,%lld,%lld,%lld)\n", l, w, p, s, r, x1.l, x1.w, x1.p, x1.s, x1.r);
                }
            }
        }
        printf("mode B: states=%lld no-gain-no-exit counterexamples=%lld max returns to first gain=%lld "
               "first-assertion cases=%lld box violations=%lld\n", tested, cex, maxret, firstbox, boxviol);
    }
    return 0;
}
