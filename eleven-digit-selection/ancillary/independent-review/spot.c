/* Post-read helper (written after the pre-read hashes were recorded).
 * Centered paired-level map, same primitive as cen.c.
 * Input lines: B G0 G1 G2 G3 G4 S1 S2 S3 S4 S5
 * Output: B G0..G4 HIT MU T  where HIT = first t with X_t == S (or -1 if not
 * reached before the first cyclic state), MU = first cyclic difference time,
 * T = MU + 1.
 */
#include <stdio.h>
#include <string.h>

typedef long long ll;

static void step(const ll *X, ll *Y) {
    ll L[5] = {X[0], X[1], X[2], X[3], X[4] - 6};
    ll R[5] = {-X[4] - 6, -X[3] - 6, -X[2] - 6, -X[1] - 6, -X[0]};
    ll hi[5], lo[5];
    for (int k = 0; k < 5; k++) {
        if (L[k] >= R[k]) { hi[k] = L[k]; lo[k] = R[k]; } else { hi[k] = R[k]; lo[k] = L[k]; }
    }
    Y[0] = -6 - lo[4]; Y[1] = hi[0] - hi[4]; Y[2] = lo[0] - lo[3];
    Y[3] = hi[1] - hi[3]; Y[4] = lo[1] - lo[2];
}

static ll resid(const ll *X) {
    ll g[6] = {-X[0], X[0] - X[1], X[1] - X[2], X[2] - X[3], X[3] - X[4], X[4]};
    ll m = g[0];
    for (int i = 1; i < 6; i++) if (g[i] < m) m = g[i];
    return -m;
}

static int eq(const ll *a, const ll *b) { return memcmp(a, b, 5 * sizeof(ll)) == 0; }

int main(void) {
    ll B, G[5], S[5];
    while (scanf("%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld", &B, &G[0], &G[1], &G[2], &G[3], &G[4],
                 &S[0], &S[1], &S[2], &S[3], &S[4]) == 11) {
        ll X[5], Y[5], s = 0;
        for (int i = 0; i < 5; i++) { s += G[i]; X[i] = (ll)(i + 1) * B - 6 * s; }
        ll t = 0, hit = -1;
        while (resid(X) > 504) {
            if (hit < 0 && eq(X, S)) hit = t;
            step(X, Y); memcpy(X, Y, sizeof X); t++;
        }
        ll t0 = t;
        ll tort[5], hare[5], tmp[5], power = 1, lam = 1;
        memcpy(tort, X, sizeof tort); step(X, hare);
        while (!eq(tort, hare)) {
            if (power == lam) { memcpy(tort, hare, sizeof tort); power *= 2; lam = 0; }
            step(hare, tmp); memcpy(hare, tmp, sizeof hare); lam++;
        }
        memcpy(tort, X, sizeof tort); memcpy(hare, X, sizeof hare);
        for (ll i = 0; i < lam; i++) { step(hare, tmp); memcpy(hare, tmp, sizeof hare); }
        ll mu = 0;
        while (!eq(tort, hare)) {
            if (hit < 0 && eq(tort, S)) hit = t0 + mu;
            step(tort, tmp); memcpy(tort, tmp, sizeof tort);
            step(hare, tmp); memcpy(hare, tmp, sizeof hare); mu++;
        }
        if (mu == 0) { printf("FAIL\n"); fflush(stdout); continue; }
        printf("%lld %lld %lld %lld %lld %lld %lld %lld %lld\n", B, G[0], G[1], G[2], G[3], G[4], hit, t0 + mu, t0 + mu + 1);
        fflush(stdout);
    }
    return 0;
}
