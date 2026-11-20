// Independent actual clock: sort eleven digits and subtract column by column.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using Word = std::array<int, 11>;
using Q = std::array<int, 5>;
using U = std::uint64_t;
using Wide = __uint128_t;

struct Rng {
    U state;
    U next() {
        U v = (state += 0x9e3779b97f4a7c15ULL);
        v = (v ^ (v >> 30)) * 0xbf58476d1ce4e5b9ULL;
        v = (v ^ (v >> 27)) * 0x94d049bb133111ebULL;
        return v ^ (v >> 31);
    }
    int digit(U b) {
        const U bound = -b % b;
        Wide x;
        do { x = Wide(next())*b; } while (U(x) < bound);
        return int(x >> 64);
    }
};

Word sorted(Word w) {
    for (std::size_t i = 1; i < w.size(); ++i) {
        int v = w[i];
        std::size_t j = i;
        while (j && w[j-1] < v) {
            w[j] = w[j-1];
            --j;
        }
        w[j] = v;
    }
    return w;
}

Q differences(const Word &w) {
    Q q{};
    for (int i = 0; i < 5; ++i) q[i] = w[i]-w[10-i];
    return q;
}

Word update(const Word &ordered, int base) {
    Word result{};
    int borrow = 0;
    for (int i = 10; i >= 0; --i) {
        int v = ordered[i]-ordered[10-i]-borrow;
        borrow = v < 0;
        result[i] = v+borrow*base;
    }
    assert(borrow == 0);
    return result;
}

bool leading_pi(Q q, int b) {
    for (int t = 0; t < 5; ++t) {
        if (!(q[4] > 0 && q[0] < b)) return false;
        for (int i = 0; i < 5; ++i) {
            if (q[i]+q[4-i] <= b) return false;
            if (i < 4 && b-q[4-i] <= q[i+1]) return false;
        }
        q = {q[0], q[0]-q[4], q[1]-q[4], q[1]-q[3], q[2]-q[3]};
    }
    return true;
}

struct Pairing {
    bool regular;
    int loss; // 0 paired, 1..4 interlevel inversion, 5 zero innermost.
    int slack;
};
Pairing pairing(const Q &q, int b) {
    if (!q[4]) return {false, 5, 0};
    bool regular = true;
    int loss = 0, slack = 0;
    for (int i = 0; i < 5; ++i) {
        const int a = q[i] - (i == 4), c = b-q[4-i]-(i < 4);
        regular &= a >= c;
        if (i < 4) {
            const int upper = std::max(q[i+1]-(i == 3), b-q[3-i]-(i < 3));
            const int s = std::min(a,c)-upper;
            if (s < slack) { loss = i+1; slack = s; }
        }
    }
    return {regular && !loss, loss, slack};
}

void observe(const Word &w0, int base, std::vector<Word> &path) {
    int checkpoint = 0, power = 1, period = 0;
    int run = 0, entry = -1;
    int first_loss = -1, last_loss = -1, loss_count = 0, mechanism = 0, slack = 0;
    Q at_loss{}, original = differences(sorted(w0));
    Word w = w0;
    path.clear();
    for (int t = 0; t < 48*base+10000; ++t) {
        if (t && w == path[checkpoint]) {
            period = t-checkpoint;
            break;
        }
        if (t-checkpoint == power) {
            checkpoint = t;
            power *= 2;
        }
        path.push_back(w);
        const Word s = sorted(w);
        const Q q = differences(s);
        if (entry < 0) {
            run = pairing(q, base).regular ? run+1 : 0;
            if (run == 5) entry = t-4;
        } else {
            const Pairing p = pairing(q, base);
            if (p.loss) {
                ++loss_count;
                last_loss = t;
                if (first_loss < 0) {
                    first_loss = t;
                    mechanism = p.loss;
                    slack = p.slack;
                    at_loss = q;
                }
            }
        }
        w = update(s, base);
    }

    int capture = -1;
    Q terminal{};
    if (period) {
        int lo = 0, hi = checkpoint;
        while (lo < hi) {
            const int mid = (lo+hi)/2;
            if (path[mid] == path[mid+period]) hi = mid;
            else lo = mid+1;
        }
        capture = lo;
        terminal = differences(sorted(path[capture]));
    }
    std::cout << capture << ' ' << period << ' ' << entry << ' ' << first_loss << ' '
              << last_loss << ' ' << loss_count << ' ' << mechanism << ' ' << slack;
    for (int x : original) std::cout << ' ' << x;
    for (int x : terminal) std::cout << ' ' << x;
    for (int x : at_loss) std::cout << ' ' << x;
    std::cout << '\n';
}

Word normalized_box(int b, Rng &r) {
    // Rejection in a rectangular superset in q coordinates. Exact integer
    // tests for the five defining chart inequalities follow each draw.
    constexpr Q lo{94040,73770,57940,41685,23730};
    constexpr Q hi{94130,73850,58015,41860,23845};
    Q lower{}, upper{}, q{};
    for (int i = 0; i < 5; ++i) {
        lower[i] = int(std::int64_t(lo[i])*b/100000);
        upper[i] = int(std::int64_t(hi[i])*b/100000)+1;
    }
    for (;;) {
        for (int i = 0; i < 5; ++i) q[i] = lower[i]+r.digit(upper[i]-lower[i]+1);
        const Q x{q[0],q[1],q[2],2*(q[1]-q[4])-b,b-q[2]-q[3]};
        constexpr Q left{94040,73770,57940,20,200}, right{94130,73850,58015,80,300};
        bool inside = true;
        for (int i = 0; i < 5; ++i)
            inside &= std::int64_t(left[i])*b <= std::int64_t(x[i])*100000 &&
                      std::int64_t(x[i])*100000 <= std::int64_t(right[i])*b;
        if (!inside) continue;
        // Original limiting density by exact rejection against a product
        // bound for the six gap factors. The integer test uses 128 bits.
        const Q previous{b,q[0],q[1],q[2],q[3]};
        const Q previous_hi{b,upper[0],upper[1],upper[2],upper[3]};
        Wide num = q[4], den = upper[4];
        for (int i = 0; i < 5; ++i) {
            num *= previous[i]-q[i];
            den *= previous_hi[i]-lower[i];
        }
        const U scale = U(1) << 40;
        if (Wide(r.next() % scale)*den >= scale*num) continue;
        Word w{};
        std::copy(q.begin(), q.end(), w.begin());
        return w;
    }
}
int main(int argc, char **argv) {
    if (argc != 5) {
        std::cerr << "usage: actual base count seed original|uniform-pi|box\n";
        return 2;
    }
    const int base = std::stoi(argv[1]);
    const int count = std::stoi(argv[2]);
    Rng rng{std::stoull(argv[3])};
    const std::string mode = argv[4];
    std::vector<Word> path;
    for (int row = 0; row < count; ++row) {
        Word w{};
        if (mode == "box") w = normalized_box(base, rng);
        else {
            do {
                w.fill(0);
                if (mode == "uniform-pi") {
                    for (int i = 0; i < 5; ++i) w[i] = rng.digit(base);
                    std::sort(w.begin(), w.begin()+5, std::greater<int>());
                } else {
                    for (int &v : w) v = rng.digit(base);
                }
                const Q q = differences(sorted(w));
                if (q[0] && (mode != "uniform-pi" || leading_pi(q, base))) break;
            } while (true);
        }
        observe(w, base, path);
    }
}

