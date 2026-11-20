// Exact stage transport on the current affine lattice. See README.md for the
// renewal argument and the distinction between a sampled graph and closure.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using Z = std::int64_t;
using W = __int128_t;
using U = std::uint64_t;
using V = std::array<Z, 5>;
using WV = std::array<W, 5>;

Z narrow(W a) {
    if (a < INT64_MIN || a > INT64_MAX) throw std::overflow_error("stage integer");
    return Z(a);
}
Z mod(Z a, Z b) {
    assert(b > 0);
    Z r = a % b;
    return r < 0 ? r + b : r;
}
W absw(W a) { return a < 0 ? -a : a; }
W gcdw(W a, W b) {
    while (b) {
        const W r = a % b;
        a = b;
        b = r;
    }
    return absw(a);
}
V mix(Z a, const V &x, Z b, const V &y) {
    V r{};
    for (int i = 0; i < 5; ++i) r[i] = narrow(W(a)*x[i] + W(b)*y[i]);
    return r;
}
Z dot(const V &x, const V &y) {
    W r = 0;
    for (int i = 0; i < 5; ++i) r += W(x[i])*y[i];
    return narrow(r);
}
std::ostream &operator<<(std::ostream &out, const V &v) {
    for (Z z : v) out << ' ' << z;
    return out;
}

struct Ratio {
    W n = 0, d = 1;
    Ratio() = default;
    Ratio(W numerator, W denominator) : n(numerator), d(denominator) {
        assert(d > 0);
        const W g = gcdw(n, d);
        n /= g;
        d /= g;
    }
};
Ratio operator+(const Ratio &a, const Ratio &b) {
    const W g = gcdw(a.d, b.d);
    return {a.n*(b.d/g) + b.n*(a.d/g), a.d*(b.d/g)};
}
bool operator<(const Ratio &a, const Ratio &b) {
    return a.n*b.d < b.n*a.d;
}
std::string numstr(W n) {
    if (n == 0) return "0";
    std::string out;
    const bool negative = n < 0;
    n = absw(n);
    while (n) {
        out.push_back(char('0' + n % 10));
        n /= 10;
    }
    if (negative) out.push_back('-');
    std::reverse(out.begin(), out.end());
    return out;
}
std::ostream &operator<<(std::ostream &out, const Ratio &x) {
    return out << numstr(x.n) << '/' << numstr(x.d);
}

struct Point {
    WV n{}; // X/B as a common-denominator exact rational point.
    W den = 1;
    bool zero() const { return n == WV{}; }
    void reduce() {
        W g = den;
        for (W k : n) g = gcdw(g, k);
        den /= g;
        for (W &k : n) k /= g;
    }
};

struct Form {
    V row{};
    Z b = 0, c = 0; // row*X + b*B + c.
    Z eval(const V &x, Z base) const { return narrow(W(dot(row, x)) + W(b)*base + c); }
    W lead(const Point &x) const {
        W r = W(b)*x.den;
        for (int j = 0; j < 5; ++j) r += W(row[j])*x.n[j];
        return r;
    }
};
Form operator+(const Form &a, const Form &b) {
    return {mix(1, a.row, 1, b.row), narrow(W(a.b)+b.b), narrow(W(a.c)+b.c)};
}
Form operator-(const Form &a, const Form &b) {
    return {mix(1, a.row, -1, b.row), narrow(W(a.b)-b.b), narrow(W(a.c)-b.c)};
}
Form operator*(Z k, const Form &f) {
    return {mix(k, f.row, 0, {}), narrow(W(k)*f.b), narrow(W(k)*f.c)};
}

using State = std::array<Form, 5>;

State identity() {
    State result{};
    for (int i = 0; i < 5; ++i) result[i].row[i] = 1;
    return result;
}
V eval(const State &f, const V &x, Z base) {
    V result{};
    for (int i = 0; i < 5; ++i) result[i] = f[i].eval(x, base);
    return result;
}
V linear(const State &f, const V &v) {
    V result{};
    for (int i = 0; i < 5; ++i) result[i] = dot(f[i].row, v);
    return result;
}
Point lead(const State &f, const Point &x) {
    Point result{{}, x.den};
    for (int i = 0; i < 5; ++i) result.n[i] = f[i].lead(x);
    return result;
}

struct Lattice {
    V p{};
    std::vector<V> cols;
    bool varies(const Form &f) const {
        return std::any_of(cols.begin(), cols.end(), [&](const V &v) {
            return dot(f.row, v) != 0;
        });
    }
    // Integer column elimination. Zero columns are genuine kernel directions.
    void reduce() {
        std::size_t pivot = 0;
        for (int i = 0; i < 5 && pivot < cols.size(); ++i) {
            auto it = std::find_if(cols.begin()+pivot, cols.end(), [&](const V &v) {return v[i] != 0;});
            if (it == cols.end()) continue;
            std::swap(cols[pivot], *it);
            V &v = cols[pivot];
            for (std::size_t j = pivot+1; j < cols.size(); ++j) {
                V &w = cols[j];
                while (w[i]) {
                    v = mix(1, v, -v[i]/w[i], w);
                    std::swap(v, w);
                }
            }
            if (v[i] < 0) v = mix(-1, v, 0, {});
            for (std::size_t j = 0; j < pivot; ++j) {
                V &w = cols[j];
                w = mix(1, w, -(w[i]-mod(w[i], v[i]))/v[i], v);
            }
            p = mix(1, p, -(p[i]-mod(p[i], v[i]))/v[i], v);
            ++pivot;
        }
        cols.resize(pivot);
    }
};

Lattice root_lattice(Z base) {
    Lattice l;
    for (int i = 0; i < 5; ++i) {
        l.p[i] = mod((i-5)*(base%6), 6);
        V col{};
        col[i] = 6;
        l.cols.push_back(col);
    }
    return l;
}
Lattice push(const State &f, const Lattice &l, Z base) {
    Lattice out{eval(f, l.p, base), {}};
    for (const V &v : l.cols) out.cols.push_back(linear(f, v));
    out.reduce();
    return out;
}

struct Split {
    Lattice kernel;
    V shift{};
    Z index = 1;
    Z stride = 0;
};
Split split(const Lattice &l, const Form &f, Z s) {
    assert(s > 0);
    Lattice q = l;
    assert(!q.cols.empty());
    for (std::size_t j = 1; j < q.cols.size(); ++j) {
        V &v = q.cols[0], &w = q.cols[j];
        while (const Z value = dot(f.row, w)) {
            v = mix(1, v, -dot(f.row, v)/value, w);
            std::swap(v, w);
        }
    }
    V shift = q.cols[0];
    const Z stride = dot(f.row, shift);
    const Z index = s / std::gcd(s, stride);
    if (!stride) shift = {};
    q.cols[0] = mix(index, q.cols[0], 0, {});
    return {q, shift, index, stride};
}
V div(const V &v, Z d) {
    V out{};
    for (int i = 0; i < 5; ++i) {
        assert(v[i] % d == 0);
        out[i] = v[i]/d;
    }
    return out;
}

struct Rng {
    U state;
    U next() {
        U v = (state += 0x9e3779b97f4a7c15ULL);
        v = (v ^ (v >> 30)) * 0xbf58476d1ce4e5b9ULL;
        v = (v ^ (v >> 27)) * 0x94d049bb133111ebULL;
        return v ^ (v >> 31);
    }
    U below(U b) {
        const U min = -b % b;
        U x;
        do { x = next(); } while (x < min);
        return x % b;
    }
};

struct Ambiguous {};
struct ZeroWord {};

struct Context {
    const Point &x;
    const Lattice &l;
    Z base;
    bool nonnegative(const Form &f) const {
        const W a = f.lead(x);
        if (a) return a > 0;
        if (l.varies(f)) throw Ambiguous{};
        return f.eval(l.p, base) >= 0;
    }
};

struct Phase {
    State end;
    std::vector<Form> guards;
    bool paired;
    int height;
    U word;
};

// The labels are positions in the borrowed list, with the central position
// at r. Ties take increasing label. Each difference sort guard is a multiple
// of six on the physical lattice, so >0 is equivalent to >=6.
Phase primitive(const State &f, const Context &ctx) {
    std::array<Form, 5> q{};
    for (int i = 0; i < 5; ++i) q[i] = f[i] + Form{{}, 5-i, 0};
    int r = 5;
    std::vector<Form> guards;
    while (r && !ctx.nonnegative(q[r-1]-Form{{}, 0, 6})) {
        guards.push_back(-1*q[r-1]);
        --r;
    }
    if (!r) throw ZeroWord{};
    guards.push_back(q[r-1]-Form{{}, 0, 6});
    std::array<Form, 11> digits{};
    int j = 0;
    for (int i = 0; i < r; ++i) digits[j++] = q[i] - Form{{}, 0, i == r-1 ? 6 : 0};
    for (int i = 0; i < 11-2*r; ++i) digits[j++] = Form{{}, 6, -6};
    for (int i = r-1; i >= 0; --i) digits[j++] = Form{{}, 6, i ? -6 : 0} - q[i];
    assert(j == 11);
    std::array<int, 11> order{0,1,2,3,4,5,6,7,8,9,10};
    std::stable_sort(order.begin(), order.end(), [&](int i, int k) {
        return ctx.nonnegative(digits[i] - digits[k] - Form{{}, 0, i < k ? 0 : 6});
    });
    U word = 0;
    for (int i = 0; i < 11; ++i) {
        word = (word << 4) | order[i];
        if (i) guards.push_back(digits[order[i-1]] - digits[order[i]]
                              - Form{{}, 0, order[i-1] < order[i] ? 0 : 6});
    }
    bool paired = r == 5;
    if (paired) {
        std::array<int, 11> pos{};
        for (int i = 0; i < 11; ++i) pos[order[i]] = i;
        for (int i = 0; i < 4; ++i) {
            const int a = pos[i] > pos[6+i] ? i : 6+i;
            const int b = pos[i+1] < pos[7+i] ? i+1 : 7+i;
            const Form gap = digits[a]-digits[b];
            const bool ordered = ctx.nonnegative(gap);
            guards.push_back(ordered ? gap : -1*gap-Form{{}, 0, 6});
            paired &= ordered;
        }
    }
    State next{};
    for (int i = 0; i < 5; ++i)
        next[i] = digits[order[i]] - digits[order[10-i]] - Form{{}, 5-i, 0};
    return {next, guards, paired, r, word};
}

W energy(const Point &p) {
    W prev = 0, sum = 0;
    for (W x : p.n) {
        const W g = prev-x;
        sum += g*g;
        prev = x;
    }
    return sum+prev*prev;
}

constexpr int max_word = 120;
constexpr int max_short = 500;
constexpr int max_stage = 15;

struct Stage {
    Phase first;
    V d{};
    Form exit{};
    Z s = 0;
    std::vector<U> words;
    int failures = 0, zero = 0;
};

Stage stage(const Context &ctx) {
    State f = identity();
    const W e = energy(ctx.x);
    std::vector<Form> guards;
    Stage out{};
    int failures = 0, zero = 0;
    std::vector<U> words;
    for (int k = 1; k <= max_word; ++k) {
        Phase phase = primitive(f, ctx);
        if (k == 1) out.first = phase;
        words.push_back(phase.word);
        failures += !phase.paired;
        zero += phase.height < 5;
        guards.insert(guards.end(), phase.guards.begin(), phase.guards.end());
        f = phase.end;
        const Point y = lead(f, ctx.x);
        // Positive-inner leading subtraction majorizes the gaps. Restrict
        // this shortcut to strictly positive leading inner coordinates.
        if (ctx.x.n[4]+ctx.x.den > 0 && energy(y) < e) break;
        if (y.n != ctx.x.n) continue;
        bool translation = true;
        for (const V &col : ctx.l.cols) translation &= linear(f, col) == col;
        if (!translation) continue;
        const V d = mix(1, eval(f, ctx.l.p, ctx.base), -1, ctx.l.p);
        if (d == V{} || linear(f, d) != d) continue;

        std::vector<std::pair<Z, Form>> falling;
        for (const Form &g : guards) {
            const Z s = -dot(g.row, d);
            if (s > 0) falling.push_back({s, g});
        }
        if (falling.empty()) continue;
        auto before = [&](const auto &a, const auto &b) {
            const Form comparison = b.first*a.second - a.first*b.second;
            const W lead_gap = comparison.lead(ctx.x);
            if (lead_gap) return lead_gap < 0;
            if (ctx.l.varies(comparison)) throw Ambiguous{};
            return comparison.eval(ctx.l.p, ctx.base) < 0;
        };
        const auto [s, g] = *std::min_element(falling.begin(), falling.end(), before);
        if (g.lead(ctx.x) <= 0) continue;
        out.d = d;
        out.exit = g;
        out.s = s;
        out.words = words;
        out.failures = failures;
        out.zero = zero;
        return out;
    }
    return out;
}

Point projected(const Point &x, const Stage &st) {
    Point y{{}, x.den*st.s};
    const W h = st.exit.lead(x);
    for (int i = 0; i < 5; ++i) y.n[i] = x.n[i]*st.s + h*st.d[i];
    y.reduce();
    return y;
}

struct Edge {
    Lattice next;
    Z remainder;
};
Edge exit_layer(const Stage &st, const Split &split, Z lane, Z base) {
    Lattice l = split.kernel;
    l.p = mix(1, l.p, lane, split.shift);
    const Z a = st.exit.eval(l.p, base);
    const Z r = mod(a, st.s);
    l.p = mix(1, l.p, 1+(a-r)/st.s, st.d);
    for (V &v : l.cols) {
        v = div(mix(st.s, v, dot(st.exit.row,v), st.d), st.s);
    }
    assert(st.exit.eval(l.p, base) == r-st.s);
    l.reduce();
    return {l, r};
}

using Digits = std::array<Z, 11>;
V differences(Digits w) {
    std::sort(w.begin(), w.end());
    V q{};
    for (int i = 0; i < 5; ++i) q[i] = w[10-i]-w[i];
    return q;
}
V leading_step(const V &q, Z base) {
    Digits w{};
    w[5] = base;
    for (int i = 0; i < 5; ++i) {
        w[i] = q[i];
        w[10-i] = base-q[i];
    }
    return differences(w);
}
bool in_pi(V q, Z b) {
    for (int k = 0; k < 5; ++k) {
        if (q[4] <= 0 || q[0] >= b) return false;
        for (int i = 0; i < 5; ++i) {
            if (q[i] <= b-q[4-i]) return false;
            if (i < 4 && b-q[4-i] <= q[i+1]) return false;
        }
        q = {q[0], q[0]-q[4], q[1]-q[4], q[1]-q[3], q[2]-q[3]};
    }
    return true;
}
int entry(V &q, Z b, int cap) {
    for (int j = 0; j <= cap; ++j) {
        if (in_pi(q, b)) return j;
        const V next = leading_step(q, b);
        if (next == q) return -1;
        q = next;
    }
    return -1;
}
Point centered(const V &q, Z base) {
    Point p{{}, base};
    for (int i = 0; i < 5; ++i) p.n[i] = W(6)*q[i] - W(5-i)*base;
    p.reduce();
    return p;
}

struct Call {
    Ratio elapsed;
    int stages = 0;
    int shorts = 0;
    int first_loss = -1;
    bool any_zero = false;
    int keeps_rank = 0;
};
struct Trail {
    int status = 0; // Complete, no long return, accidental tie, stage cap, or zero word.
    Call count;
    int rank = 0;
    Ratio radius;
    V original{};
    std::vector<std::tuple<int, int, int, int, int, V, Form, Z>> stages;
};
Ratio radius(const Point &x) {
    W prev = 0, r = 0;
    for (W v : x.n) {
        r = std::max(r, v-prev);
        prev = v;
    }
    return {std::max(r,-prev), x.den};
}

Trail trace(Point x, Lattice l, Z base, Rng &rng, bool keep_log) {
    Call c{};
    Trail result;
    for (; c.stages < max_stage && !x.zero();) {
        result.rank = int(l.cols.size());
        try {
            const Stage st = stage({x, l, base});
            if (!st.s) {
                if (!st.first.paired && c.first_loss < 0) c.first_loss = c.stages;
                c.any_zero |= st.first.height < 5;
                if (++c.shorts > max_short) { result.status = 1; break; }
                x = lead(st.first.end, x);
                l = push(st.first.end, l, base);
                continue;
            }
            if (st.failures && c.first_loss < 0) c.first_loss = c.stages;
            c.any_zero |= st.zero != 0;
            if (keep_log) result.stages.emplace_back(c.stages, c.shorts, int(l.cols.size()),
                int(st.words.size()), st.failures, st.d, st.exit, st.s);
            c.elapsed = c.elapsed+Ratio{W(st.words.size())*st.exit.lead(x), x.den*st.s};
            const Split lanes = split(l, st.exit, st.s);
            const Edge edge = exit_layer(st, lanes, Z(rng.below(lanes.index)), base);
            c.keeps_rank += edge.next.cols.size() == l.cols.size();
            x = projected(x, st);
            l = edge.next;
            ++c.stages;
        } catch (const Ambiguous &) {
            result.status = 2;
            break;
        } catch (const ZeroWord &) {
            result.status = 4;
            break;
        }
    }
    if (result.status == 0 && !x.zero()) result.status = 3;
    result.count = c;
    result.radius = radius(x);
    return result;
}

int main(int argc, char **argv) {
    if (argc != 5 && argc != 11) {
        std::cerr << "usage: graph base count seed uniform-pi|entry\n"
                     "       graph base count seed point grid q1 q2 q3 q4 q5\n";
        return 2;
    }
    const Z base = std::stoll(argv[1]);
    const int count = std::stoi(argv[2]);
    Rng sources{std::stoull(argv[3])}, lanes{std::stoull(argv[3]) ^ 0xb6ac13b6111cc387ULL};
    const std::string mode = argv[4];
    const Z grid = argc == 11 ? std::stoll(argv[5]) : Z(1) << 38;
    V point{};
    if (mode == "point") for (int j = 0; j < 5; ++j) point[j] = std::stoll(argv[6+j]);
    for (int j = 0; j < count; ++j) {
        V q{};
        do {
            if (mode == "point") q = point;
            else if (mode == "uniform-pi") {
                for (Z &v : q) v = Z(sources.below(grid));
                std::sort(q.begin(), q.end(), std::greater<Z>());
            } else {
                Digits w{};
                for (Z &v : w) v = Z(sources.below(grid));
                q = differences(w);
            }
        } while (mode == "uniform-pi" && !in_pi(q, grid));
        const V original = q;
        const int h = entry(q, grid, 500);
        if (h < 0) {
            std::cout << "5 -1 -1 0 0 0 0 0 0/1 0/1" << original << '\n';
            continue;
        }
        Trail t = trace(centered(q, grid), root_lattice(base), base, lanes, mode == "point");
        const Call &c = t.count;
        std::cout << t.status << ' ' << h << ' ' << c.first_loss << ' ' << c.any_zero
                  << ' ' << c.stages << ' ' << c.shorts << ' ' << t.rank << ' ' << c.keeps_rank
                  << ' ' << c.elapsed << ' ' << t.radius << original << '\n';
        for (const auto &[id, short_steps, rank, length, failures, drift, g, s] : t.stages)
            std::cout << "stage " << id << ' ' << short_steps << ' ' << rank << ' ' << length
                      << ' ' << failures << " drift" << drift << " guard" << g.row
                      << ' ' << g.b << ' ' << g.c << " decrement " << s << '\n';
    }
}

