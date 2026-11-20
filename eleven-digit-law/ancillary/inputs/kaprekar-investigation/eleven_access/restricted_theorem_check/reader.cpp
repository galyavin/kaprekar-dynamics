#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <deque>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

using Z = __int128_t;
using Vec = std::array<Z, 5>;
using Form = std::array<Z, 7>;
using Linear = std::array<Z, 6>;
using Map = std::array<Form, 5>;

Z gcd(Z a, Z b) {
    while (b) {
        const auto remainder = a % b;
        a = b;
        b = remainder;
    }
    return a < 0 ? -a : a;
}

template <class T> T mix(Z a, const T& x, Z b, const T& y) {
    T output{};
    for (size_t i = 0; i < output.size(); ++i) output[i] = a*x[i] + b*y[i];
    return output;
}

template <class T, class U> Z dot(const T& x, const U& y) {
    return std::inner_product(y.begin(), y.end(), x.begin(), Z(0));
}

template <class T> void simplify(T& row) {
    auto common = std::accumulate(row.begin(), row.end(), Z(0), gcd);
    if (common > 1) for (auto& entry : row) entry /= common;
}

Z read_int() {
    std::string s;
    std::cin >> s;
    Z value = 0;
    size_t j = !s.empty() && s[0] == '-';
    for (size_t k = j; k < s.size(); ++k) value = value * 10 + (s[k] - '0');
    return j ? -value : value;
}

void write_int(Z n) {
    if (n < 0) { std::cout.put('-'); n = -n; }
    char digits[50];
    size_t length = 0;
    do { digits[length++] = '0' + n % 10; n /= 10; } while (n);
    while (length) std::cout.put(digits[--length]);
}

template <class T> void write_row(const T& r) {
    for (auto v : r) { write_int(v); std::cout.put(' '); }
    std::cout.put('\n');
}

constexpr Form UNIT = {0, 0, 0, 0, 0, 0, 1};
constexpr Form BASE = {0, 0, 0, 0, 0, 1, 0};
constexpr Map identity() {
    Map m{};
    for (size_t i = 0; i < m.size(); ++i) m[i][i] = 1;
    return m;
}

struct Sheet {
    Form origin;
    Linear shift;
    std::vector<Vec> basis;

    void advance(const Map& m) {
        for (auto& v : basis) {
            const auto before = v;
            for (size_t i = 0; i < 5; ++i) v[i] = dot(m[i], before);
        }
        const auto before = origin;
        const auto earlier = shift;
        for (size_t i = 0; i < 5; ++i) {
            origin[i] = dot(m[i], before);
            shift[i] = dot(m[i], earlier);
        }
    }
};

struct Shape {
    std::array<Linear, 5> source;
    Z denominator;
    Linear point;

    Linear pull(const Form& f) const {
        Linear row{};
        row[5] = f[5] * denominator;
        for (size_t i = 0; i < 5; ++i) row = mix(1, row, f[i], source[i]);
        return row;
    }

    void advance(const Map& word) {
        std::array<Linear, 5> rows;
        for (size_t i = 0; i < 5; ++i) rows[i] = pull(word[i]);
        source = rows;
    }

    Linear position() const {
        Linear x{};
        x[5] = point[5] * denominator;
        for (size_t i = 0; i < 5; ++i) x[i] = dot(source[i], point);
        simplify(x);
        return x;
    }
};

struct Phase {
    Map image;
    std::vector<Form> guards;
    bool paired;
    std::array<int, 11> order;
};

Phase step(const Map& current, const Linear& lead, const Form& exact) {
    std::array<Form, 11> digits{};
    for (int i = 0; i < 5; ++i) {
        digits[i] = mix(1, current[i], -(i == 4), UNIT);
        digits[10-i] = mix(-1, current[i], 1, mix(1, BASE, -(i > 0), UNIT));
    }
    digits[5] = mix(1, BASE, -1, UNIT);
    std::array<int, 11> order{};
    std::iota(order.begin(), order.end(), 0);
    auto sign = [&](const Form& f) {
        Z v = dot(f, lead);
        if (!v) v = dot(f, exact);
        return (v > 0) - (v < 0);
    };
    std::sort(order.begin(), order.end(), [&](int i, int j) {
        const auto z = sign(mix(1, digits[i], -1, digits[j]));
        return z ? z > 0 : i < j;
    });
    std::array<int, 11> pos;
    for (int i = 0; i < 11; ++i) pos[order[i]] = i;
    Phase result{{}, {mix(1, current[4], -1, UNIT)}, true, order};
    for (int i = 0; i < 5; ++i) result.image[i] = mix(1, digits[order[i]], -1, digits[order[10-i]]);
    for (int i = 0; i < 10; ++i) {
        int a = order[i], b = order[i+1];
        result.guards.push_back(mix(1, mix(1, digits[a], -1, digits[b]), -(a > b), UNIT));
    }
    for (int i = 0; i < 4; ++i) {
        int upper = pos[i] > pos[6+i] ? i : 6+i;
        int lower = pos[i+1] < pos[7+i] ? i+1 : 7+i;
        Form guard = mix(1, digits[upper], -1, digits[lower]);
        if (sign(guard) < 0) {
            result.paired = false;
            guard = mix(-1, guard, -1, UNIT);
        }
        result.guards.push_back(guard);
    }
    return result;
}

Z energy(const Linear& x) {
    Z value = x[4]*x[4];
    for (size_t i = 0; i < 5; ++i) {
        const Z gap = (i ? x[i-1] : x[5]) - x[i];
        value += gap*gap;
    }
    return value;
}

struct Exit {
    Form guard;
    Z step;
};

using RForm = std::array<Z, 7>;

void summarize(const std::set<Linear>& guards) {
    std::cout << guards.size() << '\n';
    for (const auto& g : guards) write_row(g);
    std::cout.flush();
}

// A positive macroscopic word must be visibly present twice before its first exit.
// Store the primitive path once so rejecting a period at successive short steps
// need not reconstruct the same two copies of every candidate.
class LookAhead {
    struct Frame {
        Shape shape;
        Sheet sheet;
        Linear point;
        Phase phase;
    };
    std::deque<Frame> frames;
    size_t verified = 0;

public:
    bool excludes_period(const Shape& start, const Sheet& input, size_t horizon,
                         std::set<Linear>& guards) {
        while (frames.size() < 2*horizon+1) {
            Shape shape = frames.empty() ? start : frames.back().shape;
            Sheet sheet = frames.empty() ? input : frames.back().sheet;
            if (!frames.empty()) {
                shape.advance(frames.back().phase.image);
                sheet.advance(frames.back().phase.image);
            }
            const auto x = shape.position();
            auto phase = step(identity(), x, sheet.origin);
            frames.push_back(Frame{shape, std::move(sheet), x, std::move(phase)});
        }
        for (size_t period = 1; period <= horizon; ++period) {
            if (frames[period].point != frames[0].point) continue;
            bool repeated = true;
            for (size_t k = 0; k < period && repeated; ++k) {
                const auto& a = frames[k].phase;
                const auto& b = frames[k+period].phase;
                repeated = a.paired == b.paired && a.order == b.order;
            }
            if (repeated) return false;
        }
        for (; verified < frames.size(); ++verified) {
            const auto& f = frames[verified];
            for (const auto& guard : f.phase.guards) {
                auto row = f.shape.pull(guard);
                if (row == Linear{}) {
                    assert(dot(guard, f.sheet.origin) >= 0);
                    assert(dot(guard, f.sheet.shift) == 0);
                    for (const auto& v : f.sheet.basis) assert(dot(guard, v) == 0);
                } else {
                    if (dot(guard, f.point) <= 0) {
                        frames.clear();
                        verified = 0;
                        return false;
                    }
                    simplify(row);
                    guards.insert(row);
                }
            }
        }
        return true;
    }

    Map first() const { return frames.front().phase.image; }
    void advance() {
        if (frames.empty()) return;
        frames.pop_front();
        if (verified) --verified;
    }
};

void read_chart(Shape shape, Sheet sheet, int horizon, int budget) {
    std::set<Linear> guards;
    LookAhead lookahead;
    for (int shorts = 0; shorts <= budget; ++shorts) {
        const auto x = shape.position();
        const auto initial_energy = energy(x);
        if (6 * initial_energy == x[5]*x[5]) {
            const auto normalized = shape.source;
            bool central = true;
            for (size_t i = 0; i < 5; ++i) {
                const Linear expected{0, 0, 0, 0, 0, (5-Z(i))*shape.denominator};
                central &= mix(6, normalized[i], -1, expected) == Linear{};
            }
            std::cout << (central ? "central\n" : "wall\n");
            summarize(guards);
            return;
        }
        if (shorts >= 4 && budget == 250 && sheet.basis.size() == 1 &&
            lookahead.excludes_period(shape, sheet, horizon, guards)) {
            const auto next = lookahead.first();
            shape.advance(next);
            sheet.advance(next);
            lookahead.advance();
            continue;
        }
        Map current = identity();
        Map first{};
        std::set<Form> primitive_guards;
        int failures = 0;
        for (int n = 1; n <= horizon; ++n) {
            Phase phase = step(current, x, sheet.origin);
            current = phase.image;
            if (n == 1) first = current;
            failures += !phase.paired;
            for (const auto& g : phase.guards) {
                if (!primitive_guards.insert(g).second) continue;
                auto projected = shape.pull(g);
                if (projected == Linear{}) {
                    assert(dot(g, sheet.origin) >= 0);
                    assert(dot(g, sheet.shift) == 0);
                    for (const auto& v : sheet.basis) assert(dot(g, v) == 0);
                } else {
                    if (dot(g, x) <= 0) {
                        std::cout << "wall\n";
                        summarize(guards);
                        return;
                    }
                    simplify(projected);
                    guards.insert(projected);
                }
            }
            Linear y = x;
            for (size_t i = 0; i < 5; ++i) y[i] = dot(current[i], x);
            if (energy(y) < initial_energy) break;
            if (y != x) continue;
            Vec drift{};
            bool identity_return = true;
            for (size_t i = 0; i < 5; ++i) {
                identity_return &= (shape.pull(current[i]) == shape.source[i]);
                drift[i] = dot(current[i], sheet.origin) - sheet.origin[i];
            }
            if (!identity_return || drift == Vec{}) continue;
            for (size_t i = 0; i < 5; ++i) identity_return &= dot(current[i], drift) == drift[i];
            if (!identity_return) continue;
            for (const auto& c : sheet.basis)
                for (size_t i = 0; i < 5; ++i) assert(dot(current[i], c) == c[i]);
            for (size_t i = 0; i < 5; ++i) assert(dot(current[i], sheet.shift) == sheet.shift[i]);
            std::map<RForm, Exit> exits;
            bool positive = true;
            for (const auto& g : primitive_guards) {
                Z speed = -dot(g, drift);
                if (speed <= 0) continue;
                if (dot(g, x) == 0) { positive = false; break; }
                const auto f = shape.pull(g);
                RForm key{};
                std::copy(f.begin(), f.end(), key.begin());
                key[6] = speed * shape.denominator;
                simplify(key);
                auto [it, fresh] = exits.emplace(key, Exit{g, speed});
                if (!fresh) {
                    const auto candidate = dot(g, sheet.origin) * it->second.step;
                    const auto old = dot(it->second.guard, sheet.origin) * speed;
                    if (candidate < old) it->second = {g, speed};
                }
            }
            if (!positive || exits.empty()) continue;
            std::cout << "stage " << shorts << ' ' << n << ' ' << failures << '\n';
            write_row(drift);
            for (const auto& f : shape.source) write_row(f);
            write_int(shape.denominator); std::cout.put('\n');
            write_row(sheet.origin);
            write_row(sheet.shift);
            std::cout << sheet.basis.size() << '\n';
            for (const auto& c : sheet.basis) write_row(c);
            std::cout << exits.size() << '\n';
            for (const auto& [form, exit] : exits) {
                write_row(form);
                write_row(exit.guard);
                write_int(exit.step); std::cout.put('\n');
            }
            summarize(guards);
            return;
        }
        shape.advance(first);
        sheet.advance(first);
        lookahead.advance();
    }
    std::cout << "exhausted\n";
    summarize(guards);
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    int horizon, budget, rank;
    while (std::cin >> horizon >> budget >> rank) {
        Shape shape;
        Sheet sheet;
        shape.denominator = read_int();
        for (auto& row : shape.source) for (auto& value : row) value = read_int();
        for (auto& value : shape.point) value = read_int();
        for (auto& value : sheet.origin) value = read_int();
        for (auto& value : sheet.shift) value = read_int();
        sheet.basis.resize(rank);
        for (auto& row : sheet.basis) for (auto& value : row) value = read_int();
        read_chart(shape, sheet, horizon, budget);
    }
}
