// Exact guarded route and population certificate on the small gap rectangle.
// Build: clang++ -O3 -std=c++17 -Wall -Wextra verify.cpp -o verify
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using Z = int64_t;
using Wide = __int128_t;
using V = std::array<Z, 5>;
using Basis = std::array<V, 5>;
using Word = std::vector<int>;
constexpr V center_gaps{26293, 27557, 28481, 40447, 29865};
constexpr Z center_base = 180000, scale = 9000000;
constexpr int short_limit = 1000, word_limit = 65, block_limit = 12;
constexpr int early_returns = 30, cycle_envelope = 504;
const V center = [] {
    V x{};
    Z sum = 0;
    for (int i = 0; i < 5; ++i) {
        sum += center_gaps[i];
        x[i] = (i + 1) * center_base - 6 * sum;
    }
    return x;
}();

void require(bool condition, const std::string &reason) {
    if (!condition) throw std::runtime_error(reason);
}

Z narrow(Wide value) {
    if (value < -Wide(std::numeric_limits<Z>::max())
        || value > std::numeric_limits<Z>::max())
        throw std::overflow_error("integer overflow");
    return static_cast<Z>(value);
}

Z dot(const V &a, const V &b) {
    Wide sum = 0;
    for (int i = 0; i < 5; ++i) sum += Wide(a[i]) * b[i];
    return narrow(sum);
}

V mix(Z a, const V &x, Z b, const V &y) {
    V out{};
    for (int i = 0; i < 5; ++i) out[i] = narrow(Wide(a) * x[i] + Wide(b) * y[i]);
    return out;
}

std::ostream &operator<<(std::ostream &out, const V &v) {
    for (Z z : v) out << ' ' << z;
    return out;
}

struct Form {
    V row{};
    Z c = 0;
    Z den = 1;
    Form() = default;
    Form(Z constant) : c(constant) {}
    Form(V a, Z b, Z d) : row(a), c(b), den(d) {
        require(d > 0, "invalid denominator");
        Z g = std::gcd(c, den);
        for (Z n : row) g = std::gcd(g, n);
        if (g == 1) return;
        for (Z &n : row) n /= g;
        c /= g;
        den /= g;
    }
    bool small() const { return row == V{}; }
    Form leading() const { return Form(row, 0, den); }
    Z integer(const V &x) const {
        Z n = narrow(Wide(dot(row, x)) + c);
        require(n % den == 0, "noninteger affine evaluation");
        return n / den;
    }
    auto key() const { return std::tie(row, c, den); }
    bool operator<(const Form &b) const { return key() < b.key(); }
    bool operator==(const Form &b) const { return key() == b.key(); }
};

std::ostream &operator<<(std::ostream &out, const Form &f) {
    return out << f.row << ' ' << f.c << ' ' << f.den;
}

Form operator-(const Form &f) { return Form(mix(-1, f.row, 0, {}), -f.c, f.den); }

Form operator+(const Form &a, const Form &b) {
    Z gcd = std::gcd(a.den, b.den);
    Z p = b.den / gcd, q = a.den / gcd;
    return Form(mix(p, a.row, q, b.row), narrow(Wide(p) * a.c + Wide(q) * b.c),
                narrow(Wide(a.den) * p));
}

Form operator-(const Form &a, const Form &b) { return a + -b; }
Form operator*(Z n, const Form &f) {
    return Form(mix(n, f.row, 0, {}), narrow(Wide(n) * f.c), f.den);
}
Form operator/(const Form &f, Z n) {
    return Form(f.row, f.c, narrow(Wide(f.den) * n));
}

using State = std::array<Form, 5>;
bool constant(const State &x) {
    return std::all_of(x.begin(), x.end(), [](const Form &f) { return f.small(); });
}
V plain(const State &x) {
    require(constant(x), "expected a constant state");
    V result{};
    for (int i = 0; i < 5; ++i) result[i] = x[i].integer({});
    return result;
}

template <typename T> using Coordinates = std::array<T, 5>;

template <typename T>
std::array<std::array<T,2>,5> pairs(const Coordinates<T> &x) {
    return {{{x[0], -x[4] - 6}, {x[1], -x[3] - 6}, {x[2], -x[2] - 6},
             {x[3], -x[1] - 6}, {x[4] - 6, -x[0]}}};
}

template <typename T>
Coordinates<T> step(const Coordinates<T> &x, int guard) {
    auto levels = pairs(x);
    std::array<T, 11> ranked{};
    ranked[0] = -6;
    for (int i = 0; i < 5; ++i) {
        bool left_first = guard & (1 << i);
        ranked[2*i+1] = levels[i][!left_first];
        ranked[2*i+2] = levels[i][left_first];
    }
    Coordinates<T> out{};
    for (int i = 0; i < 5; ++i) out[i] = ranked[i] - ranked[10-i];
    return out;
}

int integer_guard(const V &x) {
    auto levels = pairs(x);
    int guard = 0;
    for (int i = 0; i < 5; ++i) guard |= (levels[i][0] >= levels[i][1]) << i;
    return guard;
}

template <typename T>
Coordinates<T> slacks(const Coordinates<T> &x, int guard) {
    auto levels = pairs(x);
    Coordinates<T> result{};
    for (int i = 0; i < 5; ++i) {
        T d = levels[i][0] - levels[i][1];
        result[i] = guard & (1 << i) ? d : -d - 6;
    }
    return result;
}

Z threshold = 1;

bool nonnegative(const Form &f) {
    if (f.small()) return f.c >= 0;
    Z radius = 0, suffix = 0;
    for (int i = 4; i >= 0; --i) {
        suffix = narrow(Wide(suffix) + f.row[i]);
        radius = narrow(Wide(radius) + Wide(54) * std::abs(suffix));
    }
    Z middle = narrow(Wide(50) * dot(f.row, center));
    Z margin = narrow(Wide(std::abs(middle)) - radius);
    if (margin <= 0) {
        std::cerr << "undecided affine form" << f << " scaled interval "
                  << narrow(Wide(middle) - radius) << ' ' << narrow(Wide(middle) + radius) << '\n';
        throw std::runtime_error("leading guard meets cell");
    }
    Wide top = Wide(scale) * (Wide(std::abs(f.c)) + f.den);
    threshold = std::max(threshold, narrow((top + margin - 1) / margin));
    return middle > 0;
}

int choose(const State &state) {
    auto levels = pairs(state);
    int result = 0;
    for (int i = 0; i < 5; ++i) result |= nonnegative(levels[i][0] - levels[i][1]) << i;
    return result;
}

struct Block {
    Word word;
    V drift{};
    Z divisor = 0;
    Form numerator{};
};

bool center_order(const Form &a, const Form &b) {
    Wide aa = Wide(dot(a.row, center)) * b.den;
    Wide bb = Wide(dot(b.row, center)) * a.den;
    if (aa != bb) return aa < bb;
    return Wide(a.c) * b.den < Wide(b.c) * a.den;
}

std::map<State, Block> block_cache;
const Block &macro(const State &start) {
    auto old = block_cache.find(start);
    if (old != block_cache.end()) return old->second;
    State current = start;
    Word word;
    std::vector<Form> pulled;
    for (int length = 1; length <= word_limit; ++length) {
        int code = choose(current);
        word.push_back(code);
        auto gs = slacks(current, code);
        pulled.insert(pulled.end(), gs.begin(), gs.end());
        current = step(current, code);
        State diff{};
        for (int i = 0; i < 5; ++i) diff[i] = current[i] - start[i];
        if (!constant(diff)) continue;
        V drift = plain(diff), delta = drift;
        require(drift != V{}, "cycle with nonzero leading state");
        std::vector<std::pair<Z, Form>> falling;
        int k = 0;
        for (int guard : word) {
            V slopes = mix(1, slacks(delta, guard), -1, slacks(V{}, guard));
            for (Z slope : slopes) {
                if (slope < 0) falling.emplace_back(-slope, pulled[k]);
                ++k;
            }
            delta = mix(1, step(delta, guard), -1, step(V{}, guard));
        }
        if (delta != drift) continue;
        require(!falling.empty(), "no declining guard for translation");
        auto limiting = std::min_element(falling.begin(), falling.end(), [](const auto &a, const auto &b) {
            return center_order(a.second/a.first, b.second/b.first);
        });
        const auto [divisor, numerator] = *limiting;
        for (const auto &[div, form] : falling)
            require(nonnegative(form/div - numerator/divisor), "earliest exit changes on cell");
        if (numerator.small()) continue;
        return block_cache.emplace(start, Block{word, drift, divisor, numerator}).first->second;
    }
    return block_cache.emplace(start, Block{{word.front()}}).first->second;
}

Z mod(Z x, Z m) { return (x % m + m) % m; }

Z bezout(Z a, Z b, Z &s, Z &t) {
    if (!b) { s = 1; t = 0; return a; }
    Z x, y;
    Z g = bezout(b, a % b, x, y);
    s = y;
    t = narrow(Wide(x) - Wide(a/b)*y);
    return g;
}

struct Split {
    Basis basis;
    V shift;
    Z index, g, residue, inverse;
};

Split split(const Form &num, Z m, const V &point, Basis basis) {
    Z residue = mod(num.integer(point), m);
    V v{};
    for (int j = 0; j < 5; ++j) v[j] = mod(num.leading().integer(basis[j]), m);
    for (int j = 1; j < 5; ++j) {
        if (!v[j]) continue;
        Z s, t, g = bezout(v[0], v[j], s, t);
        V a = basis[0], b = basis[j];
        basis[0] = mix(s, a, t, b);
        basis[j] = mix(-v[j]/g, a, v[0]/g, b);
        v[0] = g;
    }
    Z g = std::gcd(v[0], m), index = m/g, s, t;
    require(bezout(v[0]/g, index, s, t) == 1, "modular inverse");
    V shift = basis[0];
    basis[0] = mix(index, shift, 0, {});
    return {basis, shift, index, g, residue, mod(s, index)};
}

struct Speed {
    Z period = 0;
    Z a = 0, b = 0;
    auto key() const { return std::tie(period, a, b); }
    bool operator<(const Speed &v) const { return key() < v.key(); }
    bool operator==(const Speed &v) const { return key() == v.key(); }
};

std::ostream &operator<<(std::ostream &o, const Speed &v) {
    return o << ' ' << v.period << ' ' << v.a << ' ' << v.b;
}
Speed speed(Z length, Z a, Z b = 0) {
    Z g = std::gcd(length, std::gcd(a,b));
    return {length/g, a/g, b/g};
}

struct Lane {
    Speed sheet{}, ray{};
    Form arrival{}, clock{};
    auto key() const { return std::tie(sheet, ray, arrival, clock); }
    bool operator<(const Lane &v) const { return key() < v.key(); }
};

struct Summary {
    std::map<Z,Z> by_index;
    Form low, high;
    void add(const Form &clock, Z index) {
        if (by_index.empty() || Wide(clock.c)*low.den < Wide(low.c)*clock.den) low = clock;
        if (by_index.empty() || Wide(clock.c)*high.den > Wide(high.c)*clock.den) high = clock;
        ++by_index[index];
    }
};

const Form A{{1,55,-30,-55,95}, 0, 6};
const Form H = A/2;
const Form prefix{{5,-55,30,55,-95}, 0, 36};
const State sheet{A,H,0,A,H};

struct FiberData {
    std::set<V> fibers;
    std::set<Word> words;
    std::map<Z,Z> by_index;
    std::map<Z,Z> phase_ties, primitive_ties;
};

struct RouteData {
    int beta;
    std::map<Lane, Summary> lanes;
    std::map<Speed,FiberData> fibers;
    std::array<Z, block_limit> entries{};
    std::map<V,Z> tails;
    std::vector<std::pair<Word,Z>> path;
    bool reported = false;
    Z leaves = 0;
    Z min_tail = std::numeric_limits<Z>::max(), max_tail = 0;

    explicit RouteData(int beta) : beta(beta) {}

    Z tail(V x) {
        auto old = tails.find(x);
        if (old != tails.end()) return old->second;
        Z residue = std::max({x[0], x[1]-x[0], x[2]-x[1], x[3]-x[2], x[4]-x[3], -x[4]});
        require(residue > cycle_envelope, "preterminal point inside cyclic envelope");
        std::map<V,Z> seen;
        auto start = x;
        for (int t = 0; t < 32768; ++t) {
            auto [it, added] = seen.emplace(x,t);
            if (!added) {
                require(t-it->second > 0, "zero period");
                min_tail = std::min(min_tail, it->second);
                max_tail = std::max(max_tail, it->second);
                tails.emplace(start, it->second);
                return it->second;
            }
            x = step(x, integer_guard(x));
        }
        throw std::runtime_error("preterminal orbit did not close");
    }

    void classify(const State &start, const Block &block, int depth, const Form &elapsed,
                  Z index, Lane &lane) {
        State current = start;
        V d = block.drift;
        std::set<Speed> found;
        bool center_tie = false, phase_zero_tie = false;
        for (size_t j = 0; j < block.word.size(); ++j) {
            State head{};
            for (int i = 0; i < 5; ++i) head[i] = current[i].leading();
            if (depth == 4 && head[2].small() && d[2] == 0
                && current[2].integer({}) == -3) center_tie = true;
            if (depth == 4 && head == sheet) {
                require(elapsed.leading() == prefix, "sheet prefix clock changed");
                phase_zero_tie |= current[2].integer({}) == -3;
                require(d[2] == 0 && d[0] == d[3] && d[1] == d[4], "sheet drift");
                const V offsets{(current[0]-2*current[1]).integer({}), current[2].integer({}),
                                (current[3]-current[0]).integer({}),
                                (current[4]-current[1]).integer({}), 0};
                for (int i = 1; i <= 3; ++i)
                    require(mod(offsets[i],6) == (3*beta)%6, "sheet fiber lattice");
                Speed v = speed(block.word.size(), -d[0], -d[1]);
                found.insert(v);
                auto &data = fibers[v];
                data.fibers.insert(offsets);
                Word rotation = block.word;
                std::rotate(rotation.begin(), rotation.begin()+j, rotation.end());
                data.words.insert(rotation);
            }
            if (depth >= 5 && head[0] == head[3] && !head[0].small()
                && head[1].small() && head[2].small() && head[4].small()) {
                require(d[0] == d[3] && d[1] == 0 && d[2] == 0 && d[4] == 0, "ray drift");
                found.insert(speed(block.word.size(), -d[0]));
            }
            int code = block.word[j];
            current = step(current, code);
            d = mix(1, step(d,code), -1, step(V{},code));
        }
        if (depth >= 4) require(found.size() == 1, "unclassified macroscopic return");
        if (depth == 4) {
            lane.sheet = *found.begin();
            require(beta%2 || !center_tie, "even center tie");
            auto &data = fibers[lane.sheet];
            ++data.by_index[index];
            if (center_tie) ++data.primitive_ties[index];
            if (phase_zero_tie) ++data.phase_ties[index];
        }
        if (depth >= 5) {
            lane.ray = *found.begin();
            lane.arrival = elapsed.leading();
        }
    }

    void explore(State state, Form elapsed, const V &point, const Basis &basis,
                 Z index, int depth, Lane lane) {
        try {
            require(depth < block_limit, "excess macroscopic return depth");
            ++entries[depth];
            const Block *block = nullptr;
            for (int j = 0; j < short_limit; ++j) {
                block = &macro(state);
                if (block->divisor) break;
                state = step(state, block->word[0]);
                elapsed = elapsed + 1;
            }
            require(block && block->divisor, "missing substantial repeat");
            classify(state, *block, depth, elapsed, index, lane);
            Z m = block->divisor, length = block->word.size();
            auto remainders = split(block->numerator, m, point, basis);
            Z next_index = narrow(Wide(index) * remainders.index);
            auto leading_count = block->numerator.leading() / m;
            State leading_out{};
            for (int i = 0; i < 5; ++i)
                leading_out[i] = state[i].leading() + block->drift[i]*leading_count;
            Form leading_time = elapsed.leading() + length*leading_count;
            bool final = constant(leading_out);
            if (final) {
                require(lane.sheet.period, "missing sheet lane before constant state");
                // A boundary lane with a=2b and mean a-decrease twice the
                // b-decrease can exhaust both coordinates in the sheet itself.
                if (!lane.ray.period)
                    require(depth == 4 && lane.sheet.a == 2*lane.sheet.b,
                            "unexplained direct sheet exhaustion");
                lane.clock = leading_time;
            }
            for (Z r = remainders.residue % remainders.g; r < m; r += remainders.g) {
                Form count = (block->numerator - r)/m + 1;
                require(nonnegative(count - (final ? early_returns : 1)), "too few valid copies");
                Form t = elapsed + length*count;
                State end{};
                for (int i = 0; i < 5; ++i) end[i] = state[i] + block->drift[i]*count;
                if (final) {
                    V before = mix(1, plain(end), -early_returns, block->drift);
                    Z remainder_time = tail(before);
                    lanes[lane].add(t + (remainder_time + 1 - early_returns*length), next_index);
                    ++leaves;
                } else {
                    Z offset = mod(narrow(Wide((r - remainders.residue)/remainders.g)
                                          * remainders.inverse), remainders.index);
                    V next_point = mix(1, point, offset, remainders.shift);
                    require(mod(block->numerator.integer(next_point),m) == r, "wrong child remainder");
                    path.emplace_back(block->word, r);
                    explore(end, t, next_point, remainders.basis, next_index, depth+1, lane);
                    path.pop_back();
                }
            }
        } catch (const std::exception &error) {
            if (!reported) {
                reported = true;
                std::cerr << "UNCERTIFIED beta " << beta << " depth " << depth
                          << " parent share 1/" << index << " because " << error.what() << '\n';
                std::cerr << "representative" << point << '\n';
                for (const auto &v : basis) std::cerr << "basis" << v << '\n';
                for (const auto &[word, rem] : path) {
                    std::cerr << "path remainder " << rem << " word";
                    for (int code : word) std::cerr << ' ' << code;
                    std::cerr << '\n';
                }
                for (const Form &f : state) std::cerr << "state" << f << '\n';
            }
            throw;
        }
    }

    void emit() const {
        std::cout << "beta " << beta << " entries";
        for (Z n : entries) if (n) std::cout << ' ' << n;
        std::cout << " leaves " << leaves << " tails " << tails.size() << " minmax "
                  << min_tail << ' ' << max_tail << " threshold " << threshold << '\n';
        for (const auto &[lane, data] : lanes) {
            std::cout << "lane" << lane.sheet << lane.ray << " arrival" << lane.arrival
                      << " clock" << lane.clock << " offsets " << data.low.c << '/' << data.low.den
                      << ' ' << data.high.c << '/' << data.high.den << " indices";
            for (const auto &[idx,n] : data.by_index) std::cout << ' ' << idx << ':' << n;
            std::cout << '\n';
        }
        for (const auto &[lane,data] : fibers) {
            std::cout << "fiberweights" << lane;
            for (const auto &[idx,n] : data.by_index) std::cout << ' ' << idx << ':' << n;
            std::cout << '\n';
            std::cout << "phase_ties" << lane;
            for (const auto &[idx,n] : data.phase_ties) std::cout << ' ' << idx << ':' << n;
            std::cout << '\n';
            std::cout << "primitive_ties" << lane;
            for (const auto &[idx,n] : data.primitive_ties) std::cout << ' ' << idx << ':' << n;
            std::cout << '\n';
            for (const auto &f : data.fibers) std::cout << "fiber" << lane << f << '\n';
            for (const auto &word : data.words) {
                std::cout << "word" << lane;
                for (int k : word) std::cout << ' ' << k;
                std::cout << '\n';
            }
        }
        std::cout << std::flush;
    }
};

int main(int argc, char **argv) {
    try {
        State start{};
        Basis basis{};
        for (int i = 0; i < 5; ++i) {
            V row{};
            row[i] = 1;
            start[i] = Form(row,0,1);
            basis[i][i] = 6;
        }
        int first = argc > 1 ? std::stoi(argv[1]) : 0;
        int last = argc > 2 ? std::stoi(argv[2]) : 5;
        require(0 <= first && first <= last && last < 6, "expected residue interval in 0..5");
        for (int beta = first; beta <= last; ++beta) {
            // Every gap on the closure exceeds (1/7 + 66/B) for B >= 21000.
            threshold = 21000;
            V point{};
            for (int i = 0; i < 5; ++i) point[i] = (i+1)*beta;
            RouteData result(beta);
            result.explore(start,0,point,basis,1,0,{});
            result.emit();
            block_cache.clear();
        }
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
