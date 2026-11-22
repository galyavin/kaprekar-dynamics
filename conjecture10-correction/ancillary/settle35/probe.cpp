// Task 245 finite diagnostics. All state, loss and radius computations are integer.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

using Int = std::int64_t;
constexpr int r35 = 12;
constexpr std::array<int, r35> root_slots{1,2,4,8,16,3,6,12,11,13,9,17};

struct Random {
    std::uint64_t seed;
    std::uint64_t next() {
        std::uint64_t z = (seed += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    Int offset(Int half_width) {
        const auto span = static_cast<std::uint64_t>(2*half_width+1);
        const auto threshold = -span % span;
        std::uint64_t z;
        do z = next(); while (z < threshold);
        return static_cast<Int>(z % span) - half_width;
    }
};

int sign(Int x) { return (x > 0) - (x < 0); }
Int loss(Int x) { return std::max<Int>(0, -1-x); }

using Root = std::array<Int, r35>;

Int root_spread(const Root& x) {
    const auto [lo, hi] = std::minmax_element(x.begin(), x.end());
    return std::max<Int>(*hi, 0) - std::min<Int>(*lo, 0);
}

Int secondary(const Root& x) { return x[4] + x.back() - x.front(); }

Root root_step(const Root& x) {
    const Int y = x.back(), h1 = loss(y), h2 = loss(secondary(x));
    const int e = sign(y);
    Root next;
    next[0] = y - 3*e + 2*h1;
    next[1] = x[0] + 2*h1 - h2;
    next[2] = x[1] - e + h1 + h2;
    next[5] = x[4] - e + h2;
    for (int j : {3,4,6,7,8,9,10,11}) next[j] = x[j-1] - e + h1;
    return next;
}

struct Cycle { Int entry, period, detected; };

template<typename State, typename Step>
Cycle detect(const State& origin, const Step& step, Int cap) {
    State first = origin;
    State second = step(origin);
    Int power = 1, period = 1, used = 1;
    while (first != second && used < cap) {
        if (power == period) {
            first = second;
            power *= 2;
            period = 0;
        }
        second = step(second);
        ++period;
        ++used;
    }
    if (first != second) return {-1, 0, used};
    first = origin;
    second = origin;
    for (Int j = 0; j < period; ++j) second = step(second);
    Int entry = 0;
    while (first != second) {
        first = step(first);
        second = step(second);
        ++entry;
    }
    return {entry, period, used};
}

void print_root(const Root& origin, Int L) {
    const auto hit = detect(origin, root_step, 512*L);
    auto x = origin;
    Int spread = 0, before = 0, recurrent = 0, peak_loss = 0;
    const Int length = hit.period ? hit.entry + hit.period : hit.detected;
    for (Int t = 0; t <= length; ++t) {
        spread = std::max(spread, root_spread(x));
        if (t == length) break;
        const Int h = loss(secondary(x));
        peak_loss = std::max(peak_loss, h);
        if (hit.entry < 0 || t < hit.entry) before += h;
        else recurrent += h;
        x = root_step(x);
    }
    const bool budget = hit.period && !recurrent && spread <= 32*L && before <= 512*L;
    std::cout << hit.entry << '\t' << hit.period << '\t' << hit.detected << '\t'
              << spread << '\t' << before << '\t' << recurrent << '\t' << peak_loss
              << '\t' << budget;
    for (auto value : origin) std::cout << '\t' << value;
}

template<int n>
struct Physical {
    static constexpr int k = (n-1)/2, m = k+1, p = m/2;
    using Gap = std::array<Int, m>;
    Int base;
    std::array<int, m> permutation{};
    std::array<int, m> orbit{};
    int orbit_count = 0;

    explicit Physical(Int b) : base(b) {
        for (int j = 1; j < m; ++j) permutation[j] = std::min(2*j, n-2*j);
        orbit.fill(-1);
        orbit[0] = 0;
        for (int j = 1; j < m; ++j) {
            if (orbit[j] >= 0) continue;
            int s = j;
            do { orbit[s] = orbit_count; s = permutation[s]; } while (s != j);
            ++orbit_count;
        }
    }

    Gap make_source(const std::vector<Int>& center, Int width, Random& random) const {
        assert(static_cast<int>(center.size()) == m);
        Gap g{};
        g[0] = base;
        for (int j = 1; j < m; ++j) {
            g[j] = center[j] + random.offset(width);
            g[0] -= g[j];
        }
        assert(g[0] >= 1 && *std::min_element(g.begin()+1, g.end()) >= 0);
        return g;
    }

    Gap suffix(const Gap& g) const {
        Gap q{};
        for (int j = k; j >= 1; --j) q[j] = g[j] + (j == k ? 0 : q[j+1]);
        return q;
    }

    Gap literal(const Gap& g) const {
        const Gap q = suffix(g);
        // D-A has coefficients +q_i and -q_i at the opposite powers.
        std::array<Int, n> digits{};
        for (int i = 1; i <= k; ++i) {
            digits[n-i] = q[i];
            digits[i-1] = -q[i];
        }
        for (int j = 0; j < n-1; ++j) {
            if (digits[j] < 0) {
                digits[j] += base;
                --digits[j+1];
            }
            assert(0 <= digits[j] && digits[j] < base);
        }
        assert(0 <= digits.back() && digits.back() < base);
        std::sort(digits.begin(), digits.end(), std::greater<Int>{});
        Gap result{};
        Int previous = base;
        for (int i = 1; i <= k; ++i) {
            const Int difference = digits[i-1] - digits[n-i];
            result[i-1] = previous - difference;
            previous = difference;
        }
        result[k] = previous;
        return result;
    }

    Int pairing_margin(const Gap& g) const {
        const auto q = suffix(g);
        Int margin = std::numeric_limits<Int>::max(), previous = base - 1;
        for (int i = 1; i <= k; ++i) {
            const Int a = q[i] - (i == k);
            const Int b = base - q[k+1-i] - (i < k);
            margin = std::min(margin, previous - std::max(a,b));
            previous = std::min(a,b);
        }
        return margin;
    }

    std::array<Int,p+1> rows(const Gap& g) const {
        std::array<Int,p+1> out{};
        for (int i = 1; i <= p; ++i) out[i] = out[i-1] + g[m-i] - g[i-1];
        return out;
    }

    Int peak_numerator(const Gap& g) const {
        const auto q = suffix(g);
        Int value = 0;
        for (int i = 1; i <= k; ++i) value = std::max(value, std::abs(m*q[i]-(m-i)*base));
        return value;
    }
};

template<int n>
void print_physical(const Physical<n>& map, const typename Physical<n>::Gap& origin, Int cap, int admitted_row) {
    constexpr int k = Physical<n>::k, m = Physical<n>::m, p = Physical<n>::p;
    const auto literal = [&](const auto& g) { return map.literal(g); };
    const auto hit = detect(origin, literal, cap);
    const Int length = hit.period ? hit.entry + hit.period : hit.detected;
    auto g = origin;
    Root x{};
    if constexpr (n == 35) {
        for (int i = 0; i < r35; ++i) x[i] = g[root_slots[i]] - g[0];
    }
    std::array<int,m> labels{};
    std::iota(labels.begin(), labels.end(), 0);
    std::array<Int,m> export_by_label{};
    std::array<Int,p+1> minimum;
    minimum.fill(std::numeric_limits<Int>::max());
    Int first_height = -1, first_pair = -1, first_guard = -1;
    Int first_leave = -1, first_return = -1, first_secondary = -1, first_row_back = -1;
    Int first_full_leave = -1, first_full_return = -1;
    Int cycle_peak = 0, prefix_pair = std::numeric_limits<Int>::max(), last_compared = -1;
    bool cycle_c = true, cycle_wide = true;
    for (Int t = 0; t < length; ++t) {
        const auto rows = map.rows(g);
        const Int pairing = map.pairing_margin(g);
        bool active_ok = true, all_ok = true;
        for (int i = 2; i <= p; ++i) {
            if (rows[i] >= -1) continue;
            all_ok = false;
            if (i != admitted_row) active_ok = false;
        }
        const bool in_wide = g[k] > 0 && pairing >= 0 && active_ok;
        const bool in_c = in_wide && all_ok;
        if (g[k] <= 0 && first_height < 0) first_height = t;
        if (pairing < 0 && first_pair < 0) first_pair = t;
        if (!active_ok && first_guard < 0) first_guard = t;
        if (!in_wide && first_leave < 0) first_leave = t;
        if (in_wide && first_leave >= 0 && first_return < 0) first_return = t;
        if (!in_c && first_full_leave < 0) first_full_leave = t;
        if (in_c && first_full_leave >= 0 && first_full_return < 0) first_full_return = t;
        if (admitted_row > 0 && rows[admitted_row] < -1 && first_secondary < 0) first_secondary = t;
        if (in_c && first_secondary >= 0 && first_row_back < 0) first_row_back = t;
        if (hit.period && t >= hit.entry) {
            cycle_peak = std::max(cycle_peak, map.peak_numerator(g));
            cycle_c = cycle_c && in_c;
            cycle_wide = cycle_wide && in_wide;
        }
        const auto next = literal(g);
        if (first_leave < 0) {
            for (int i = 1; i <= p; ++i) minimum[i] = std::min(minimum[i], rows[i]);
            prefix_pair = std::min(prefix_pair, pairing);
            last_compared = t;
            if constexpr (n == 35) {
                assert(x.back() == rows[1] && secondary(x) == rows[2]);
                auto projected = g;
                auto next_labels = labels;
                for (int j = 0; j < m; ++j) {
                    projected[map.permutation[j]] = g[j];
                    next_labels[map.permutation[j]] = labels[j];
                }
                const Int e = sign(rows[1]), h1 = loss(rows[1]), h2 = loss(rows[2]);
                projected[0] += e - h1; projected[1] += -2*e + h1;
                projected[2] += e + h1 - h2; projected[3] += -h1 + h2;
                projected[4] += h2; projected[5] -= h2;
                assert(projected == next);
                x = root_step(x);
                for (int i = 0; i < r35; ++i) assert(x[i] == next[root_slots[i]] - next[0]);
                labels = next_labels;
                export_by_label[labels[5]] += h2;
                for (int j = 0; j < m; ++j) if (map.orbit[j] != 0) {
                    assert(next[j] == origin[labels[j]] - export_by_label[labels[j]]);
                }
            }
        }
        g = next;
    }
    std::cout << hit.entry << '\t' << hit.period << '\t' << hit.detected << '\t'
              << cycle_peak << '\t' << (hit.period && cycle_c) << '\t' << (hit.period && cycle_wide)
              << '\t' << first_height << '\t' << first_pair << '\t' << first_guard
              << '\t' << first_leave << '\t' << first_return << '\t' << first_secondary
              << '\t' << first_row_back << '\t' << first_full_leave << '\t' << first_full_return
              << '\t' << last_compared << '\t' << (last_compared < 0 ? -1 : prefix_pair);
    for (int i = 1; i <= p; ++i) std::cout << '\t' << (last_compared < 0 ? 0 : minimum[i]);
    if constexpr (n == 35) for (int j : {5,7,10,14,15}) std::cout << '\t' << export_by_label[j];
    for (auto value : origin) std::cout << '\t' << value;
}

int main(int argc, char** argv) {
    assert(argc >= 3);
    const std::string mode = argv[1];
    const Int scale = std::stoll(argv[2]);
    const int first = argc > 3 ? std::stoi(argv[3]) : 0;
    const int last = argc > 4 ? std::stoi(argv[4]) : (mode == "cube" ? 96 : 24);
    if (mode == "cube") {
        Random random{static_cast<std::uint64_t>(350000+scale)};
        for (int i = 0; i < last; ++i) {
            Root x;
            for (auto& value : x) value = random.offset(scale);
            if (i < first) continue;
            std::cout << scale << '\t' << i << '\t';
            print_root(x, scale);
            std::cout << std::endl;
        }
    } else if (mode == "box35") {
        const Int L = scale/10000;
        const Physical<35> map(scale);
        Random random{static_cast<std::uint64_t>(351000+scale)};
        std::vector<Int> center(map.m);
        for (int j = 1; j < map.m; ++j) center[j] = scale/20 + (map.orbit[j] ? scale/50 : 0);
        for (int i = 0; i < last; ++i) {
            const auto g = map.make_source(center, L, random);
            if (i < first) continue;
            Root x;
            for (int t = 0; t < r35; ++t) x[t] = g[root_slots[t]]-g[0];
            std::cout << scale << '\t' << i << '\t';
            print_root(x, L);
            std::cout << '\t';
            print_physical(map, g, 512*L, 2);
            std::cout << std::endl;
        }
    } else {
        return 2;
    }
    return 0;
}
