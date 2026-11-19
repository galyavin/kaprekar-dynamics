#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using Int = std::int64_t;
using State = std::array<Int, 5>;

constexpr std::array<Int, 12> factorial{
    1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800
};

State step(State q, Int base) {
    std::array<Int, 11> digit{};
    Int borrow = 0;
    for (int i = 10; i >= 0; --i) {
        Int difference = -borrow;
        if (i < 5) difference += q[i];
        if (i > 5) difference -= q[10 - i];
        borrow = difference < 0;
        digit[i] = difference + base * borrow;
    }
    std::sort(digit.begin(), digit.end());
    for (int i = 0; i < 5; ++i) q[i] = digit[10 - i] - digit[i];
    return q;
}

Int differences_time(State input, Int base) {
    // Brent detects the cycle. The saved states locate its first entrance
    // with only one short replay after finding a repeated state.
    constexpr Int stride = 512;
    std::vector<State> checkpoints{input};
    State anchor = input;
    State current = step(input, base);
    Int time = 1, period = 1, power = 1;
    while (anchor != current) {
        if (power == period) {
            anchor = current;
            power *= 2;
            period = 0;
        }
        current = step(current, base);
        ++time;
        ++period;
        if (time % stride == 0) checkpoints.push_back(current);
    }

    const State cycle_start = current;
    std::vector<State> cycle{cycle_start};
    for (current = step(current, base); current != cycle_start; current = step(current, base)) {
        cycle.push_back(current);
    }
    if (static_cast<Int>(cycle.size()) != period) throw std::logic_error("period mismatch");
    std::sort(cycle.begin(), cycle.end());
    auto on_cycle = [&cycle](const State &q) {
        return std::binary_search(cycle.begin(), cycle.end(), q);
    };
    auto first = std::partition_point(checkpoints.begin(), checkpoints.end(),
                                      [&on_cycle](const State &q) { return !on_cycle(q); });
    if (first == checkpoints.begin()) return 0;
    time = (first - checkpoints.begin() - 1) * stride;
    current = *(first - 1);
    while (!on_cycle(current)) {
        current = step(current, base);
        ++time;
    }
    return time;
}

std::array<Int, 6> gaps(const State &q, Int base) {
    return {base - q[0], q[0] - q[1], q[1] - q[2], q[2] - q[3], q[3] - q[4], q[4]};
}

std::uint64_t below(std::uint64_t bound, std::mt19937_64 &random) {
    const auto threshold = -bound % bound;
    for (;;) {
        const auto candidate = random();
        if (candidate >= threshold) return candidate % bound;
    }
}

bool accept_word_weight(const std::array<Int, 6> &g, const std::array<Int, 6> &maxima,
                        std::mt19937_64 &random) {
    // Propose a least digit and five symmetric spacings from a common envelope.
    // Each valid sorted word has the same proposal probability given its q.
    // Its ordered-word multiplicity is 11! divided by the factorials of its
    // tied runs. Accept with that inverse factorial product, since 11! is shared.
    if (below(maxima[0], random) >= static_cast<std::uint64_t>(g[0])) return false;
    std::array<bool, 10> tied{};
    for (int i = 0; i < 5; ++i) {
        const auto split = below(maxima[i + 1] + 1, random);
        const auto gap = static_cast<std::uint64_t>(g[i + 1]);
        if (split > gap) return false;
        tied[i] = split == 0;
        tied[9 - i] = split == gap;
    }
    Int denominator = 1, run = 1;
    for (bool tie : tied) {
        if (tie) {
            ++run;
        } else {
            denominator *= factorial[run];
            run = 1;
        }
    }
    denominator *= factorial[run];
    return denominator == 1 || below(denominator, random) == 0;
}

constexpr State lower{94040, 73770, 57940, 20, 200};
constexpr State upper{94130, 73850, 58015, 80, 300};
constexpr State shift{0, -12, 18, -6, -18};

State chart(const State &q, Int base) {
    return {q[0], q[1], q[2], 2 * (q[1] - q[4]) - base, base - q[2] - q[3]};
}

State invert(const State &x, Int base) {
    return {x[0], x[1], x[2], base - x[2] - x[4], x[1] - (base + x[3]) / 2};
}

bool in_box(State x, Int base) {
    for (int i = 0; i < 5; ++i) {
        if (100000 * x[i] < base * lower[i] || 100000 * x[i] > base * upper[i]) return false;
    }
    return true;
}

class BoxSampler {
    Int base;
    State low{}, high{};
    std::array<Int, 6> maxima;
    std::mt19937_64 random;

public:
    BoxSampler(Int b, std::uint64_t seed) : base(b), random(seed) {
        for (int i = 0; i < 5; ++i) {
            low[i] = (base * lower[i] + 99999) / 100000;
            high[i] = base * upper[i] / 100000;
        }
        low[3] += (base - low[3]) & 1;
        high[3] -= (high[3] - base) & 1;
        maxima = gaps(invert(low, base), base);
        for (int mask = 1; mask < 32; ++mask) {
            State corner{};
            for (int i = 0; i < 5; ++i) corner[i] = mask & (1 << i) ? high[i] : low[i];
            auto g = gaps(invert(corner, base), base);
            for (int i = 0; i < 6; ++i) maxima[i] = std::max(maxima[i], g[i]);
        }
    }

    State draw() {
        for (;;) {
            State x{};
            for (int i = 0; i < 5; ++i) {
                const Int stride = i == 3 ? 2 : 1;
                x[i] = low[i] + stride * below((high[i] - low[i]) / stride + 1, random);
            }
            auto q = invert(x, base);
            if (accept_word_weight(gaps(q, base), maxima, random)) return q;
        }
    }
};

std::array<double, 7> graphs(State chart_numerator, Int base) {
    double x1 = chart_numerator[0], x2 = chart_numerator[1], x3 = chart_numerator[2];
    double x4 = chart_numerator[3], x5 = chart_numerator[4];
    double n = (x1 + x3 - 1.5 * base - 5.5 * x4 - 5.0 * x5) / 22;
    double p = x3 - x4 / 2 - 2 * n;
    double r = x2 - 1.5 * x4 - x5 - 4 * n;
    return {
        x1 - 11.0 * base / 324 - 133.0 * p / 162 - 35.0 * r / 81,
        x1 - 7.0 * base / 408 - 57.0 * p / 68 - 15.0 * r / 34,
        x1 + base / 576.0 - 247.0 * p / 288 - 65.0 * r / 144,
        x1 + 31.0 * base / 96 - 19.0 * p / 16 - 5.0 * r / 8,
        x1 + 29.0 * base / 54 - 38.0 * p / 27 - 20.0 * r / 27,
        x1 + 139.0 * base / 144 - 133.0 * p / 72 - 35.0 * r / 36,
        x1 + 863.0 * base / 918 - 836.0 * p / 459 - 440.0 * r / 459
    };
}

void box_run(Int base, Int draws, Int shift_draws, std::uint64_t seed) {
    BoxSampler sample(base, seed);
    std::array<Int, 7> count{};
    std::array<Int, 3> shifts{}, boundary{}, failures{};
    double worst = 0, margin = std::numeric_limits<double>::infinity();
    double separation = std::numeric_limits<double>::infinity();
    for (Int i = 0; i < std::max(draws, shift_draws); ++i) {
        State q = sample.draw();
        const Int mu = differences_time(q, base);
        if (mu < 1) throw std::logic_error("initial differences already on a cycle");
        const Int word_time = mu + 1;
        if (i < draws) {
            const auto graph = graphs(chart(q, base), base);
            std::array<double, 7> delta{};
            for (int j = 0; j < 7; ++j) delta[j] = std::abs(word_time - graph[j]);
            int best = std::min_element(delta.begin(), delta.end()) - delta.begin();
            count[best]++;
            worst = std::max(worst, delta[best]);
            for (int j = 0; j < 7; ++j) {
                if (j == best) continue;
                margin = std::min(margin, delta[j] - delta[best]);
                separation = std::min(separation, std::abs(graph[j] - graph[best]));
            }
        }
        if (i < shift_draws) {
            int cohort = (q[1] - q[4]) % 3;
            State moved = q;
            for (int j = 0; j < 5; ++j) moved[j] += shift[j];
            if (!in_box(chart(moved, base + 12), base + 12)) {
                boundary[cohort]++;
                continue;
            }
            const Int next_mu = differences_time(moved, base + 12);
            shifts[cohort]++;
            if (next_mu != mu - 10) failures[cohort]++;
        }
    }

    std::cout << "base " << base << " draws " << draws << " seed " << seed << "\ngraph_counts";
    for (auto n : count) std::cout << ' ' << n;
    std::cout << "\nmax_nearest_error " << worst << " min_runner_up_margin " << margin
              << " min_nearest_graph_separation " << separation << '\n';
    for (int cohort = 0; cohort < 3; ++cohort) {
        std::cout << "cohort " << cohort << " overlap " << shifts[cohort]
                  << " discarded " << boundary[cohort] << " failures " << failures[cohort] << '\n';
    }
}

void u_run(Int base, Int draws, std::uint64_t seed) {
    static constexpr std::array<Int, 5> center{26293, 27557, 28481, 40447, 29865};
    static constexpr Int denominator = 9000000;
    std::array<Int, 6> low{}, high{};
    for (int i = 0; i < 5; ++i) {
        low[i] = ((50 * center[i] - 9) * base) / denominator + 1;
        high[i] = ((50 * center[i] + 9) * base - 1) / denominator;
    }
    high[5] = base;
    for (int i = 0; i < 5; ++i) high[5] -= low[i];
    std::mt19937_64 random(seed);
    Int below_ramp = 0, on_ramp = 0, above_ramp = 0;
    for (Int i = 0; i < draws;) {
        std::array<Int, 6> g{};
        g[5] = base;
        for (int j = 0; j < 5; ++j) {
            g[j] = low[j] + below(high[j] - low[j] + 1, random);
            g[5] -= g[j];
        }
        if (!accept_word_weight(g, high, random)) continue;
        State q{};
        q[0] = base - g[0];
        for (int j = 1; j < 5; ++j) q[j] = q[j - 1] - g[j];
        const Int mu = differences_time(q, base);
        if (mu < 1) throw std::logic_error("initial differences already on a cycle");
        const Int t = mu + 1;
        if (500 * t <= 13 * base) ++below_ramp;
        else if (5000 * t < 131 * base) ++on_ramp;
        else ++above_ramp;
        ++i;
    }
    std::cout << "base " << base << " draws " << draws << " seed " << seed
              << " below_ramp " << below_ramp << " on_ramp " << on_ramp
              << " above_ramp " << above_ramp << '\n';
}

int main(int argc, char **argv) {
    if (argc == 6 && std::string(argv[1]) == "box") {
        box_run(std::stoll(argv[2]), std::stoll(argv[3]), std::stoll(argv[4]), std::stoull(argv[5]));
    } else if (argc == 5 && std::string(argv[1]) == "u") {
        u_run(std::stoll(argv[2]), std::stoll(argv[3]), std::stoull(argv[4]));
    } else if (argc == 8 && std::string(argv[1]) == "clock") {
        State q{};
        const Int base = std::stoll(argv[2]);
        for (int i = 0; i < 5; ++i) q[i] = std::stoll(argv[3 + i]);
        std::cout << differences_time(q, base) + 1 << '\n';
    } else {
        std::cerr << "Usage: direct box base draws shift_draws seed | direct u base draws seed"
                     " | direct clock base q1 q2 q3 q4 q5\n";
        return EXIT_FAILURE;
    }
}
