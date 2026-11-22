#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>

using Integer = std::int64_t;
constexpr unsigned digits = 35;
constexpr unsigned m = 18;
constexpr unsigned k = m - 1;
constexpr unsigned period_bound = 66;
using Word = std::array<Integer, digits>;
using Gaps = std::array<Integer, m>;
using Differences = std::array<Integer, k>;

void check(bool condition, const char* message, Integer time = 0) {
    if (!condition) {
        std::cerr << "FAIL at step " << time << ": " << message << '\n';
        std::exit(1);
    }
}

unsigned fold(unsigned j) {
    return std::min(2 * j, digits - 2 * j);
}

Word literal_subtract(const Word& sorted, Integer base) {
    Word out{};
    Integer borrow = 0;
    for (unsigned j = digits; j-- > 0;) {
        Integer column = sorted[digits - 1 - j] - sorted[j] - borrow;
        borrow = column < 0;
        out[j] = column + borrow * base;
    }
    check(borrow == 0, "signed subtraction ended with a borrow");
    return out;
}

Word original(const Gaps& gaps) {
    std::array<Integer, digits - 1> adjacent{};
    for (unsigned j = 0; j < k; ++j) {
        adjacent[j] = gaps[j + 1] / 2;
        adjacent[digits - 2 - j] = gaps[j + 1] - adjacent[j];
    }
    Word result{};
    result[0] = gaps[0] / 3;
    for (unsigned j = 0; j < adjacent.size(); ++j) {
        result[j + 1] = result[j] + adjacent[j];
    }
    return result;
}

Differences symmetric_differences(const Word& sorted) {
    Differences q{};
    for (unsigned j = 0; j < k; ++j) {
        q[j] = sorted[digits - 1 - j] - sorted[j];
    }
    return q;
}

Gaps gaps_from(const Differences& q, Integer base) {
    Gaps gaps{};
    gaps[0] = base - q[0];
    for (unsigned j = 1; j < k; ++j) {
        gaps[j] = q[j - 1] - q[j];
    }
    gaps[k] = q[k - 1];
    return gaps;
}

struct Geometry {
    Integer least_pair_margin;
    Integer least_outer_row;
    Integer second_row;
    Integer radius_numerator;
};

Geometry measure(const Gaps& gaps, const Differences& q, Integer base) {
    Integer pair_margin = base;
    Integer outer_row = base;
    Integer radius = 0;
    Integer previous_min = base - 1;
    for (unsigned j = 0; j < k; ++j) {
        const Integer a = q[j] - (j == k - 1);
        const Integer b = base - q[k - 1 - j] - (j < k - 1);
        pair_margin = std::min(pair_margin, previous_min - std::max(a, b));
        previous_min = std::min(a, b);
        const Integer centred = m * q[j] - (m - j - 1) * base;
        radius = std::max(radius, std::abs(centred));
    }
    Integer row = gaps[k] - gaps[0];
    row += gaps[k - 1] - gaps[1];
    const Integer d2 = row;
    for (unsigned j = 2; j < m / 2; ++j) {
        row += gaps[k - j] - gaps[j];
        outer_row = std::min(outer_row, row);
    }
    return {pair_margin, outer_row, d2, radius};
}

Gaps projected_update(const Gaps& gaps, Integer h2) {
    Gaps predicted{};
    for (unsigned j = 0; j < m; ++j) {
        predicted[fold(j)] = gaps[j];
    }
    const Integer d1 = gaps[k] - gaps[0];
    const Integer sign = (d1 > 0) - (d1 < 0);
    const Integer h1 = std::max<Integer>(0, -1 - d1);
    predicted[0] += sign - h1;
    predicted[1] += -2 * sign + h1;
    predicted[2] += sign + h1 - h2;
    predicted[3] += -h1 + h2;
    predicted[4] += h2;
    predicted[5] -= h2;
    return predicted;
}

struct Record {
    Word word;
    Gaps gaps;
    Geometry geometry;
};

void run(Integer base, Integer scale, const Gaps& start) {
    Word word = original(start);
    check(word.back() < base && word.front() >= 0, "out-of-range original word");
    check(std::adjacent_find(word.begin(), word.end()) == word.end(), "repeated original digit");
    std::cout << "base " << base << " L " << scale << " original";
    for (auto digit : word) std::cout << ' ' << digit;
    std::cout << '\n' << std::flush;
    Gaps expected = start;
    std::array<Record, period_bound> ring{};
    std::array<unsigned, 5> locations{5, 7, 10, 14, 15};
    const auto labels = locations;
    std::array<Integer, 5> loss{};
    Integer first_gap_entry = -1;
    Integer total_export = 0, breached_phases = 0;
    Integer least_pair_margin = base, least_outer_row = base;
    const Integer stop = 3785731 * scale + period_bound + 1;
    for (Integer t = 0; t <= stop; ++t) {
        Word sorted = word;
        std::sort(sorted.begin(), sorted.end());
        const Differences q = symmetric_differences(sorted);
        const Gaps gaps = gaps_from(q, base);
        check(gaps == expected, "literal gap state differs from reconstructed update", t);
        check(gaps.back() > 0, "lost full height", t);
        const Geometry geo = measure(gaps, q, base);
        check(geo.least_pair_margin >= 0, "borrowed levels overlap", t);
        check(geo.least_outer_row >= -1, "lost a complementary-row guard", t);
        least_pair_margin = std::min(least_pair_margin, geo.least_pair_margin);
        least_outer_row = std::min(least_outer_row, geo.least_outer_row);
        for (unsigned i = 0; i < locations.size(); ++i) {
            check(gaps[locations[i]] == start[labels[i]] - loss[i], "labelled export disagrees", t);
        }
        Record& old = ring[t % period_bound];
        if (t >= period_bound && first_gap_entry < 0 && gaps == old.gaps) {
            first_gap_entry = t - period_bound;
        }
        if (t >= period_bound && word == old.word) {
            const Integer word_entry = t - period_bound;
            unsigned period = 1;
            while (period < period_bound && ring[(word_entry + period) % period_bound].word != word) {
                ++period;
            }
            check(period_bound % period == 0 && period % 11 == 0, "unexpected least period", t);
            check(first_gap_entry == word_entry || first_gap_entry == word_entry - 1,
                  "word and gap entries differ by more than expected", t);
            Integer min_d2 = base, min_radius = m * base, max_radius = 0;
            for (const Record& phase : ring) {
                check(phase.geometry.second_row >= -1, "recurrent phase outside original C", t);
                check(396 * phase.geometry.radius_numerator > m * base, "recurrent phase below cut", t);
                min_d2 = std::min(min_d2, phase.geometry.second_row);
                min_radius = std::min(min_radius, phase.geometry.radius_numerator);
                max_radius = std::max(max_radius, phase.geometry.radius_numerator);
            }
            check(total_export >= 461994 * scale && total_export <= 2940551 * scale,
                  "literal lifetime export violates certified bound", t);
            std::cout << "word entry " << word_entry << " gap entry " << first_gap_entry
                      << " least period " << period << " all-path min pairing margin " << least_pair_margin
                      << " all-path min outer row " << least_outer_row << " terminal min d2 " << min_d2
                      << " total export " << total_export << " transient d2 breach phases " << breached_phases
                      << " terminal radius numerators " << min_radius << ' ' << max_radius
                      << " denominator " << m * base << " label losses";
            for (auto amount : loss) std::cout << ' ' << amount;
            std::cout << '\n' << std::flush;
            return;
        }
        old = {word, gaps, geo};
        const Integer h2 = std::max<Integer>(0, -1 - geo.second_row);
        total_export += h2;
        breached_phases += h2 > 0;
        for (unsigned i = 0; i < locations.size(); ++i) {
            locations[i] = fold(locations[i]);
            if (locations[i] == 5) loss[i] += h2;
        }
        expected = projected_update(gaps, h2);
        word = literal_subtract(sorted, base);
    }
    check(false, "no word recurrence by the certificate's horizon", stop);
}

int main() {
    std::string header;
    std::getline(std::cin, header);
    Integer base, scale;
    unsigned controls = 0;
    while (std::cin >> base >> scale) {
        Gaps gaps{};
        for (auto& gap : gaps) check(bool(std::cin >> gap), "incomplete input");
        check(std::accumulate(gaps.begin(), gaps.end(), Integer{0}) == base, "gap sum mismatch");
        run(base, scale, gaps);
        ++controls;
    }
    check(controls == 2, "expected two controls");
}
