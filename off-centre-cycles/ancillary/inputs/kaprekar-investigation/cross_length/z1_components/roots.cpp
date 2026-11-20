// Exact recurrence and guarded projected-word acceleration for orders 3..20.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using Integer = std::int64_t;
using Root = std::array<Integer, 20>;
using Mask = std::uint32_t;

struct Stepper {
    int order;
    Mask high_bit;
    Mask all_bits;

    explicit Stepper(int r) : order(r), high_bit(1U << (r - 1)), all_bits((1U << r) - 1) {}

    Mask rotate(Mask mask) const {
        return ((mask << 1) & all_bits) | (mask >> (order - 1));
    }

    Integer step(Root &x, bool forced = false) const {
        const Integer z = x[order - 1];
        const Integer sign = forced ? 1 : (z > 0) - (z < 0);
        const Integer loss = forced ? 0 : std::max<Integer>(0, -1 - z);
        for (int i = order - 1; i >= 2; --i) x[i] = x[i - 1] - sign + loss;
        x[1] = x[0] + 2 * loss;
        x[0] = z - 3 * sign + 2 * loss;
        return loss;
    }

    bool equal(const Root &a, const Root &b, Mask omitted) const {
        for (int i = 0; i < order; ++i) {
            if (!(omitted & (1U << i)) && a[i] != b[i]) return false;
        }
        return true;
    }
};

struct Mark {
    int zero_loss;  // -1 if censored
    Integer work;
    Integer clock;
    Integer loss;
    Integer longest_projected_search;
};

struct Arrival {
    Root section{};
    Integer period = 0;
    Integer descent = 0;
    Integer residue = 0;
};

Mark classify(const Stepper &map, Root x, Integer cap, bool trace = false,
              bool accelerate = true, Arrival *arrival = nullptr) {
    Mask active = accelerate ? map.all_bits : 0;
    Integer work = 0;
    Integer clock = 0;
    Integer loss = 0;
    Integer longest_search = 0;
    constexpr Integer coordinate_limit = 1'000'000'000'000'000LL;

    while (work < cap) {
        Root reference = x;
        Integer reference_loss = loss;
        Integer period = 0;
        Integer span = map.order;
        Integer stage_start = clock;
        for (;;) {
            if (active & map.high_bit && x[map.order - 1] <= 0) {
                if (trace) {
                    std::cout << "retire " << (map.order - 1 - clock % map.order)
                              << " at " << clock << " endpoint " << x[map.order - 1] << '\n';
                }
                active ^= map.high_bit;
                longest_search = std::max(longest_search, clock - stage_start);
                break;  // Test the newly inactive endpoint using the ordinary map.
            }
            if (work == cap) return {-1, work, clock, loss, std::max(longest_search, clock - stage_start)};
            loss += map.step(x);
            ++work;
            ++clock;
            ++period;
            active = map.rotate(active);
            if (x[0] <= -coordinate_limit || x[0] >= coordinate_limit ||
                loss >= coordinate_limit || clock >= coordinate_limit) {
                throw std::runtime_error("integer budget exceeded");
            }
            if (!active && map.equal(x, reference, 0)) {
                return {loss == reference_loss, work, clock, loss, longest_search};
            }
            if (active && period % map.order == 0 && map.equal(x, reference, active)) {
                longest_search = std::max(longest_search, clock - stage_start);
                // Read minima relative to the initial height in each active slot.
                Root walk = reference;
                Root minimum;
                minimum.fill(coordinate_limit);
                const bool last_active = arrival && !arrival->period && !(active & (active - 1));
                Root canonical{};
                bool canonical_found = false;
                Integer canonical_height = 0;
                for (Integer i = 0; i < period; ++i) {
                    const int slot = map.order - 1 - i % map.order;
                    if (active & (1U << slot)) {
                        minimum[slot] = std::min(minimum[slot], walk[map.order - 1]);
                        assert(walk[map.order - 1] > 0);
                        if (last_active) {
                            Root projected = walk;
                            projected[map.order - 1] = 0;
                            if (!canonical_found || projected < canonical) {
                                canonical = projected;
                                canonical_height = walk[map.order - 1];
                                canonical_found = true;
                            }
                        }
                    }
                    if (work == cap) return {-1, work, clock, loss, longest_search};
                    map.step(walk);
                    ++work;
                }
                assert(map.equal(walk, x, 0));
                Integer repeats = coordinate_limit;
                for (int i = 0; i < map.order; ++i) {
                    if (!(active & (1U << i))) continue;
                    const Integer descent = reference[i] - x[i];
                    if (descent > 0) {
                        repeats = std::min(repeats, (minimum[i] - 1) / descent);
                        if (last_active) {
                            assert(canonical_found);
                            *arrival = {canonical, period, descent, canonical_height % descent};
                        }
                    }
                }
                if (repeats == coordinate_limit) {
                    // No forced label descends, so ordinary recurrence suffices.
                    active = 0;
                    break;
                }
                if (trace) {
                    std::cout << "return " << clock - period << ' ' << clock
                              << " skip " << repeats << " active " << active << '\n';
                }
                if (clock + __int128(repeats) * period >= coordinate_limit ||
                    loss + __int128(repeats) * (loss - reference_loss) >= coordinate_limit) {
                    throw std::runtime_error("skipped clock or loss exceeds integer budget");
                }
                for (int i = 0; i < map.order; ++i) {
                    const __int128 shifted = x[i] + __int128(repeats) * (x[i] - reference[i]);
                    if (shifted <= -coordinate_limit || shifted >= coordinate_limit) {
                        throw std::runtime_error("skipped height exceeds integer budget");
                    }
                    x[i] = static_cast<Integer>(shifted);
                }
                loss += repeats * (loss - reference_loss);
                clock += repeats * period;
                // At most one further repeated period remains before a retirement.
                stage_start = clock;
                reference = x;
                reference_loss = loss;
                period = 0;
                span = 2 * map.order;
            } else if (period == span) {
                reference = x;
                reference_loss = loss;
                period = 0;
                span *= 2;
            }
        }
    }
    return {-1, work, clock, loss, longest_search};
}

Root read_root(int order) {
    Root x{};
    for (int i = 0; i < order; ++i) {
        if (!(std::cin >> x[i])) throw std::runtime_error("missing root coordinate");
    }
    return x;
}

int main(int argc, char **argv) {
    if (argc < 3) return 2;
    const std::string mode = argv[1];
    const int order = std::stoi(argv[2]);
    if (order < 3 || order > 20) return 2;
    const Stepper map(order);
    if (mode == "panel" || mode == "trace" || mode == "control" || mode == "mix") {
        int count;
        Integer cap;
        if (!(std::cin >> count >> cap)) return 2;
        for (int sample = 0; sample < count; ++sample) {
            const Root root = read_root(order);
            Arrival arrival;
            const Mark marked = classify(map, root, cap, mode == "trace", true,
                                         mode == "mix" ? &arrival : nullptr);
            if (mode == "control") {
                const Mark plain = classify(map, root, cap, false, false);
                if (plain.zero_loss == -1 || plain.zero_loss != marked.zero_loss) {
                    throw std::runtime_error("accelerated and plain marks differ");
                }
            }
            std::cout << marked.zero_loss << ' ' << marked.work << ' ' << marked.clock
                      << ' ' << marked.loss << ' ' << marked.longest_projected_search;
            if (mode == "mix") {
                std::cout << ' ' << arrival.period << ' ' << arrival.descent << ' ' << arrival.residue;
                for (int i = 0; i < order; ++i) std::cout << ' ' << arrival.section[i];
            }
            std::cout << '\n';
        }
        return 0;
    }
    if (mode != "component") return 2;
    // Input marker in original slot, forced prefix length, search cap, and root.
    int active;
    Integer prefix, search_cap;
    if (!(std::cin >> active >> prefix >> search_cap)) return 2;
    Root x = read_root(order);
    const int endpoint = order - 1;
    x[active] = 0;
    for (Integer i = 0; i < prefix; ++i) {
        map.step(x, active == endpoint);
        active = (active + 1) % order;
        x[active] = 0;
    }
    const Root section = x;
    const int section_slot = active;
    struct Read { Root ordinary; Integer offset; };
    std::vector<Read> reads;
    Integer period = 0;
    do {
        if (active == endpoint) reads.push_back({x, x[endpoint]});
        map.step(x, active == endpoint);
        active = (active + 1) % order;
        ++period;
        if (period > search_cap) throw std::runtime_error("projected return not found");
    } while (active != section_slot || !map.equal(x, section, 1U << active));
    const Integer descent = -x[active];
    Integer minimum = std::numeric_limits<Integer>::max();
    for (const auto &read : reads) minimum = std::min(minimum, read.offset);
    const Integer first = 1 - minimum;
    Integer height = first;
    Integer success = 0;
    if (descent <= 0) throw std::runtime_error("nonpositive active descent");
    std::cout << "period " << period << " descent " << descent << " min_read " << minimum
              << " section_slot " << section_slot << " section";
    for (int i = 0; i < order; ++i) std::cout << ' ' << section[i];
    std::cout << '\n';
    for (const auto &read : reads) {
        const Integer upper = std::min(first + descent - 1, descent - read.offset);
        while (height <= upper) {
            Root exit = read.ordinary;
            exit[endpoint] = height - descent + read.offset;
            const Mark marked = classify(map, exit, 10'000'000, false, false);
            if (marked.zero_loss == -1) throw std::runtime_error("censored exact component exit");
            success += marked.zero_loss;
            std::cout << ((height % descent + descent) % descent) << '\t' << marked.zero_loss;
            for (int i = 0; i < order; ++i) std::cout << '\t' << exit[i];
            std::cout << '\n';
            ++height;
        }
    }
    if (height != first + descent) throw std::runtime_error("unclassified incoming residue");
    std::cout << "good " << success << " total " << descent << '\n';
}
