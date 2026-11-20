// Exact generic-point stage reader for a rational chart partition.
#define LOTTERY_NO_MAIN
#include "lottery.cpp"
#include <deque>
#include <set>

struct Sheet {
    Z base, period;
    V slope;
    Lattice at_base;
};

// A shift by period in B adds slope, modulo the fixed-base lattice.
Sheet advance(const State &step, const Sheet &in) {
    return {in.base, in.period, mix(1, eval(step, in.slope, in.period), -1, eval(step, {}, 0)),
            push(step, in.at_base, in.base)};
}

V normalized_slope(V t, const Lattice &at_base) {
    Lattice shifted{t, at_base.cols};
    shifted.reduce();
    return shifted.p;
}

Z bezout_inverse(Z a, Z modulus) {
    Z m = modulus, u = 1, v = 0;
    while (m) {
        const Z q = a / m;
        const Z r = a % m;
        const Z c = narrow(W(u) - W(q) * v);
        a = m;
        m = r;
        u = v;
        v = c;
    }
    assert(a == 1 || a == -1);
    return mod(a * u, modulus);
}

struct Child {
    Sheet sheet;
    Z remainder, index;
};

std::vector<Child> children(const Stage &st, const Sheet &in) {
    const Split reference = split(in.at_base, st.exit, st.s);
    const Z h = st.s / reference.index;
    const Z delta = narrow(W(dot(st.exit.row, in.slope)) + W(st.exit.b) * in.period);
    const Z refinement = h / std::gcd(h, delta);
    const Z next_period = narrow(W(in.period) * refinement);
    const Z target = narrow(W(refinement) * delta);
    const Z shift = reference.index == 1 ? 0 : mod(
        narrow(-W(target / h) * bezout_inverse(reference.stride / h, reference.index)),
        reference.index);
    const Z full_delta = narrow(W(target) + W(shift) * reference.stride);
    assert(full_delta % st.s == 0);
    const V advance_on_lane = mix(refinement, in.slope, shift, reference.shift);
    const V advance_after_exit = mix(1, advance_on_lane, full_delta / st.s, st.d);

    std::vector<Child> result;
    for (Z k = 0; k < refinement; ++k) {
        const Z base = in.base + k * in.period;
        Lattice start = in.at_base;
        start.p = mix(1, start.p, k, in.slope);
        const Split lanes = split(start, st.exit, st.s);
        for (Z j = 0; j < lanes.index; ++j) {
            const Edge edge = exit_layer(st, lanes, j, base);
            result.push_back({{base, next_period, normalized_slope(advance_after_exit, edge.next),
                               edge.next}, edge.remainder, lanes.index});
        }
    }
    return result;
}

// phase_guards() uses rows expressed in the initial leading chart. Structural
// equalities are permitted only where the row is constant on the joint sheet.
void phase_guards(const Context &ctx, const Sheet &sheet, const Phase &phase,
                  const Chart &chart, std::set<Linear> &guards) {
    for (const Form &g : phase.guards) {
        const Linear row = pull(g, chart);
        if (row.row == V{} && row.b == 0) {
            assert(!sheet.at_base.varies(g));
            assert(dot(g.row, sheet.slope) + g.b * sheet.period == 0);
            assert(g.eval(sheet.at_base.p, sheet.base) >= 0);
        } else {
            if (g.lead(ctx.x) <= 0) throw Ambiguous{};
            guards.insert(row);
        }
    }
}

struct ShortFrame {
    Point point;
    Chart chart;
    Sheet sheet;
    Phase phase;
};

struct ShortLookahead {
    std::deque<ShortFrame> frames;
    int verified = 0;

    bool disproves_every_return(const Point &point, const Chart &chart, const Sheet &sheet,
                                std::set<Linear> &guards) {
        if (sheet.at_base.cols.size() > 2) return false;
        try {
            while (frames.size() < 2 * max_word) {
                Point p = point;
                Chart c = chart;
                Sheet s = sheet;
                if (!frames.empty()) {
                    const auto &f = frames.back();
                    p = lead(f.phase.end, f.point);
                    c = advance(f.phase.end, f.chart);
                    s = advance(f.phase.end, f.sheet);
                }
                if (p.zero()) return false;
                const Phase ph = primitive(identity(), {p, s.at_base, s.base});
                frames.push_back({p, c, s, ph});
            }
            for (int length = 1; length <= max_word; ++length) {
                if (frames[length].point.n != point.n) continue;
                bool repeated = true;
                for (int k = 0; k < length && repeated; ++k) {
                    const Phase &a = frames[k].phase, &b = frames[k + length].phase;
                    repeated = std::tie(a.word, a.paired, a.height) == std::tie(b.word, b.paired, b.height);
                }
                if (repeated) return false;
            }
            // A positive macroscopic return must repeat its labelled word twice
            // before reaching any declining comparison, for all sufficiently large B.
            for (; verified < 2 * max_word; ++verified) {
                const auto &f = frames[verified];
                phase_guards({f.point, f.sheet.at_base, f.sheet.base}, f.sheet, f.phase, f.chart, guards);
            }
            return true;
        } catch (const Ambiguous &) {
            frames.clear();
            verified = 0;
            return false;
        } catch (const ZeroWord &) {
            frames.clear();
            verified = 0;
            return false;
        }
    }

    void discard_first() {
        if (frames.empty()) return;
        frames.pop_front();
        verified = std::max(0, verified - 1);
    }
};

void examine(Point point, Chart chart, Sheet sheet, int short_cap) {
    std::set<Linear> guards;
    int shorts = 0;
    int max_probe = 0;
    bool found = false;
    ShortLookahead short_lookahead;
    try {
        while (!point.zero() && shorts <= short_cap) {
            if (shorts >= 5 && shorts < short_cap &&
                short_lookahead.disproves_every_return(point, chart, sheet, guards)) {
                const State next = short_lookahead.frames.front().phase.end;
                chart = advance(next, chart);
                point = lead(next, point);
                sheet = advance(next, sheet);
                short_lookahead.discard_first();
                ++shorts;
                continue;
            }
            const Context context{point, sheet.at_base, sheet.base};
#ifdef STAGE_GUARD_PROBE
            std::vector<Phase> probe;
            recorded_stage_phases = &probe;
            const Stage st = stage(context);
            recorded_stage_phases = nullptr;
#else
            const Stage st = stage(context);
#endif
            State composed = identity();
            std::map<Linear, std::pair<Z, Form>> exit_options;
            const auto old_energy = energy(point);
#ifdef STAGE_GUARD_PROBE
            const int phases = int(probe.size());
#else
            const int phases = st.s ? int(st.words.size()) : max_word;
#endif
            for (int k = 0; k < phases; ++k) {
#ifdef STAGE_GUARD_PROBE
                const Phase &ph = probe[k];
#else
                const Phase ph = primitive(composed, context);
#endif
                if (st.s) assert(ph.word == st.words[k]);
                phase_guards(context, sheet, ph, chart, guards);
                if (st.s) for (const Form &g : ph.guards) {
                    const Z s = -dot(g.row, st.d);
                    if (s <= 0) continue;
                    const Linear distance = pull(g, chart) / s;
                    if (g.lead(point) <= 0) throw Ambiguous{};
                    auto [option, fresh] = exit_options.try_emplace(distance, s, g);
                    if (fresh) continue;
                    const auto &[other_s, other] = option->second;
                    const Form diff = other_s * g - s * other;
                    assert(!sheet.at_base.varies(diff));
                    assert(dot(diff.row, sheet.slope) + diff.b * sheet.period == 0);
                    if (diff.eval(sheet.at_base.p, sheet.base) < 0) option->second = {s, g};
                }
                composed = ph.end;
                max_probe = std::max(max_probe, k + 1);
                if (!st.s && point.n[4] + point.den > 0 && energy(lead(composed, point)) < old_energy)
                    break;
            }
            if (st.s) {
                const Chart returned = advance(composed, chart);
                for (int i = 0; i < 5; ++i) {
                    const Linear a = returned[i], b = chart[i];
                    if (a < b || b < a) throw Ambiguous{};
                }
                std::cout << "stage " << shorts << ' ' << max_probe << ' ' << st.words.size()
                          << ' ' << st.failures << ' ' << sheet.at_base.cols.size() << '\n';
                std::cout << exit_options.size() << '\n';
                for (const auto &[distance, option] : exit_options) {
                    Stage chosen = st;
                    std::tie(chosen.s, chosen.exit) = option;
                    const Chart out = project(chosen, chart);
                    const bool central = std::all_of(out.begin(), out.end(), [](const Linear &row) {
                        return row.row == V{} && row.b == 0;
                    });
                    if (central) {
                        for (const V &col : sheet.at_base.cols) {
                            assert(mix(chosen.s, col, dot(chosen.exit.row, col), chosen.d) == V{});
                        }
                        const Z shift_guard = narrow(W(dot(chosen.exit.row, sheet.slope)) +
                                                     W(chosen.exit.b) * sheet.period);
                        assert(mix(chosen.s, sheet.slope, shift_guard, chosen.d) == V{});
                    }
                    const auto next = central && sheet.at_base.cols.size() == 1 ?
                                      std::vector<Child>{} : children(chosen, sheet);
                    std::cout << distance << '\n';
                    std::cout << chosen.s << chosen.d << chosen.exit.row
                              << ' ' << chosen.exit.b << ' ' << chosen.exit.c << '\n';
                    for (const Linear &f : out) std::cout << ' ' << f;
                    std::cout << '\n';
                    std::cout << next.size() << '\n';
                    for (const auto &[a, rem, index] : next) {
                        std::cout << a.base << ' ' << a.period << ' ' << rem << ' ' << index
                                  << ' ' << a.at_base.cols.size() << a.slope << a.at_base.p;
                        for (const V &c : a.at_base.cols) std::cout << c;
                        std::cout << '\n';
                    }
                }
                found = true;
                break;
            }
            if (shorts == short_cap) break;
            chart = advance(st.first.end, chart);
            point = lead(st.first.end, point);
            sheet = advance(st.first.end, sheet);
            short_lookahead.discard_first();
            ++shorts;
        }
        if (point.zero()) {
            for (const Linear &f : chart) assert(f.row == V{} && f.b == 0);
            std::cout << "terminal " << shorts << ' ' << max_probe << '\n';
        } else if (!found && shorts == short_cap) {
            std::cout << "cap " << shorts << ' ' << max_probe << '\n';
        }
    } catch (const Ambiguous &) {
#ifdef STAGE_GUARD_PROBE
        recorded_stage_phases = nullptr;
#endif
        std::cout << "ambiguous " << shorts << ' ' << max_probe << '\n';
    } catch (const ZeroWord &) {
#ifdef STAGE_GUARD_PROBE
        recorded_stage_phases = nullptr;
#endif
        std::cout << "zero_word " << shorts << ' ' << max_probe << '\n';
    }
    std::cout << guards.size() << '\n';
    for (const Linear &f : guards) std::cout << f << '\n';
    std::cout << std::flush;
}

int main() {
    Z den;
    while (std::cin >> den) {
        Point point{{}, den};
        for (W &x : point.n) { Z value; std::cin >> value; x = value; }
        Chart chart;
        for (Linear &f : chart) {
            for (Z &x : f.row) std::cin >> x;
            std::cin >> f.b >> f.den;
        }
        Sheet sheet;
        int n, short_cap;
        std::cin >> sheet.base >> sheet.period >> n >> short_cap;
        for (Z &x : sheet.slope) std::cin >> x;
        for (Z &x : sheet.at_base.p) std::cin >> x;
        sheet.at_base.cols.resize(n);
        for (V &c : sheet.at_base.cols) for (Z &x : c) std::cin >> x;
        examine(point, chart, sheet, short_cap);
    }
}
