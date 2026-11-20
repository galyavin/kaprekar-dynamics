#define S1_NO_MAIN
#include "model.cpp"

V centered_integer(const Q &q, Z base) {
    V x{};
    for (int i=0; i<5; ++i) x[i] = 6*q[i]-(5-i)*base;
    return x;
}

Q uncenter(const V &x, Z base) {
    Q q{};
    for (int i=0; i<5; ++i) {
        const Z n = x[i]+(5-i)*base;
        assert(n%6==0);
        q[i] = int(n/6);
    }
    return q;
}

struct Repeat {
    V after{}, before{};
    int steps = 0, loss = 0;
};

Repeat physical_repeat(V actual, int base, int work_limit) {
    // A zero leading form forces the exact affine comparison on this singleton.
    const Point zero_lead{{},0};
    const Lattice singleton{actual,{}};
    const Context ctx{zero_lead,singleton,base};
    State state = identity();
    std::vector<Form> guards;
    int loss = 0;
    for (int length=1; length<=40; ++length) {
        Form outside = -1*state[4];
        for (int i=0; i<5; ++i) {
            const Form r = i ? state[i]-state[i-1] : state[0];
            if (r.eval(actual,base)>outside.eval(actual,base)) outside = r;
        }
        outside.c -= 505;
        if (outside.eval(actual,base)<0) break;
        guards.push_back(outside);
        const State before = state;
        const Phase phase = primitive(state,ctx);
        state = phase.end;
        guards.insert(guards.end(),phase.guards.begin(),phase.guards.end());
        if (!phase.paired) loss |= phase.height==5 ? 1 : 2;
        const V drift = mix(1,eval(state,actual,base),-1,actual);
        if (drift==V{}) break;
        if (std::any_of(drift.begin(),drift.end(),[](Z n) { return std::abs(n)>720; })) continue;
        if (linear(state,drift)!=drift) continue;
        Z returns = work_limit/length;
        for (const Form &guard : guards) {
            const Z decrement = -dot(guard.row,drift);
            const Z slack = guard.eval(actual,base);
            assert(slack>=0);
            if (decrement>0) returns = std::min(returns,1+slack/decrement);
        }
        if (returns<8) continue;
        return {mix(1,actual,returns,drift),eval(before,mix(1,actual,returns-1,drift),base),
                int(length*returns),loss};
    }
    return {};
}

struct Accelerated {
    Physical result;
    int skipped = 0;
    bool fallback = false;
};

Accelerated physical_clock(const Word &word, int base, std::vector<Q> &path) {
    auto fallback = [&]() { return Accelerated{observe(word,base,path),0,true}; };
    if (base<6222) return fallback();
    Q q = word_differences(word), previous{};
    Physical result;
    int t = 0, skipped = 0, searches = 0;
    while (t<4*base && searches<150) {
        const V actual = centered_integer(q,base);
        W r = radius(Point{{actual[0],actual[1],actual[2],actual[3],actual[4]},1}).n;
        if (r<=504) {
            if (result.strict<0 || result.five<0) return fallback();
            const Physical suffix = observe(t ? schoolbook(previous,base) : word,base,path);
            if (suffix.clock<0) return fallback();
            result.clock = t+suffix.clock;
            return {result,skipped,false};
        }
        const auto phase = physical_phase(q,base);
        result.read(q,phase,t,base);
        if (result.strict>=0 && result.five>=0 && q[4]>0) {
            ++searches;
            Repeat rep;
            try { rep = physical_repeat(actual,base,4*base-t); }
            catch (const ZeroWord &) {}
            if (rep.steps) {
                result.strict_loss |= rep.loss;
                result.five_loss |= rep.loss;
                q = uncenter(rep.after,base);
                previous = uncenter(rep.before,base);
                t += rep.steps;
                skipped += rep.steps;
                continue;
            }
        }
        previous = q;
        q = phase.next;
        ++t;
    }
    return fallback();
}

Word physical_word(const Ranks &ranks, U row, int base, U seed) {
    Rng retry{indexed_seed(seed ^ (U(base)<<32),row)};
    const U threshold = -U(base)%U(base);
    Word result{};
    auto digit = [&](U rank) {
        __uint128_t p = __uint128_t(rank)*base;
        while (U(p)<threshold) p = __uint128_t(retry.next())*base;
        return int(p>>64);
    };
    for (int i=0; i<11; ++i) result[i] = digit(ranks[i]);
    while (std::all_of(result.begin(),result.end(),[&](int d) { return d==result[0]; }))
        for (int &d : result) d = digit(retry.next());
    return result;
}

void verify_clock() {
    constexpr U validation_seed = 14052172411ULL;
    Rng source{validation_seed};
    std::vector<Q> fast_path, slow_path;
    int checked = 0, accelerated = 0;
    Z skipped = 0;
    for (const int base : {3534,12007,192003,768005,3072023}) {
        for (int i=0; i<24; ++i) {
            Word word{};
            for (int &digit : word) digit = int(source.below(base));
            const Accelerated fast = physical_clock(word,base,fast_path);
            const Physical slow = observe(word,base,slow_path);
            assert(fast.result.clock==slow.clock && fast.result.strict==slow.strict &&
                   fast.result.five==slow.five && fast.result.strict_loss==slow.strict_loss &&
                   fast.result.five_loss==slow.five_loss);
            checked++;
            accelerated += fast.skipped>0;
            skipped += fast.skipped;
        }
    }
    std::cout << checked << " independent physical clock and entry comparisons, " << accelerated
              << " accelerated, " << skipped << " exact updates skipped\n";
}

#ifndef S1_PHYSICAL_NO_MAIN
int main(int argc, char **argv) {
    if (argc==2 && std::string(argv[1])=="verify") { verify_clock(); return 0; }
    if (argc!=4) return 2;
    const bool premise = std::string(argv[1])=="premise";
    const U seed = premise ? source_seed : 14089004531ULL;
    const U lane = premise ? lane_seed : 14089004547ULL;
    const int anchor = premise ? 192000 : std::stoi(argv[1]);
    const U start = std::stoull(argv[2]), n = std::stoull(argv[3]);
    const std::vector<int> residues = premise ? std::vector<int>{3,5,15,17}
                                              : std::vector<int>{3,5,9,11,15,17,0,1};
    std::vector<Q> path;
    std::cout.precision(17);
    for (U i=start; i<start+n; ++i) {
        const Ranks ranks = indexed_ranks(seed,i);
        int entered;
        const Point p = model_input(ranks,entered);
        std::cout << i << ' ' << entered;
        for (const int r : residues) {
            const int base = anchor+r;
            const Word word = physical_word(ranks,i,base,seed);
            const Accelerated a = physical_clock(word,base,path);
            const auto &v = a.result;
            Rng rng{indexed_seed(lane,i)};
            const bool model_active = premise || (i<12000 && r!=9 && r!=11);
            const Trail g = entered<0 || !model_active ? Trail{} : trace(p,root_lattice(base),base,rng,false);
            std::cout << ' ' << v.clock << ' ' << v.strict << ' ' << v.strict_loss
                      << ' ' << v.five << ' ' << v.five_loss << ' ' << a.skipped << ' ' << a.fallback
                      << ' ' << (entered<0 ? 5 : !model_active ? 6 : g.status) << ' ' << (g.count.first_loss>=0)
                      << ' ' << clock_score(g.count.elapsed);
        }
        std::cout << '\n';
    }
}
#endif
