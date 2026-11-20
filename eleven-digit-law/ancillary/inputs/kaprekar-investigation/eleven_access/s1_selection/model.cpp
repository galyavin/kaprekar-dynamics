#define R2_NO_MAIN
#include "diagnostic.cpp"

struct Passage {
    Point after;
    Ratio time;
    V incoming, outgoing;
    std::vector<V> columns, next_columns;
    Form guard;
    Z length, decrement, index, residue;
    bool lost;
};

struct Model {
    Ratio time;
    int status = 0;
    std::vector<Passage> passages;
};

bool equal(Ratio a, Ratio b) { return a.n == b.n && a.d == b.d; }
bool equal(const Point &a, const Point &b) { return a.n == b.n && a.den == b.den; }

double clock_score(Ratio t) { return std::clamp(double(t.n)/double(t.d), 0.0, 1.0); }

Model model(Point x, Z base, Rng &rng) {
    Lattice l = root_lattice(base);
    Model out;
    int shorts = 0;
    try {
        while (!x.zero() && out.passages.size() < max_stage) {
            const Stage st = stage({x,l,base});
            if (!st.s) {
                if (++shorts > max_short) { out.status = 1; return out; }
                x = lead(st.first.end, x);
                l = push(st.first.end, l, base);
                continue;
            }
            const Z length = Z(st.words.size());
            const Ratio t{W(length)*st.exit.lead(x), x.den*st.s};
            const Split lanes = split(l, st.exit, st.s);
            const Edge edge = exit_layer(st, lanes, Z(rng.below(lanes.index)), base);
            x = projected(x, st);
            out.time = out.time+t;
            out.passages.push_back({x,t,l.p,edge.next.p,l.cols,edge.next.cols,st.exit,
                                    length,st.s,lanes.index,edge.remainder,st.failures!=0});
            l = edge.next;
        }
    } catch (const Ambiguous &) { out.status = 2; }
      catch (const ZeroWord &) { out.status = 4; }
    if (!out.status && !x.zero()) out.status = 3;
    return out;
}

Ranks indexed_ranks(U seed, U row) {
    Rng source{indexed_seed(seed,row)};
    Ranks ranks{};
    for (U &rank : ranks) rank = source.next();
    return ranks;
}

Point model_input(const Ranks &ranks, int &leading_entry) {
    Digits word{};
    std::transform(ranks.begin(),ranks.end(),word.begin(),[](U rank) { return Z(rank>>26); });
    V q = differences(word);
    leading_entry = entry(q,grid,500);
    return centered(q,grid);
}

#ifndef S1_NO_MAIN
int main(int argc, char **argv) {
    if (argc!=4) return 2;
    const std::string mode = argv[1];
    const U begin = std::stoull(argv[2]), count = std::stoull(argv[3]);
    std::cout.precision(17);
    for (U i=begin; i<begin+count; ++i) {
        int entered;
        const bool premise = mode=="premise";
        const U seeds = premise ? source_seed : 14078212001ULL;
        const U stream = premise ? lane_seed : 14078212007ULL;
        const Point p = model_input(indexed_ranks(seeds,i),entered);
        std::cout << i << ' ' << entered;
        if (premise) {
            for (int scale : {12000,48000,192000}) for (int r : {3,5,15,17}) {
                Rng rng{indexed_seed(stream,i)};
                const Trail a = entered<0 ? Trail{} : trace(p,root_lattice(scale+r),scale+r,rng,false);
                std::cout << ' ' << (entered<0 ? 5 : a.status) << ' ' << clock_score(a.count.elapsed);
            }
        } else if (entered>=0) {
            for (int r=0; r<24; ++r) {
                const int base = 192000+r;
                Rng rng{indexed_seed(stream,i)};
                const Model a = model(p,base,rng);
                std::cout << ' ' << a.status << ' ' << clock_score(a.time)
                          << ' ' << a.passages.size();
                for (const Passage &s : a.passages) {
                    std::cout << ' ' << double(s.time.n)/double(s.time.d)
                              << ' ' << s.columns.size() << ' ' << s.next_columns.size()
                              << ' ' << s.length << ' ' << s.guard.b << ' ' << s.decrement
                              << ' ' << s.index << ' ' << s.residue << ' ' << s.lost
                              << ' ' << s.incoming << ' ' << s.outgoing;
                    for (const auto &n : s.after.n) std::cout << ' ' << numstr(n);
                    std::cout << ' ' << numstr(s.after.den);
                }
            }
        }
        std::cout << '\n';
    }
}
#endif
