// Enumerate exact lanes on one open cell of the frozen stage graph.
#define main simulation_main
#include "graph.cpp"
#undef main

#include <set>

struct Linear {
    V row{};
    Z b = 0, den = 1;
    void reduce() {
        Z g = std::gcd(den, b);
        for (Z c : row) g = std::gcd(g,c);
        den /= g;
        b /= g;
        for (Z &c : row) c /= g;
    }
    bool operator<(const Linear &rhs) const {
        return std::tie(row,b,den) < std::tie(rhs.row,rhs.b,rhs.den);
    }
};
Linear operator+(const Linear &a, const Linear &b) {
    const Z g = std::gcd(a.den,b.den), s = b.den/g, t = a.den/g;
    Linear f{mix(s,a.row,t,b.row),narrow(W(s)*a.b+W(t)*b.b),narrow(W(s)*a.den)};
    f.reduce();
    return f;
}
Linear operator*(Z k, const Linear &x) {
    Linear f{mix(k,x.row,0,{}),narrow(W(k)*x.b),x.den};
    f.reduce();
    return f;
}
Linear operator/(const Linear &x, Z k) {
    Linear f{x.row,x.b,narrow(W(k)*x.den)};
    f.reduce();
    return f;
}
std::ostream &operator<<(std::ostream &out, const Linear &c) {
    return out << c.row << ' ' << c.b << ' ' << c.den;
}
Ratio operator*(const Ratio &a, const Ratio &b) {
    const W g = gcdw(a.n,b.d), h = gcdw(b.n,a.d);
    return {a.n/g*(b.n/h),a.d/h*(b.d/g)};
}

using Chart = std::array<Linear,5>;
Chart root() {
    Chart c{};
    for (int i = 0; i < 5; ++i) c[i].row[i] = 1;
    return c;
}
Linear pull(const Form &f, const Chart &c) {
    Linear result{{},f.b,1};
    for (int i = 0; i < 5; ++i) if (f.row[i]) result = result + f.row[i]*c[i];
    return result;
}
Chart advance(const State &step, const Chart &c) {
    Chart out{};
    for (int i = 0; i < 5; ++i) out[i] = pull(step[i],c);
    return out;
}
Chart project(const Stage &st, const Chart &c) {
    const Linear unit = pull(st.exit,c)/st.s;
    Chart out{};
    for (int i = 0; i < 5; ++i) out[i] = c[i]+st.d[i]*unit;
    return out;
}

struct Key {
    V p;
    std::vector<V> cols;
    Chart chart;
    int stage, short_budget;
    bool operator<(const Key &rhs) const {
        return std::tie(p,cols,chart,stage,short_budget) <
               std::tie(rhs.p,rhs.cols,rhs.chart,rhs.stage,rhs.short_budget);
    }
};
struct EdgeType {
    int id, n, rank_in, rank_out, failures;
    Z s, index;
    V d;
    Linear distance;
    bool operator<(const EdgeType &x) const {
        return std::tie(id,n,rank_in,rank_out,failures,s,index,d,distance) <
               std::tie(x.id,x.n,x.rank_in,x.rank_out,x.failures,x.s,x.index,x.d,x.distance);
    }
};
struct Answer {
    std::map<Linear,Ratio> clocks;
    std::map<int,Ratio> missed;
    std::map<EdgeType,Ratio> edges;
};
using Next = std::pair<V,std::vector<V>>;

struct Ledger {
    Z base;
    Point origin;
    std::map<Key,Answer> cache;
    std::set<Linear> guards;
    Ratio radius{1,1};

    void impose(const Form &test, const Chart &c) {
        Linear f = pull(test,c);
        if (!guards.insert(f).second) return;
        const W v = Form{f.row,f.b,0}.lead(origin);
        assert(v >= 0);
        Z slope = 0;
        for (Z a : f.row) slope += std::abs(a);
        if (!slope) return;
        assert(v > 0);
        const Ratio limit{v,origin.den*slope};
        if (limit < radius) radius = limit;
    }
    void phase_guards(const Phase &p, const Chart &c) {
        for (const Form &g : p.guards) impose(g,c);
    }
    void certify(const Context &ctx, const Stage &st, const Chart &c) {
        State state = identity();
        for (U w : st.words) {
            const Phase phase = primitive(state,ctx);
            assert(phase.word == w);
            phase_guards(phase,c);
            for (const Form &g : phase.guards) {
                const Z s = -dot(g.row,st.d);
                if (s > 0) impose(st.s*g-s*st.exit,c);
            }
            state = phase.end;
        }
        impose(st.exit,c);
    }

    Answer read(Point u, Lattice l, Chart c, int stage_id, int short_budget) {
        const Key key{l.p,l.cols,c,stage_id,short_budget};
        auto found = cache.find(key);
        if (found != cache.end()) return found->second;
        Answer ans;
        try {
            while (!u.zero() && stage_id < max_stage && short_budget > 0) {
                const Context ctx{u,l,base};
                const Stage st = stage(ctx);
                if (!st.s) {
                    phase_guards(st.first,c);
                    c = advance(st.first.end,c);
                    u = lead(st.first.end,u);
                    l = push(st.first.end,l,base);
                    --short_budget;
                    continue;
                }
                certify(ctx,st,c);
                const Split p = split(l,st.exit,st.s);
                std::map<Next,Z> children;
                for (Z r = 0; r < p.index; ++r) {
                    const auto next = exit_layer(st,p,r,base).next;
                    ++children[{next.p,next.cols}];
                }
                const Point after = projected(u,st);
                const Chart next_chart = project(st,c);
                const Linear distance = pull(st.exit,c);
                const Linear time = Z(st.words.size())*distance/st.s;
                for (const auto &[state,count] : children) {
                    const EdgeType type{stage_id+1,int(st.words.size()),int(l.cols.size()),
                                        int(state.second.size()),st.failures,st.s,p.index,
                                        st.d,distance};
                    const Ratio w{count,p.index};
                    ans.edges[type] = ans.edges[type]+w;
                    const auto downstream = read(after,{state.first,state.second},
                                                 next_chart,stage_id+1,short_budget);
                    for (const auto &[child_clock,m] : downstream.clocks)
                        ans.clocks[time+child_clock] = ans.clocks[time+child_clock]+m*w;
                    for (const auto &[miss,m] : downstream.missed)
                        ans.missed[miss] = ans.missed[miss]+m*w;
                    for (const auto &[edge,m] : downstream.edges)
                        ans.edges[edge] = ans.edges[edge]+m*w;
                }
                cache.emplace(key,ans);
                return ans;
            }
        } catch (const Ambiguous &) {
            ans.missed[2] = {1,1};
            return ans;
        } catch (const ZeroWord &) {
            ans.missed[4] = {1,1};
            return ans;
        }
        if (u.zero()) {
            // Radius zero at one point can be a genuine chart endpoint only
            // when all five leading forms vanish throughout the open cell.
            for (const Linear &f : c) assert((f.row == V{} && f.b == 0));
            ans.clocks[Linear{}] = {1,1};
        } else {
            ans.missed[stage_id >= max_stage ? 3 : 1] = {1,1};
        }
        cache.emplace(key,ans);
        return ans;
    }
};

#ifndef LOTTERY_NO_MAIN
int main(int argc, char **argv) {
    if (argc != 8) {
        std::cerr << "usage: lottery base grid q1 q2 q3 q4 q5\n";
        return 2;
    }
    const Z base = std::stoll(argv[1]), grid = std::stoll(argv[2]);
    V q{};
    for (int i = 0; i < 5; ++i) q[i] = std::stoll(argv[i+3]);
    assert(in_pi(q,grid));
    const Point x = centered(q,grid);
    Ledger ledger{base,x,{},{},{1,1}};
    const Answer a = ledger.read(x,root_lattice(base),root(),0,max_short);
    std::cout << "base " << base << " states " << ledger.cache.size()
              << " inequalities " << ledger.guards.size()
              << " strict-radius " << ledger.radius << '\n';
    Ratio total;
    for (const auto &[c,w] : a.clocks) {
        std::cout << "clock" << c << ' ' << w << '\n';
        total = total+w;
    }
    for (const auto &[kind,w] : a.missed) {
        std::cout << "gap " << kind << ' ' << w << '\n';
        total = total+w;
    }
    assert(total.n == total.d);
    for (const auto &[e,w] : a.edges) {
        std::cout << "edge " << e.id << ' ' << e.n << ' ' << e.rank_in << ' '
                  << e.rank_out << ' ' << e.failures << ' ' << e.s << ' '
                  << e.index << e.d << e.distance << ' ' << w << '\n';
    }
}
#endif

