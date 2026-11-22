// Search for a flag of projected returns under the exact T_35 recurrence.
// The search's modular span decisions are independently checked over Q by the certificate verifier.
#define main settle35_probe_main
#include "probe.cpp"
#undef main
#include <map>

constexpr int R = r35, Forms = 2*R;
using Vector = std::array<Int,R>;
using Bits = std::uint32_t;

constexpr int wrap(int i) { return i >= R ? i-R : i; }
bool bit(Bits mask, int f) { return (mask >> f) & 1u; }

Int evaluate(int f, const Vector& v) {
    if (f < R) return v[f];
    const int i = f-R;
    return v[i] + v[wrap(i+7)] - v[wrap(i+11)];
}

struct State {
    Vector offsets{};
    Int gauge = 0;
    int phase = 0;

    Vector values() const {
        Vector v = offsets;
        for (auto& i : v) i += gauge;
        return v;
    }

    Int read(int f) const { return evaluate(f, offsets) + gauge; }

    void step(Bits mask) {
        const Int y = bit(mask,phase) ? read(phase) : 1;
        const Int h1 = loss(y);
        const Int h2 = bit(mask,R+phase) ? loss(read(R+phase)) : 0;
        const int e = sign(y);
        offsets[phase] += h1-2*e;
        offsets[wrap(phase+11)] += h1+e-h2;
        offsets[wrap(phase+10)] += h2;
        offsets[wrap(phase+7)] += h2-h1;
        gauge += h1-e;
        phase = wrap(phase+1);
    }

    void shift(const Vector& ray, Int count) {
        for (int t=0; t<R; ++t) offsets[t] += count*ray[t];
    }
};

struct Span {
    static constexpr int prime = 1009;
    std::array<Vector,R> reduced{};
    Bits pivots = 0;

    Vector residual(int form) const {
        Vector row{};
        const int i = form % R;
        row[i] = 1;
        if (form >= R) {
            row[wrap(i+7)] = 1;
            row[wrap(i+11)] = -1;
        }
        for (int j=0; j<R; ++j) {
            if (!bit(pivots,j) || !row[j]) continue;
            const Int value = row[j];
            for (int k=j; k<R; ++k) {
                row[k] = (row[k]-value*reduced[j][k]) % prime;
            }
        }
        return row;
    }

    static Int inverse(Int x) {
        Int answer = 1;
        for (int e=prime-2; e; e/=2, x=x*x%prime) if (e%2) answer=answer*x%prime;
        return answer;
    }

    void insert(int form) {
        auto row = residual(form);
        for (int j=0; j<R; ++j) {
            if (!row[j]) continue;
            assert(!bit(pivots,j));
            const Int inv = inverse(row[j]);
            for (auto& cell : row) cell=cell*inv%prime;
            reduced[j]=row;
            pivots |= 1u<<j;
            return;
        }
        assert(false);
    }

    Bits closure() const {
        Bits mask = 0;
        for (int f=0; f<Forms; ++f) if (residual(f)==Vector{}) mask |= 1u<<f;
        return mask;
    }
};


// The recession cone keeps every forced primary read positive and every
// forced secondary read inactive. Intersect the nonnegative orthant with
// the twelve secondary halfspaces, keeping only genuine extremal rays.
std::vector<Bits> cone_rays() {
    struct Ray { Vector coordinates; Bits zeros; };
    std::vector<Ray> rays;
    constexpr Bits primary = (1u<<R)-1;
    for (int j=0; j<R; ++j) {
        Vector x{};
        x[j]=1;
        rays.push_back({x,primary ^ (1u<<j)});
    }
    for (int f=R; f<Forms; ++f) {
        std::vector<Ray> next;
        for (const auto& a : rays) {
            const Int A=evaluate(f,a.coordinates);
            if (A>=0) next.push_back({a.coordinates,a.zeros | (A==0 ? 1u<<f : 0)});
            if (A<=0) continue;
            for (const auto& b : rays) {
                const Int B=evaluate(f,b.coordinates);
                if (B>=0) continue;
                const Bits common=a.zeros & b.zeros;
                bool adjacent=true;
                for (const auto& c : rays) {
                    if (&c==&a || &c==&b) continue;
                    if ((c.zeros & common)==common) { adjacent=false; break; }
                }
                if (!adjacent) continue;
                Vector x;
                Int divisor=0;
                for (int j=0; j<R; ++j) {
                    x[j]=A*b.coordinates[j]-B*a.coordinates[j];
                    divisor=std::gcd(divisor,x[j]);
                }
                for (Int& v : x) v/=divisor;
                next.push_back({x,common | (1u<<f)});
            }
        }
        rays=std::move(next);
    }
    std::vector<Bits> zeros;
    for (const auto& ray : rays) zeros.push_back(ray.zeros);
    return zeros;
}

struct Faces {
    std::vector<Bits> cones=cone_rays();
    std::map<Bits,bool> cache;
    bool includes(Bits mask) {
        if (const auto pos=cache.find(mask); pos!=cache.end()) return pos->second;
        constexpr Bits all=(1u<<Forms)-1;
        Bits common=all;
        for (Bits ray : cones) {
            if ((mask & ray) != mask) continue;
            common &= ray;
            if (common==mask) break;
        }
        return cache.emplace(mask,common==mask).first->second;
    }
};

struct Projection {
    Int prefix = 0, period = 0;
    State start, end;
    Vector ray{};
    std::array<Int,Forms> minima{};
};

Projection project(const State& s, Bits tracked, Int cap) {
    auto same = [&](const State& a, const State& b) {
        for (Bits pending=tracked; pending; pending &= pending-1) {
            int f=__builtin_ctz(pending);
            if (a.read(f)!=b.read(f)) return false;
        }
        return true;
    };
    auto turn = [&](State t) {
        for (int j=0; j<R; ++j) t.step(tracked);
        return t;
    };
    Projection p;
    State a=s, b=turn(s);
    Int power=1, period=1, used=1;
    while (!same(a,b) && used<cap) {
        if (power==period) { a=b; power*=2; period=0; }
        b=turn(b);
        ++period; ++used;
        // Prune divergent projections comfortably before the arithmetic range.
        if (std::abs(b.gauge) >= (Int{1}<<38)) return p;
        for (Int x : b.offsets) if (std::abs(x)>=(Int{1}<<38)) return p;
    }
    if (!same(a,b)) return p;
    a=s; b=s;
    for (Int j=0; j<period; ++j) b=turn(b);
    while (!same(a,b)) { a=turn(a); b=turn(b); ++p.prefix; }
    p.prefix *= R;
    p.period = period*R;
    p.start=a;
    p.end=a;
    p.minima.fill(std::numeric_limits<Int>::max());
    for (Int t=0; t<p.period; ++t) {
        const int phase=p.end.phase;
        for (int f : {phase,R+phase}) p.minima[f]=std::min(p.minima[f],p.end.read(f));
        p.end.step(tracked);
    }
    assert(same(p.start,p.end));
    const auto start=p.start.values(), end=p.end.values();
    for (int j=0; j<R; ++j) p.ray[j]=start[j]-end[j];
    return p;
}

struct Stage {
    int added;
    Int prefix, period;
    Vector shift;
    Vector ray;
};

Int ceildiv(Int a, Int b) {
    assert(b>0);
    return a>=0 ? (a+b-1)/b : a/b;
}

int main(int argc, char** argv) {
    assert(argc==5 || argc==6);
    const Int limit=std::stoll(argv[1]), cap=std::stoll(argv[2]);
    const Int seed=std::stoll(argv[3]), terminal_cap=std::stoll(argv[4]);
    Random random{static_cast<std::uint64_t>(seed)};
    int longest=0;
    Faces faces;
    std::cerr << "recession extreme rays " << faces.cones.size() << std::endl;
    for (Int attempt=0; attempt<limit; ++attempt) {
        State state;
        Span span;
        Bits mask=0;
        std::vector<Stage> stages;
        for (int rank=0; rank<R; ++rank) {
            const auto p=project(state,mask,rank>=11 ? 32*cap : (rank>=10 ? 4*cap : cap));
            if (!p.period) break;
            Int common=0;
            for (Int x : p.ray) common=std::gcd(common,std::abs(x));
            if (!common) break;
            Vector direction=p.ray;
            for (Int& x : direction) x/=common;
            struct Choice { int form; Span span; Bits mask; Int shift; };
            std::vector<Choice> choices;
            for (int f=0; f<Forms; ++f) {
                if (bit(mask,f) || evaluate(f,direction)<=0) continue;
                auto next=span;
                next.insert(f);
                const Bits enlarged=next.closure();
                if (!faces.includes(enlarged)) continue;
                const Bits fresh=enlarged & ~mask;
                bool admissible=true;
                Int shift=-(Int{1}<<62);
                for (int g=0; g<Forms; ++g) {
                    if (!bit(fresh,g)) continue;
                    const Int slope=evaluate(g,direction);
                    if (slope<=0) { admissible=false; break; }
                    shift=std::max(shift,ceildiv((g<R ? 1 : -1)-p.minima[g],slope));
                }
                if (admissible) choices.push_back({f,next,enlarged,shift});
            }
            if (choices.empty()) break;
            const auto& choice=choices[random.next()%choices.size()];
            Int shift=choice.shift;
            if (argc==6 && rank>0) shift += random.next()%std::stoll(argv[5]);
            Vector translation=direction;
            for (Int& x : translation) x*=shift;
            state=p.end;
            state.shift(direction,shift);
            stages.push_back({choice.form,p.prefix,p.period,translation,p.ray});
            mask=choice.mask;
            span=choice.span;
        }
        if (static_cast<int>(stages.size())>longest) {
            longest=stages.size();
            std::cerr << "max rank " << longest << " attempt " << attempt << ':';
            for (const auto& stage : stages) std::cerr << ' ' << stage.added;
            std::cerr << std::endl;
        }
        if (stages.size()!=R) continue;
        const auto v=state.values();
        Root terminal_start;
        for (int i=0; i<R; ++i) terminal_start[R-1-i]=v[i];
        const auto cycle=detect(terminal_start,root_step,terminal_cap);
        bool zero=cycle.period>0;
        auto runner=terminal_start;
        if (zero) {
            for (Int i=0; i<cycle.entry; ++i) runner=root_step(runner);
            for (Int i=0; i<cycle.period; ++i) {
                if (loss(secondary(runner))) zero=false;
                runner=root_step(runner);
            }
        }
        if (!zero) {
            if (attempt%64==0) std::cerr << "full rank, rejected terminal " << attempt << std::endl;
            continue;
        }
        std::cout << "found seed " << seed << " attempt " << attempt << " entry " << cycle.entry
                  << " period " << cycle.period << '\n';
        for (const auto& stage : stages) {
            std::cout << stage.added << ' ' << stage.prefix << ' ' << stage.period;
            for (Int x : stage.shift) std::cout << ' ' << x;
            std::cout << '\n';
        }
        return 0;
    }
    std::cerr << "finished " << limit << " attempts, maximum rank " << longest << '\n';
    return 1;
}
