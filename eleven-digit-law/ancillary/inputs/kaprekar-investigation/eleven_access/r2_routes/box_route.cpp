#define R2_NO_MAIN
#include "diagnostic.cpp"

bool contains(const Lattice &l, const V &physical) {
    V delta = mix(1,physical,-1,l.p);
    for (const V &col : l.cols) {
        int pivot = 0;
        while (!col[pivot]) ++pivot;
        if (delta[pivot]%col[pivot]) return false;
        delta = mix(1,delta,-delta[pivot]/col[pivot],col);
    }
    return delta==V{};
}

Q physical_differences(int base, const V &x) {
    Q q{};
    for (int i = 0; i < 5; ++i) {
        assert((x[i]+(5-i)*base)%6==0);
        q[i] = int((x[i]+(5-i)*base)/6);
    }
    return q;
}

void read_box(int base, int checkpoint, const V &q, int true_clock, Z certified_N) {
    assert(in_pi(q,base));
    Point x = centered(q,base);
    Lattice l = root_lattice(base);
    V actual{};
    for (int i = 0; i < 5; ++i) actual[i] = 6*q[i]-(5-i)*base;
    Z t = 0, failures = 0, macro_failures = 0;
    std::map<int,Z> failure_returns;
    int stages = 0, shorts = 0;
    Ratio leading;
    while (!x.zero()) {
        assert(contains(l,actual));
        const Context ctx{x,l,base};
        const Stage st = stage(ctx);
        if (!st.s) {
            for (const auto &g : st.first.guards) assert(g.eval(actual,base)>=0);
            const auto expected = physical_phase(physical_differences(base,actual),base);
            actual = eval(st.first.end,actual,base);
            assert(physical_differences(base,actual)==expected.next);
            assert(st.first.paired==(expected.loss==0));
            failures += expected.loss!=0;
            x = lead(st.first.end,x);
            l = push(st.first.end,l,base);
            ++t;
            assert(++shorts<=max_short);
            continue;
        }

        const Z a = st.exit.eval(actual,base);
        assert(a>=0);
        const Z n = 1+a/st.s;
        const auto splits = split(l,st.exit,st.s);
        Edge edge{};
        int selected = 0;
        for (Z k = 0; k < splits.index; ++k) {
            auto candidate = exit_layer(st,splits,k,base);
            if (candidate.remainder==a%st.s) {
                edge = std::move(candidate);
                ++selected;
            }
        }
        assert(selected==1);
        for (Z repeat : {Z(0),n-1}) {
            const V start = mix(1,actual,repeat,st.d);
            V qstate = start;
            State f = identity();
            for (const U word : st.words) {
                const Phase phase = primitive(f,ctx);
                assert(phase.word==word);
                for (const auto &g : phase.guards) assert(g.eval(start,base)>=0);
                const auto expected = physical_phase(physical_differences(base,qstate),base);
                qstate = eval(phase.end,start,base);
                assert(phase.paired==(expected.loss==0));
                assert(physical_differences(base,qstate)==expected.next);
                f = phase.end;
            }
            assert(qstate==mix(1,start,1,st.d));
        }
        failures += n*st.failures;
        macro_failures += n*st.failures;
        if (st.failures>0) failure_returns[int(st.words.size())] += n;
        leading = leading + Ratio(W(st.words.size())*st.exit.lead(x),x.den*st.s);
        t += n*st.words.size();
        actual = mix(1,actual,n,st.d);
        assert(st.exit.eval(actual,base)<0);
        x = projected(x,st);
        l = std::move(edge.next);
        assert(++stages<=max_stage);
    }
    assert(contains(l,actual));
    assert(macro_failures>0);
    assert(2*certified_N<=failures && failures<=2*certified_N+21);
    assert(std::abs(t-true_clock)<10000);
    std::cout << base << ' ' << checkpoint << ' ' << true_clock << ' ' << t << ' '
              << leading << ' ' << stages << ' ' << shorts << ' ' << failures << ' '
              << macro_failures << ' ' << certified_N << actual;
    for (auto [length,count] : failure_returns) std::cout << ' ' << length << ':' << count;
    std::cout << '\n';
}

int main() {
    int base, checkpoint, true_clock;
    Z certified_N;
    V q;
    while (std::cin >> base >> checkpoint >> true_clock >> certified_N >> q[0] >> q[1] >> q[2] >> q[3] >> q[4])
        read_box(base,checkpoint,q,true_clock,certified_N);
}
