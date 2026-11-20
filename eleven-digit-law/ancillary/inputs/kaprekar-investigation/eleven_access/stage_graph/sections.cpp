// Read all first-exit lanes at a generic point of one exact first winner.
#define LOTTERY_NO_MAIN
#include "lottery.cpp"

void report(const Form &g, const Chart &c, const Lattice &l, Z base, Z decrement, int kind) {
    const Linear f = pull(g,c);
    Z stride = decrement;
    for (const V &col : l.cols) stride = std::gcd(stride,dot(g.row,col));
    const Z offset = g.eval(l.p,base);
    std::cout << "guard" << f << ' ' << g.b << ' ' << decrement << ' '
              << offset << ' ' << stride << ' ' << kind << '\n';
}

int main() {
    Z base=0, den=0;
    V q{};
    while (std::cin >> base >> den >> q[0] >> q[1] >> q[2] >> q[3] >> q[4]) {
        assert(in_pi(q,den));
        const Point x = centered(q,den);
        const Lattice start = root_lattice(base);
        const Stage head = stage({x,start,base});
        assert(head.s==36 && head.words.size()==5 && head.failures==0);
        const Split first = split(start,head.exit,head.s);
        std::cout << "face " << base << ' ' << first.index << '\n';
        for (Z j=0; j<first.index; ++j) {
            const Edge edge = exit_layer(head,first,j,base);
            std::cout << "remainder " << edge.remainder << '\n';
            try {
                Point u = projected(x,head);
                Lattice l = edge.next;
                Chart chart = project(head,root());
                int steps = 0;
                for (; steps<20; ++steps) {
                    const Context ctx{u,l,base};
                    const Stage next = stage(ctx);
                    if (next.s) {
                        std::cout << "stage " << steps << ' ' << next.words.size()
                                  << ' ' << next.failures << next.d << '\n';
                        State f = identity();
                        for (U w : next.words) {
                            const Phase p = primitive(f,ctx);
                            assert(p.word==w);
                            for (std::size_t k=0; k<p.guards.size(); ++k) {
                                const Form &g = p.guards[k];
                                assert(p.height==5 && p.paired && p.guards.size()==15);
                                report(g,chart,l,base,std::max<Z>(0,-dot(g.row,next.d)),
                                       k == 0 ? 0 : k <= 10 ? 1 : 2);
                            }
                            f = p.end;
                        }
                        break;
                    }
                    for (const Form &g : next.first.guards) report(g,chart,l,base,0,1);
                    chart = advance(next.first.end,chart);
                    l = push(next.first.end,l,base);
                    u = lead(next.first.end,u);
                }
                if (steps==20) std::cout << "gap budget\n";
            } catch (const Ambiguous &) {
                std::cout << "gap tie\n";
            }
        }
    }
}

