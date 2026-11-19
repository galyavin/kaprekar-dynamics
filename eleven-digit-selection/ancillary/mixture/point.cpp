// Trace the certificate at a specified centered point. The independent control
// uses the repository's schoolbook whole-word implementation, not this map.
#define main exhaustive_route_main
#include "verify.cpp"
#undef main

int main() {
    try {
        std::array<RouteData,6> routes{RouteData(0),RouteData(1),RouteData(2),
                                       RouteData(3),RouteData(4),RouteData(5)};
        Z base;
        V gaps{};
        while (std::cin >> base >> gaps[0] >> gaps[1] >> gaps[2] >> gaps[3] >> gaps[4]) {
            State state{};
            V x{};
            Z sum = 0;
            for (int i = 0; i < 5; ++i) {
                V row{};
                row[i] = 1;
                state[i] = Form(row,0,1);
                sum += gaps[i];
                x[i] = (i+1)*base - 6*sum;
            }
            auto &data = routes[base%6];
            Form elapsed = 0;
            Lane lane;
            std::vector<Z> residues;
            for (int depth = 0; depth < block_limit; ++depth) {
                const Block *block = nullptr;
                for (int j = 0; j < short_limit; ++j) {
                    block = &macro(state);
                    if (block->divisor) break;
                    state = step(state,block->word[0]);
                    elapsed = elapsed + 1;
                }
                require(block && block->divisor, "missing substantial repeat");
                data.classify(state,*block,depth,elapsed,0,lane);
                Z r = mod(block->numerator.integer(x),block->divisor);
                residues.push_back(r);
                Form count = (block->numerator-r)/block->divisor + 1;
                elapsed = elapsed + Z(block->word.size())*count;
                for (int i = 0; i < 5; ++i) state[i] = state[i] + block->drift[i]*count;
                if (constant(state)) {
                    V spot = mix(1,plain(state),-early_returns,block->drift);
                    Z time = elapsed.integer(x) + data.tail(spot) + 1
                             - early_returns*Z(block->word.size());
                    std::cout << time << lane.sheet << lane.ray << elapsed.leading() << x;
                    for (auto rem : residues) std::cout << ' ' << rem;
                    std::cout << '\n';
                    break;
                }
                require(depth+1 < block_limit, "route block cap");
            }
        }
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
