// Literal, integer subtraction on the two comparison boxes and selected theorem lifts.
#define main settle35_probe_main
#include "probe.cpp"
#undef main

template<int n>
void comparison(Int B, Int width, std::uint64_t seed, Int cap, int admitted, const std::vector<Int>& weights) {
    const Physical<n> map(B);
    std::vector<Int> center(map.m);
    Int total = 0;
    for (int j = 0; j < map.m; ++j) total += weights[map.orbit[j]];
    for (int j = 0; j < map.m; ++j) {
        const Int numerator = B*(1000 + map.m*weights[map.orbit[j]] - total);
        const Int denominator = 1000*map.m;
        assert(numerator % denominator == 0);
        center[j] = numerator/denominator;
    }
    typename Physical<n>::Gap g{};
    std::copy(center.begin(), center.end(), g.begin());
    assert(map.literal(g)==g);
    std::cout << n << '\t' << B << "\t-1\t";
    print_physical(map, g, 2, admitted);
    std::cout << std::endl;
    Random random{seed};
    for (int i = 0; i < 24; ++i) {
        g = map.make_source(center, width, random);
        assert(g.back() > 0 && map.pairing_margin(g) >= 0);
        const auto rows = map.rows(g);
        for (int j = 2; j <= map.p; ++j) if (j!=admitted) assert(rows[j]>=0);
        std::cout << n << '\t' << B << '\t' << i << '\t';
        print_physical(map, g, cap, admitted);
        std::cout << std::endl;
    }
}

void lifts() {
    Int B, width, entry, expected_h, cut_num, cut_den, spread_limit, h_limit, time_limit;
    while (std::cin >> B >> width >> entry >> expected_h >> cut_num >> cut_den
                    >> spread_limit >> h_limit >> time_limit) {
        const Physical<35> map(B);
        Physical<35>::Gap g{};
        for (auto& value : g) std::cin >> value;
        assert(std::accumulate(g.begin(),g.end(),Int{0})==B);
        std::cout << B << '\t' << width << '\t' << entry << '\t' << expected_h << '\t';
        const auto literal = [&](const auto& state) { return map.literal(state); };
        const auto hit = detect(g,literal,4*(entry+66));
        assert(hit.entry==entry && hit.period>=11 && hit.period<=66 && 66%hit.period==0 && hit.period%11==0);
        assert(hit.entry+11<=time_limit);
        Int h = 0, spread = 0, max_radius = 0, min_radius=std::numeric_limits<Int>::max();
        Int least_pair = std::numeric_limits<Int>::max(), least_guard = std::numeric_limits<Int>::max();
        Int first_d2_breach=-1, first_c_return=-1;
        const auto origin = g;
        Root root{};
        for (int j=0;j<r35;++j) root[j]=g[root_slots[j]]-g[0];
        for (Int t = 0; t < hit.entry+hit.period; ++t) {
            assert(g[0]>0 && g.back()>0);
            const Int pairing = map.pairing_margin(g);
            const auto rows = map.rows(g);
            least_pair = std::min(least_pair,pairing);
            for (int j=3;j<=map.p;++j) least_guard=std::min(least_guard,rows[j]);
            assert(pairing>=0 && least_guard>=0);
            assert(root.back()==rows[1] && secondary(root)==rows[2]);
            if (rows[2] < -1 && first_d2_breach<0) first_d2_breach=t;
            if (rows[2] >= -1 && first_d2_breach>=0 && first_c_return<0) first_c_return=t;
            h += loss(rows[2]);
            spread = std::max(spread,root_spread(root));
            if (t>=hit.entry) {
                assert(rows[2]>=-1);
                const Int z = map.peak_numerator(g);
                assert(z*cut_den>map.m*B*cut_num);
                max_radius=std::max(max_radius,z);
                min_radius=std::min(min_radius,z);
            }
            g=literal(g);
            root=root_step(root);
            for (int j=0;j<r35;++j) assert(root[j]==g[root_slots[j]]-g[0]);
        }
        assert(h==expected_h && h<=h_limit && spread<=spread_limit);
        std::cout << hit.entry << '\t' << hit.period << '\t' << spread << '\t' << least_pair << '\t'
                  << least_guard << '\t' << first_d2_breach << '\t' << first_c_return << '\t'
                  << min_radius << '\t' << max_radius;
        for (Int x : origin) std::cout << '\t' << x;
        std::cout << std::endl;
    }
}

int main(int argc,char** argv) {
    assert(argc==2);
    const std::string mode=argv[1];
    if (mode=="97") comparison<97>(256025000,256,970256,262144,2,{0,1});
    else if (mode=="119") comparison<119>(256020000,256,1190256,262144,8,{0,1,1,1,2});
    else if (mode=="lifts") lifts();
    else return 2;
    return 0;
}
