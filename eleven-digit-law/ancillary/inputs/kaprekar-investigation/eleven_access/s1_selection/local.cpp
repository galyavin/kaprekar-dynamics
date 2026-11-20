#define S1_PHYSICAL_NO_MAIN
#include "physical.cpp"

Word representative(const Q &q) {
    Word word{};
    Q gaps{};
    for (int i=0; i<4; ++i) gaps[i] = q[i]-q[i+1];
    gaps[4] = q[4];
    for (int i=0; i<5; ++i) word[i+1] = word[i]+gaps[i]/2;
    for (int i=5; i<10; ++i) word[i+1] = word[i]+gaps[9-i]-gaps[9-i]/2;
    assert(word_differences(word)==q);
    return word;
}

long double scaled_multiplicity(const Q &q, int base) {
    // The exact 3^5 zero/positive split patterns from the ordered word weight.
    Q gaps{};
    for (int i=0; i<4; ++i) gaps[i] = q[i]-q[i+1];
    gaps[4] = q[4];
    constexpr Z fact[] = {1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800};
    long double sum = 0;
    for (int pattern=0; pattern<243; ++pattern) {
        std::array<bool,10> zero{};
        long double ways = 1;
        int code = pattern;
        for (int i=0; i<5; ++i) {
            if (code%3==0) ways *= gaps[i]-1;
            else zero[code%3==1 ? i : 9-i] = true;
            code /= 3;
        }
        int run = 1;
        Z symmetry = 1;
        for (bool same : zero) {
            if (same) ++run;
            else { symmetry *= fact[run]; run = 1; }
        }
        symmetry *= fact[run];
        sum += ways*(fact[11]/symmetry);
    }
    sum *= base-q[0];
    for (int i=0; i<6; ++i) sum /= base;
    return sum;
}

Q cube_point(const Ranks &ranks, U row, int base) {
    constexpr int center[] = {2938,3225,2886,1913,387};
    Rng retry{indexed_seed(14098603117ULL ^ (U(base)<<32),row)};
    Q q{};
    for (int i=0; i<5; ++i) {
        const Z c = (5-i)*200000+4*center[i];
        constexpr Z den = 1200000;
        const Z lo = (c-1)*base/den+1;
        const Z hi = ((c+1)*base-1)/den;
        assert(lo<=hi);
        const U width = hi-lo+1, threshold = -width%width;
        __uint128_t candidate = __uint128_t(ranks[i])*width;
        while (U(candidate)<threshold) candidate = __uint128_t(retry.next())*width;
        q[i] = int(lo+Z(candidate>>64));
    }
    return q;
}

int main(int argc, char **argv) {
    if (argc!=4) return 2;
    const int anchor = std::stoi(argv[1]);
    const U start = std::stoull(argv[2]), count = std::stoull(argv[3]);
    std::vector<Q> path;
    std::cout.precision(17);
    for (U row=start; row<start+count; ++row) {
        const Ranks ranks = indexed_ranks(14098603103ULL,row);
        std::cout << row;
        for (int residue : {3,5,9,11,15,17,21,23}) {
            const int base = anchor+residue;
            const Q q = cube_point(ranks,row,base);
            const auto a = physical_clock(representative(q),base,path);
            assert(a.result.clock>=0 && a.result.strict==0 && !a.result.strict_loss);
            std::cout << ' ' << a.result.clock << ' ' << scaled_multiplicity(q,base)
                      << ' ' << a.skipped << ' ' << a.fallback;
        }
        std::cout << '\n';
    }
}
