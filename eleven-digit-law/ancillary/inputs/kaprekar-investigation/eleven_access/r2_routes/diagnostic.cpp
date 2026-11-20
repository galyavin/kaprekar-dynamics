#define main task097_main
#include "graph.cpp"
#undef main

using Q = std::array<int, 5>;
using Word = std::array<int, 11>;
using Ranks = std::array<U, 11>;
constexpr Z grid = Z(1) << 38;
constexpr U source_seed = 13600314511317ULL;
constexpr U lane_seed = 13600318531319ULL;

U indexed_seed(U seed, U index) {
    Rng rng{seed + index * 0x9e3779b97f4a7c15ULL};
    return rng.next();
}

Q word_differences(Word word) {
    std::sort(word.begin(), word.end());
    Q q{};
    for (int i = 0; i < 5; ++i) q[i] = word[10-i]-word[i];
    return q;
}

Word schoolbook(const Q &q, int b) {
    Word w{};
    int borrow = 0;
    for (int i = 10; i >= 0; --i) {
        int v = (i<5 ? q[i] : i>5 ? -q[10-i] : 0)-borrow;
        borrow = v < 0;
        w[i] = v + borrow*b;
    }
    assert(!borrow);
    return w;
}

struct PhysicalPhase {
    Q next;
    int loss; // Bit 1 denotes interlevel inversion, bit 2 zero innermost.
    bool regular;
};

PhysicalPhase physical_phase(const Q &q, int b) {
    int height = 5;
    while (height && !q[height-1]) --height;
    if (!height) return {{}, 2, false};
    Q left{}, right{};
    for (int i = 0; i < height; ++i) {
        left[i] = q[i]-(i==height-1);
        right[i] = b-q[height-1-i]-(i<height-1);
    }
    // The last right digit carries one less borrow. Move it ahead of the
    // preceding tied q[0] entries if necessary, keeping the original levels
    // in right for the physical pairing reader.
    Q sorted_right = right;
    for (int t = height-1; t>0 && sorted_right[t]>sorted_right[t-1]; --t)
        std::swap(sorted_right[t],sorted_right[t-1]);
    // The 11-2*height central digits are all b-1.
    Word digit{};
    int k = 11-2*height;
    for (int i = 0; i < k; ++i) digit[i] = b-1;
    int i = 0, j = 0;
    while (i < height && j < height)
        digit[k++] = left[i] >= sorted_right[j] ? left[i++] : sorted_right[j++];
    while (i < height) digit[k++] = left[i++];
    while (j < height) digit[k++] = sorted_right[j++];
    PhysicalPhase out{{}, height<5 ? 2 : 0, height==5};
    for (int t = 0; t < 5; ++t) out.next[t] = digit[t]-digit[10-t];
    if (height==5) {
        for (int t = 0; t < 4; ++t)
            if (std::min(left[t],right[t]) < std::max(left[t+1],right[t+1])) out.loss = 1;
        for (int t = 0; t < 5; ++t) out.regular &= left[t] >= right[t];
        out.regular &= !out.loss;
    }
    return out;
}

struct Physical {
    int clock = -1, strict = -1, strict_loss = 0, five = -1, five_loss = 0;
    int run = 0;

    void read(const Q &q, const PhysicalPhase &phase, int t, int base) {
        if (strict < 0 && in_pi({q[0],q[1],q[2],q[3],q[4]},base)) strict = t;
        if (strict >= 0) strict_loss |= phase.loss;
        if (five >= 0) five_loss |= phase.loss;
        else {
            run = phase.regular ? run+1 : 0;
            if (run==5) five = t-4;
        }
    }
};

Physical observe(const Word &word, int base, std::vector<Q> &path) {
    Q q = word_differences(word);
    Physical result;
    path.clear();
    int checkpoint = 0, power = 1, period = 0;
    for (int t = 0; t < 48*base+10000; ++t) {
        if (t && q==path[checkpoint]) {
            period = t-checkpoint;
            break;
        }
        if (t-checkpoint==power) {
            checkpoint = t;
            power *= 2;
        }
        path.push_back(q);
        const auto phase = physical_phase(q,base);
        result.read(q,phase,t,base);
        q = phase.next;
    }
    if (!period) return result;
    int lo = 0, hi = checkpoint;
    while (lo < hi) {
        const int mid = (lo+hi)/2;
        if (path[mid]==path[mid+period]) hi = mid;
        else lo = mid+1;
    }
    const Word preceding = lo ? schoolbook(path[lo-1],base) : word;
    const Word cyclic = schoolbook(path[lo+period-1],base);
    result.clock = lo+(preceding!=cyclic);
    const int end = int(path.size())+2*period+5;
    for (int t = int(path.size()); t < end; ++t) {
        const auto phase = physical_phase(q,base);
        result.read(q,phase,t,base);
        q = phase.next;
    }
    return result;
}

Word word_at(const Ranks &ranks, U row, int base) {
    Rng retry{indexed_seed(source_seed ^ (U(base)<<32),row)};
    const U threshold = -U(base) % U(base);
    Word w{};
    auto sample = [&](U first) {
        __uint128_t prod = __uint128_t(first)*base;
        while (U(prod)<threshold) prod = __uint128_t(retry.next())*base;
        return int(prod>>64);
    };
    for (int i = 0; i < 11; ++i) w[i] = sample(ranks[i]);
    while (std::all_of(w.begin(),w.end(),[&](int digit){return digit==w[0];}))
        for (int &digit : w) digit = sample(retry.next());
    return w;
}

#ifndef R2_NO_MAIN
int main(int argc, char **argv) {
    if (argc!=3) return 2;
    const int first = std::stoi(argv[1]), count = std::stoi(argv[2]);
    std::vector<Q> path;
    std::cout.precision(17);
    for (int row = first; row < first+count; ++row) {
        Rng source{indexed_seed(source_seed,row)};
        Ranks ranks{};
        Digits continuous{};
        for (int i = 0; i < 11; ++i) {
            ranks[i] = source.next();
            continuous[i] = Z(ranks[i]>>26);
        }
        V q = differences(continuous);
        const int leading_entry = entry(q,grid,500);
        std::cout << row << ' ' << leading_entry;
        for (int scale : {12000,48000,192000}) {
            for (int r : {3,5,15,17}) {
                const int base = scale+r;
                Trail graph;
                graph.status = 5;
                if (leading_entry >= 0) {
                    Rng lanes{indexed_seed(lane_seed,row)};
                    graph = trace(centered(q,grid),root_lattice(base),base,lanes,false);
                }
                const auto p = observe(word_at(ranks,row,base),base,path);
                std::cout << ' ' << graph.status << ' ' << (graph.count.first_loss>=0)
                          << ' ' << static_cast<double>(graph.count.elapsed.n)/graph.count.elapsed.d
                          << ' ' << p.clock << ' ' << p.strict << ' ' << p.strict_loss
                          << ' ' << p.five << ' ' << p.five_loss;
            }
        }
        std::cout << '\n';
    }
}
#endif
