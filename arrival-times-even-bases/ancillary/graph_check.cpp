// Enumerate the functional graph by subtracting ascending from descending
// base-B words. Each pair (u,v) is represented by the digits (u,v,0,0).
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

struct Pair { int u, v; };

int index_of(Pair p) { return p.u * (p.u + 1) / 2 + p.v - 1; }

Pair advance(int base, Pair p) {
    std::array<int,4> digits{0, 0, p.v, p.u};
    std::int64_t high = 0, low = 0;
    for (int i = 0; i < 4; ++i) {
        low = base * low + digits[i];
        high = base * high + digits[3-i];
    }
    auto difference = high - low;
    for (int &digit : digits) {
        digit = difference % base;
        difference /= base;
    }
    std::sort(digits.begin(), digits.end());
    return {digits[3] - digits[0], digits[2] - digits[1]};
}

std::int64_t weight(int base, Pair p) {
    const std::int64_t starts = base - p.u;
    if (p.v == 0) return starts * (12*p.u - 4);
    if (p.v == p.u) return starts * 6;
    return starts * 24 * (p.u - p.v);
}

int main(int argc, char **argv) {
    if (argc != 3) throw std::runtime_error("Expected an odd part and an exponent");
    const int m = std::stoi(argv[1]), n = std::stoi(argv[2]);
    if (m < 5 || m > 53 || m % 2 == 0 || n < 0 || n > 6)
        throw std::runtime_error("Parameters outside the enumeration range");
    const int scale = 1 << n, base = m*scale;
    const int count = base*(base+1)/2-1;
    std::vector<int> depth(count,-1), cycle(count,-1), path;
    std::map<int,int> maxima;
    std::map<int,std::int64_t> histogram;
    for (int u = 1; u < base; ++u) for (int v = 0; v <= u; ++v) {
        Pair p{u,v};
        const int start = index_of(p);
        int k = start;
        path.clear();
        while (depth[k] == -1) {
            depth[k] = -2;
            path.push_back(k);
            p = advance(base,p);
            k = index_of(p);
        }
        if (depth[k] == -2) {
            auto begin = std::find(path.begin(), path.end(), k);
            const int label = *std::min_element(begin, path.end());
            const int length = path.end() - begin;
            histogram[0] += length;
            histogram[1] -= length;
            for (auto it = begin; it != path.end(); ++it) {
                cycle[*it] = label;
                depth[*it] = 0;
            }
            path.erase(begin,path.end());
        }
        for (auto it = path.rbegin(); it != path.rend(); ++it) {
            cycle[*it] = cycle[k];
            depth[*it] = 1 + depth[k];
            k = *it;
        }
        const int t = 1+depth[start];
        histogram[t] += weight(base,{u,v});
        maxima[cycle[start]] = std::max(maxima[cycle[start]],t);
    }
    for (auto [id,t] : maxima) std::cout << "cycle " << id << ' ' << t << '\n';
    for (auto [t,w] : histogram) std::cout << "weight " << t << ' ' << w << '\n';
    auto report = [&](char kind,int u,int v) {
        int k = index_of({u,v});
        std::cout << kind << ' ' << u << ' ' << cycle[k] << ' ' << depth[k] << '\n';
    };
    if (n) {
        report('Z',2*scale,0);
        report('A',(m+1)*scale/2,(m-1)*scale/2);
        for (int t=1; t<=n; ++t) report('Z',1<<t,0);
        for (int t=0; t<=n-2; ++t) report('A',base/2+(1<<t),base/2-(1<<t));
    }
}
