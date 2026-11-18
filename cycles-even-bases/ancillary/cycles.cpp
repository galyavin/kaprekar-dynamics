// Enumerate cycles by four-place positional subtraction. Each (u,v) has the
// representative word (u,v,0,0), whose outer and inner differences are u,v.
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Pair = std::pair<int, int>;
using Cycle = std::vector<Pair>;
using Counts = std::map<int, int>;

Pair step(int base, Pair p) {
    std::array<int, 4> output;
    std::array<int, 4> descending{p.first, p.second, 0, 0};
    int borrow = 0;
    for (int i = 3; i >= 0; --i) {
        int diff = descending[i] - descending[3-i] - borrow;
        borrow = diff < 0;
        output[i] = diff + borrow * base;
    }
    std::sort(output.begin(), output.end());
    return {output[3] - output[0], output[2] - output[1]};
}

int index(Pair p) {
    return p.first * (p.first + 1) / 2 + p.second;
}

std::vector<Cycle> enumerate(int base) {
    std::vector<unsigned char> color(base * (base + 1) / 2);
    std::vector<Cycle> result;
    Cycle path;
    for (int u = 1; u < base; ++u) {
        for (int v = 0; v <= u; ++v) {
            Pair q{u, v};
            if (color[index(q)] != 0) continue;
            path.clear();
            while (color[index(q)] == 0) {
                color[index(q)] = 1;
                path.push_back(q);
                q = step(base, q);
            }
            if (color[index(q)] == 1) {
                auto start = std::find(path.begin(), path.end(), q);
                result.emplace_back(start, path.end());
            }
            for (Pair p : path) color[index(p)] = 2;
        }
    }
    return result;
}

void require(bool condition, int base, const std::string& detail) {
    if (!condition) throw std::runtime_error("base " + std::to_string(base) + ": " + detail);
}

Counts small_core(int m, int n) {
    Counts expected;
    if (m == 3) {
        if (n % 2) expected[3*(n+1)] = 1;
        else {
            expected[n+1] = 1;
            expected[2*(n+1)] = 1;
        }
    } else if (n == 1) {
        expected[1] = 2;
    } else {
        // Count the entries in each row, including zero multiplicities at n=2.
        for (int i = 0; i < (n-1)/2; ++i) {
            ++expected[n];
            ++expected[n+2];
        }
        if (n % 2 == 0) {
            ++expected[n/2];
            ++expected[n/2+1];
        }
    }
    return expected;
}

void check_large(int base, int m, int n, int scale, const std::vector<Cycle>& cycles) {
    std::map<Pair, int> additional;
    if (scale % m == 1) additional[{scale, 0}] = n+1;
    if (m == 5 && n % 4 == 2) additional[{2*scale, 0}] = 2*(n+1);
    int lifted = 0;
    for (const Cycle& cycle : cycles) {
        Pair q = cycle.front();
        bool on_lift = q.first % scale == 0 && q.second % scale == 0 &&
            q.first > q.second && q.first/scale % 2 && q.second/scale % 2;
        if (on_lift) {
            for (size_t i = 0; i < cycle.size(); ++i) {
                auto [a, b] = cycle[i];
                require(a % scale == 0 && b % scale == 0, base, "scale");
                a /= scale;
                b /= scale;
                require(a > b && b > 0 && a < m && a % 2 && b % 2, base, "odd core");
                a = std::abs(2*a-m);
                b = std::abs(2*b-m);
                Pair next{scale*std::max(a,b), scale*std::min(a,b)};
                require(cycle[(i+1) % cycle.size()] == next, base, "directed lift");
            }
            lifted += static_cast<int>(cycle.size());
        } else {
            auto found = std::find_if(additional.begin(), additional.end(), [&](const auto& item) {
                return std::find(cycle.begin(), cycle.end(), item.first) != cycle.end();
            });
            require(found != additional.end(), base, "unexpected cycle");
            require(static_cast<int>(cycle.size()) == found->second, base, "extra cycle period");
            additional.erase(found);
        }
    }
    const int odds = (m-1)/2;
    require(lifted == odds*(odds-1)/2, base, "missing lifted pair");
    require(additional.empty(), base, "missing extra cycle");
}

int main() {
    std::cout << "[\n";
    bool first = true;
    for (int base = 2; base <= 4096; base += 2) {
        if (base > 1000 && (base & (base-1))) continue;
        int m = base, n = 0;
        while (m % 2 == 0) {
            m /= 2;
            ++n;
        }
        auto cycles = enumerate(base);
        Counts periods;
        for (const Cycle& cycle : cycles) ++periods[static_cast<int>(cycle.size())];
        if (m >= 5) check_large(base, m, n, base/m, cycles);
        else require(periods == small_core(m, n), base, "period counts");
        bool unique = ((m == 3 || m == 5) && n % 2) || (m == 7 && n % 3 != 0);
        require(unique == (cycles.size() == 1), base, "uniqueness");
        if (!first) std::cout << ",\n";
        first = false;
        std::cout << "  {\"base\":" << base << ",\"period_counts\":{";
        for (auto it = periods.begin(); it != periods.end(); ++it) {
            if (it != periods.begin()) std::cout << ",";
            std::cout << "\"" << it->first << "\":" << it->second;
        }
        std::cout << "}}";
    }
    std::cout << "\n]\n";
}
