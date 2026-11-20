// Integer evaluation of the same quadrature integrands as entrance.py.
// Verify against the Python evaluator. No floating-point geometry or mass.
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using Integer = __int128_t;
using Row = std::array<Integer, 6>;
using Matrix = std::array<Row, 6>;

void compositions(int total, int i, Row partial, std::vector<Row>& nodes) {
    if (i == 5) {
        partial[i] = 2 * total + 1;
        nodes.push_back(partial);
        return;
    }
    for (int k = 0; k <= total; ++k) {
        partial[i] = 2 * k + 1;
        compositions(total - k, i + 1, partial, nodes);
    }
}

Integer density(const Row& g, Integer unit, bool direct) {
    if (direct) {
        Integer product = 1;
        for (Integer value : g) product *= value;
        return product;
    }
    const std::array<Integer, 5> low{g[0], g[1], g[0] + g[2],
                                     g[1] + g[3], g[0] + g[2] + g[4]};
    Integer total = 0;
    for (unsigned mask = 0; mask < 16; ++mask) {
        Integer last = unit, product = 1;
        for (unsigned i = 0; i < 5; ++i) {
            if (mask >> i & 1) {
                const Integer next = unit - low[i];
                product *= last - next;
                last = next;
            }
        }
        for (int i = 4; i >= 0; --i) {
            if (!(mask >> i & 1)) {
                product *= last - low[i];
                last = low[i];
            }
        }
        total += 2 * product * last;
    }
    return total;
}

std::string decimal(Integer value) {
    if (!value) return "0";
    if (value < 0) return "-" + decimal(-value);
    std::string result;
    while (value) {
        result.insert(result.begin(), char('0' + value % 10));
        value /= 10;
    }
    return result;
}

int main() {
    std::array<std::vector<Row>, 4> rules;
    for (int i = 0; i < 4; ++i) compositions(3 - i, 0, {}, rules[i]);
    int mode;
    int scale;
    while (std::cin >> mode >> scale) {
        // Every input coordinate is a nonnegative gap numerator and each
        // matrix column sums to scale. Each summand is a product of six
        // nonnegative gaps whose sum is at most 12*scale. By AM-GM its
        // maximum is (2*scale)^6. Even 56 nodes times 32 branches fits
        // within signed 128 bits for scale <= 300000.
        if (scale < 1 || scale > 300000 || (mode != 0 && mode != 1)) return 1;
        Matrix matrix;
        for (auto& row : matrix) {
            for (Integer& value : row) {
                int input;
                if (!(std::cin >> input) || input < 0 || input > scale) return 1;
                value = input;
            }
        }
        for (int j = 0; j < 6; ++j) {
            Integer sum = 0;
            for (const Row& row : matrix) sum += row[j];
            if (sum != scale) return 1;
        }
        for (int i = 0; i < 4; ++i) {
            Integer result = 0;
            const Integer unit = scale * (12 - 2 * i);
            for (const Row& node : rules[i]) {
                Row g{};
                for (int j = 0; j < 6; ++j) {
                    for (int k = 0; k < 6; ++k) g[j] += matrix[j][k] * node[k];
                }
                if (!mode) {
                    if (g[0] > g[1] || g[1] > g[0] + g[2]
                        || g[0] + g[2] > g[1] + g[3]
                        || g[1] + g[3] > g[0] + g[2] + g[4]
                        || 2 * (g[0] + g[2] + g[4]) > unit) return 1;
                }
                result += density(g, unit, mode == 1);
            }
            std::cout << (i ? " " : "") << decimal(result);
        }
        std::cout << std::endl;
    }
}
