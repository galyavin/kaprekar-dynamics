// Reconstruct the bounded strict-index-two scan using ord_p(4) and F(i).
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

using std::uint16_t;
using std::uint64_t;

std::vector<uint64_t> odd_sieving_primes(uint64_t limit) {
    std::vector<uint64_t> primes;
    for (uint64_t candidate = 3; candidate * candidate < limit; candidate += 2) {
        bool is_prime = true;
        for (uint64_t divisor : primes) {
            if (divisor * divisor > candidate) {
                break;
            }
            if (candidate % divisor == 0) {
                is_prime = false;
                break;
            }
        }
        if (is_prime) {
            primes.push_back(candidate);
        }
    }
    return primes;
}

uint64_t power_of_four(uint64_t exponent, uint64_t p) {
    uint64_t base = 4;
    uint64_t result = 1;
    while (exponent != 0) {
        if (exponent % 2 != 0) {
            result = result * base % p;
        }
        base = base * base % p;
        exponent /= 2;
    }
    return result;
}

bool folded_index_two(uint64_t p, uint64_t m,
                      const std::vector<uint16_t>& smallest_divisor) {
    // For p = 8m+1 prime the order of 4 divides 2m by (2/p)=1.
    // A proper divisor of 2m divides (2m)/q for a prime divisor q of 2m.
    const uint64_t target_order = 2 * m;
    if (power_of_four(m, p) == 1) {
        return false;
    }
    while (m % 2 == 0) {
        m /= 2;
    }
    while (m > 1) {
        const uint64_t stored = smallest_divisor[m / 2];
        const uint64_t divisor = stored == 0 ? m : stored;
        if (power_of_four(target_order / divisor, p) == 1) {
            return false;
        }
        do {
            m /= divisor;
        } while (m % divisor == 0);
    }
    return true;
}

int legendre(uint64_t value, uint64_t prime) {
    // Jacobi reciprocity with prime initial denominator and 0 < value < prime.
    int sign = 1;
    while (value > 1) {
        while (value % 2 == 0) {
            value /= 2;
            if (prime % 8 == 3 || prime % 8 == 5) {
                sign = -sign;
            }
        }
        if (value == 1) {
            break;
        }
        if (value % 4 == 3 && prime % 4 == 3) {
            sign = -sign;
        }
        std::swap(value, prime);
        value %= prime;
    }
    return sign;
}

uint64_t first_strict_obstruction(uint64_t p) {
    int value = 1;  // F(1)
    for (uint64_t i = 2; i <= (p - 1) / 4; ++i) {
        value += legendre(2 * i - 1, p) - legendre(i - 1, p);
        if (value >= 1) {
            return i;
        }
    }
    return 0;
}

int main(int argc, char** argv) {
    char* end = nullptr;
    const uint64_t limit = argc == 2 ? std::strtoull(argv[1], &end, 10) : 0;
    if (limit < 18 || limit > 1000000000 || *end != '\0') {
        std::cerr << "Usage: scan_strict UPPER_EXCLUSIVE (18 to 1000000000)\n";
        return 2;
    }
    const uint64_t max_m = (limit - 2) / 8;
    std::vector<bool> composite(max_m + 1);
    std::vector<uint16_t> smallest_divisor((max_m + 1) / 2);
    for (uint64_t q : odd_sieving_primes(limit)) {
        // q^2 is 1 mod 8. This marks its multiples in the progression,
        // starting at q^2 so the prime q itself is preserved.
        for (uint64_t m = (q * q - 1) / 8; m <= max_m; m += q) {
            composite[m] = true;
        }
        for (uint64_t n = q * q; n <= max_m; n += 2 * q) {
            if (smallest_divisor[n / 2] == 0) {
                smallest_divisor[n / 2] = static_cast<uint16_t>(q);
            }
        }
    }

    uint64_t primes_in_progression = 0;
    uint64_t index_two = 0;
    uint64_t last_deciding_row = 0;
    uint64_t first_at_last_deciding_row = 0;
    for (uint64_t m = 1; m <= max_m; ++m) {
        if (composite[m]) {
            continue;
        }
        const uint64_t p = 8 * m + 1;
        ++primes_in_progression;
        if (!folded_index_two(p, m, smallest_divisor)) {
            continue;
        }
        ++index_two;
        const uint64_t row = first_strict_obstruction(p);
        if (row == 0) {
            std::cout << "strict " << p << '\n';
        }
        if (row > last_deciding_row) {
            last_deciding_row = row;
            first_at_last_deciding_row = p;
        }
    }
    std::cout << "p < " << limit << ", p = 1 mod 8 primes "
              << primes_in_progression << ", folded index two " << index_two << '\n';
    std::cout << "largest first obstruction " << last_deciding_row
              << ", first prime attaining it " << first_at_last_deciding_row << '\n';
}
