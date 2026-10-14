#pragma once

// Exact, checked arithmetic for the counters in Sections 3 and 7 of soda13.pdf.
// Overflow is an error, including overflow of intermediate counter expressions.
#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace soda2013 {

using Count = std::uint64_t;

inline Count add(Count a, Count b) {
    if (b > std::numeric_limits<Count>::max() - a)
        throw std::overflow_error("tree counter addition overflow");
    return a + b;
}

inline Count subtract(Count a, Count b) {
    if (b > a) throw std::logic_error("negative tree counter");
    return a - b;
}

inline Count multiply(Count a, Count b) {
    if (a && b > std::numeric_limits<Count>::max() / a)
        throw std::overflow_error("tree counter multiplication overflow");
    return a * b;
}

inline Count choose(Count n, Count k) {
    if (k > n) return 0;
    k = std::min(k, n - k);
    Count result = 1;
    for (Count i = 1; i <= k; ++i) {
        Count numerator = n - k + i;
        Count denominator = i;
        const Count g = std::gcd(numerator, denominator);
        numerator /= g;
        denominator /= g;
        // C(n-k+i-1,i-1) * (n-k+i) / i is integral.
        if (result % denominator != 0)
            throw std::logic_error("nonintegral binomial coefficient");
        result = multiply(result / denominator, numerator);
    }
    return result;
}

}  // namespace soda2013
