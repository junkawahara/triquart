// count_type.hpp — the integer type used for triplet / quartet totals.
//
// Counters inside a component are bounded by C(n,3) and fit in 64 bits for
// any realistic n, but the quartet totals of Fig. 3 are up to C(n,4), which
// exceeds 2^63 already for n ≈ 123 000.  All contributions and totals are
// therefore accumulated in a 128-bit integer when the compiler offers one.
#pragma once

#include <cstdint>
#include <string>

namespace soda13 {

#if defined(__SIZEOF_INT128__)
using count_t = __int128;
#else
using count_t = long long;  // quartet totals overflow for n > ~120 000
#endif

inline std::string count_to_string(count_t v) {
    if (v == 0) return "0";
    const bool neg = v < 0;
#if defined(__SIZEOF_INT128__)
    using ucount_t = unsigned __int128;
#else
    using ucount_t = unsigned long long;
#endif
    ucount_t u = neg ? static_cast<ucount_t>(-(v + 1)) + 1 : static_cast<ucount_t>(v);
    std::string s;
    while (u) {
        s += static_cast<char>('0' + static_cast<int>(u % 10));
        u /= 10;
    }
    if (neg) s += '-';
    return std::string(s.rbegin(), s.rend());
}

}  // namespace soda13
