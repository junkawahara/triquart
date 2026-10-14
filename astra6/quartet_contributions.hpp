#pragma once

// soda13.pdf, Figures 8 and 13: the six lower-triangular cases in one colored
// traversal; reversing T1 and T2 supplies the other three cases exactly once.
#include "quartet_shared_counters.hpp"

namespace soda2013 {

enum class QuartetType : std::size_t { alpha, beta, gamma };

struct QuartetCases {
    std::array<Count, 9> values{};  // T1 type is the row, T2 type the column.
    Count& at(QuartetType first, QuartetType second) {
        return values[3 * static_cast<std::size_t>(first) + static_cast<std::size_t>(second)];
    }
    Count at(QuartetType first, QuartetType second) const {
        return values[3 * static_cast<std::size_t>(first) + static_cast<std::size_t>(second)];
    }
    void accumulate(const QuartetCases& other);
    Count total() const;
};

// Complete ALL nine cases from complete forward and reverse Count traversals.
// The diagonal must agree between the traversals; it is included only once.
QuartetCases complete_quartet_cases(const QuartetCases& forward, const QuartetCases& reverse);

// Figure 13 contributions anchored at the current composition, excluding the
// children's accumulated totals. CC inputs are LOWER C1 then UPPER C2.
QuartetCases shared_quartet_cc(const QuartetSharedCounters& lower, const QuartetSharedCounters& upper);
QuartetCases shared_quartet_gg(const QuartetSharedCounters& first, const QuartetSharedCounters& second);

}  // namespace soda2013
