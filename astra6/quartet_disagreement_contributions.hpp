#pragma once

// soda13.pdf, Figure 16: all six cases for one orientation. Figure 8 completes
// the nine-case matrix using complete_quartet_cases from the reverse traversal.
#include "quartet_contributions.hpp"
#include "quartet_disagreement_counters.hpp"

namespace soda2013 {

// Contributions crossing this composition, without accumulated child totals.
// CC operands are the LOWER C1 followed by the UPPER C2 (Figure 6).
QuartetCases disagreement_quartet_cc(const QuartetCounters& lower, const QuartetCounters& upper);
QuartetCases disagreement_quartet_gg(const QuartetCounters& first, const QuartetCounters& second);

}  // namespace soda2013
