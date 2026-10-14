#pragma once

// soda13.pdf, Sections 2–7: full Figure 7 Count with extract/contract, all
// Figures 11–16 decorations, and both Figure 8 tree orientations.
#include "distance_statistics.hpp"
#include "quartet_hdt.hpp"
#include "traversal_stats.hpp"

namespace soda2013 {

struct QuartetDistanceStats {
    CountTraversalStats forward;
    CountTraversalStats reverse;
    QuartetCases agreeing_cases;
    QuartetCases disagreeing_cases;
};

// O(d n log n) structural/counter work, where d is the largest input degree.
// Labels take O(n log n) comparisons; paired-leaf access during Count is O(1).
// stats, when supplied, is reset and contains all nine cases and both passes.
DistanceCounts quartet_counts(const Tree& first, const Tree& second, QuartetDistanceStats* stats = nullptr);
Count quartet_distance(const Tree& first, const Tree& second);
long double parameterized_quartet_distance(const Tree& first, const Tree& second, long double p);

}  // namespace soda2013
