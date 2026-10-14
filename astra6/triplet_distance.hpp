#pragma once

// soda13.pdf, Sections 2–5 and 7.1–7.2: Count with extract/contract, using the
// Figure 10 triplet decoration. No calls to naive_distance.cpp are permitted.
#include "distance_statistics.hpp"
#include "triplet_counters.hpp"
#include "traversal_stats.hpp"

namespace soda2013 {

using TripletDistanceStats = CountTraversalStats;

// Section 5: O(n log n) time. Labels are matched in O(n log n) comparisons;
// the traversal uses direct paired-leaf indices, without hash assumptions.
// If supplied, stats are reset and then filled for this invocation.
DistanceCounts triplet_counts(const Tree& first, const Tree& second,
                              TripletDistanceStats* stats = nullptr);
Count triplet_distance(const Tree& first, const Tree& second);
long double parameterized_triplet_distance(const Tree& first, const Tree& second, long double p);

}  // namespace soda2013
