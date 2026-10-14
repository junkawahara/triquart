// soda13.pdf, Sections 2–5: shared A/E from Count and Section 3 marginals.
#include "triplet_distance.hpp"
#include "count_traversal.hpp"

#include <cmath>

namespace soda2013 {

DistanceCounts triplet_counts(const Tree& first, const Tree& second, TripletDistanceStats* stats) {
    validate_leaf_sets(first, second);
    const Count total = choose(first.leaf_count(), 3);
    TripletDistanceStats local_stats;
    if (!stats) stats = &local_stats;
    *stats = {};
    const auto shared = detail::CountTraversal<ColoredTripletHdt>(first, second, *stats).run();
    return from_agreeing(total, count_triplet_topologies(first).resolved,
                        count_triplet_topologies(second).resolved, shared.resolved, shared.unresolved);
}

Count triplet_distance(const Tree& first, const Tree& second) { return triplet_counts(first, second).distance(); }

long double parameterized_triplet_distance(const Tree& first, const Tree& second, long double p) {
    if (!std::isfinite(p) || p < 0 || p > 1)
        throw std::invalid_argument("distance parameter must be finite and in [0,1]");
    return triplet_counts(first, second).parameterized(p);
}

}  // namespace soda2013
