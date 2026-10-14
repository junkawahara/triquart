// soda13.pdf, Figure 8 supplies the three missing cases by swapping the trees.
// Figure 3 then recovers C,D,E from A,B and the Section 3 resolved marginals.
// Production never enumerates triplets/quartets or calls the naive oracle.
#include "quartet_distance.hpp"
#include "count_traversal.hpp"

#include <cmath>

namespace soda2013 {

DistanceCounts quartet_counts(const Tree& first, const Tree& second, QuartetDistanceStats* stats) {
    validate_leaf_sets(first, second);
    const Count total = choose(first.leaf_count(), 4);
    QuartetDistanceStats local_stats;
    if (!stats) stats = &local_stats;
    *stats = {};
    const auto forward = detail::CountTraversal<ColoredQuartetHdt>(first, second, stats->forward).run();
    const auto reverse = detail::CountTraversal<ColoredQuartetHdt>(second, first, stats->reverse).run();
    stats->agreeing_cases = complete_quartet_cases(forward.agreeing, reverse.agreeing);
    stats->disagreeing_cases = complete_quartet_cases(forward.disagreeing, reverse.disagreeing);
    return from_resolved(total, count_single_tree(first).quartets.resolved,
                         count_single_tree(second).quartets.resolved,
                         stats->agreeing_cases.total(), stats->disagreeing_cases.total());
}

Count quartet_distance(const Tree& first, const Tree& second) { return quartet_counts(first, second).distance(); }

long double parameterized_quartet_distance(const Tree& first, const Tree& second, long double p) {
    if (!std::isfinite(p) || p < 0 || p > 1)
        throw std::invalid_argument("distance parameter must be finite and in [0,1]");
    return quartet_counts(first, second).parameterized(p);
}

}  // namespace soda2013
