// Section 3 of soda13.pdf: s,p,t,q count leaves chosen from distinct children.
// Section 2 gives the identities reconstructing the five cells of Figure 3.
#include "distance_statistics.hpp"

#include <cmath>

namespace soda2013 {

namespace {
SingleTreeCounts count_topologies(const Tree& tree, bool quartets) {
    SingleTreeCounts result;
    const Count n = tree.leaf_count();
    for (NodeId v : tree.preorder()) {
        Count s = 0, p = 0, t = 0, q = 0;
        for (NodeId child : tree.node(v).children) {
            const Count size = tree.subtree_leaves(child);
            if (quartets) q = add(q, multiply(size, t));
            t = add(t, multiply(size, p));
            p = add(p, multiply(size, s));
            s = add(s, size);
        }
        result.triplets.unresolved = add(result.triplets.unresolved, t);
        if (quartets)
            result.quartets.unresolved = add(result.quartets.unresolved,
                                             add(q, multiply(t, n - tree.subtree_leaves(v))));
    }
    result.triplets.resolved = subtract(choose(n, 3), result.triplets.unresolved);
    if (quartets) result.quartets.resolved = subtract(choose(n, 4), result.quartets.unresolved);
    return result;
}
}  // namespace

SingleTreeCounts count_single_tree(const Tree& tree) { return count_topologies(tree, true); }
TopologyCounts count_triplet_topologies(const Tree& tree) { return count_topologies(tree, false).triplets; }

Count DistanceCounts::total() const { return add(add(add(add(A, B), C), D), E); }
Count DistanceCounts::distance() const { return add(add(B, C), D); }
long double DistanceCounts::parameterized(long double p) const {
    if (!std::isfinite(p) || p < 0 || p > 1)
        throw std::invalid_argument("distance parameter must be finite and in [0,1]");
    return static_cast<long double>(B) + p * static_cast<long double>(add(C, D));
}

DistanceCounts from_resolved(Count total, Count r1, Count r2, Count A, Count B) {
    const Count both = add(A, B);
    const Count C = subtract(r1, both);
    const Count D = subtract(r2, both);
    const Count E = subtract(total, add(add(both, C), D));
    return {A, B, C, D, E};
}

DistanceCounts from_agreeing(Count total, Count r1, Count r2, Count A, Count E) {
    const Count C = subtract(subtract(total, r2), E);
    const Count D = subtract(subtract(total, r1), E);
    const Count B = subtract(subtract(r1, A), C);
    if (add(add(A, B), D) != r2) throw std::logic_error("inconsistent distance marginals");
    return {A, B, C, D, E};
}

}  // namespace soda2013
