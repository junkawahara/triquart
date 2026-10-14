#pragma once

// soda13.pdf, Section 2 / Figure 3 and Section 3.
#include "tree.hpp"

namespace soda2013 {

struct TopologyCounts {
    Count resolved = 0;
    Count unresolved = 0;
};

struct SingleTreeCounts {
    TopologyCounts triplets;
    TopologyCounts quartets;
};

// O(number of vertices) dynamic program from Section 3.
SingleTreeCounts count_single_tree(const Tree& tree);
// Triplet-only counting does not evaluate quartet counters that could overflow
// on an input whose triplet counts are still representable.
TopologyCounts count_triplet_topologies(const Tree& tree);

struct DistanceCounts {
    Count A = 0;  // Resolved in both, same topology.
    Count B = 0;  // Resolved in both, different topologies.
    Count C = 0;  // Resolved only in the first tree.
    Count D = 0;  // Resolved only in the second tree.
    Count E = 0;  // Unresolved in both.

    Count total() const;
    Count distance() const;
    long double parameterized(long double p) const;
};

// Complete Figure 3 from its marginals and the joint counts, Section 2.
DistanceCounts from_agreeing(Count total, Count resolved_first,
                             Count resolved_second, Count A, Count E);
DistanceCounts from_resolved(Count total, Count resolved_first,
                             Count resolved_second, Count A, Count B);

}  // namespace soda2013
