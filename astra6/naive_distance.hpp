#pragma once

// Test oracle ONLY. Direct enumeration from soda13.pdf, Section 2 / Figures 1–3.
// This module is deliberately separate from the production distance routines.
#include "distance_statistics.hpp"

namespace soda2013::oracle {

// Topology zero is unresolved; 1,2,3 identify ab|c, ac|b, bc|a (triplets)
// or ab|cd, ac|bd, ad|bc (quartets), for the arguments in the given order.
int triplet_topology(const Tree& tree, NodeId a, NodeId b, NodeId c);
int quartet_topology(const Tree& tree, NodeId a, NodeId b, NodeId c, NodeId d);
DistanceCounts triplet_counts(const Tree& first, const Tree& second);
DistanceCounts quartet_counts(const Tree& first, const Tree& second);

}  // namespace soda2013::oracle
