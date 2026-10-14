// distance.hpp — public API: triplet and quartet distances between two
// rooted / unrooted trees of arbitrary degree (soda13.pdf).
//
//   triplet distance   O(n log n)      (Sec. 5–7, Fig. 10)
//   quartet distance   O(d n log n)    (Sec. 5–7, Fig. 13 and Fig. 16)
//   parameterized variants (Bansal et al., Sec. 2):  B + p (C + D)
//
// All routines take trees sharing the same non-empty leaf-label set and
// throw std::invalid_argument otherwise.  The returned `CaseCounts` hold
// the five cells of Fig. 3: A (resolved, agree), B (resolved, differ),
// C (resolved in T1 only), D (resolved in T2 only), E (unresolved in both).
#pragma once

#include "naive_distance.hpp"
#include "tree.hpp"

namespace soda13 {

CaseCounts triplet_case_counts(const Tree& t1, const Tree& t2);
CaseCounts quartet_case_counts(const Tree& t1, const Tree& t2);

// Totals are 128-bit (count_t): C(n,4) exceeds 2^63 for n > ~123 000.
count_t triplet_distance(const Tree& t1, const Tree& t2);
count_t quartet_distance(const Tree& t1, const Tree& t2);
double parameterized_triplet_distance(const Tree& t1, const Tree& t2, double p);
double parameterized_quartet_distance(const Tree& t1, const Tree& t2, double p);

}  // namespace soda13
