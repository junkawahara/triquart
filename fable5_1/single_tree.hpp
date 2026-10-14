// single_tree.hpp — Sec. 3 of soda13.pdf: counting resolved and unresolved
// triplets and quartets in a single tree in O(n).
//
// For a node v of degree d with children sizes n_1..n_d the paper's dynamic
// programme is
//
//     s_1 = n_1,  s_i = s_{i-1} + n_i
//     p_1 = 0,    p_i = p_{i-1} + n_i * s_{i-1}
//     t_1 = 0,    t_i = t_{i-1} + n_i * p_{i-1}
//     q_1 = 0,    q_i = q_{i-1} + n_i * t_{i-1}
//
// The number of unresolved triplets rooted at v is t_d, and the number of
// unresolved quartets rooted at v is q_d + t_d * (n - n^v) (Sec. 3, last
// paragraph).  Resolved counts are obtained by subtracting from C(n,3) and
// C(n,4).  These give the row/column sums of Fig. 3:
//
//     A + B + C = #resolved in T1,   A + B + D = #resolved in T2,
//     D + E     = #unresolved in T1, C + E     = #unresolved in T2.
#pragma once

#include <cstdint>

#include "count_type.hpp"
#include "tree.hpp"

namespace soda13 {

struct SingleTreeCounts {
    count_t n = 0;
    count_t resolved_triplets = 0;
    count_t unresolved_triplets = 0;
    count_t resolved_quartets = 0;
    count_t unresolved_quartets = 0;
};

SingleTreeCounts count_single_tree(const Tree& t);

count_t binom2(count_t n);
count_t binom3(count_t n);
count_t binom4(count_t n);

}  // namespace soda13
