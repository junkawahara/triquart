// naive_distance.hpp — O(n^3) / O(n^4) enumeration oracles.
//
// These are *test oracles only* (see CLAUDE.md): they enumerate every
// triplet / quartet, classify it into the five cases of Fig. 3 of
// soda13.pdf (A agree, B differ, C res/unres, D unres/res, E unres/unres)
// and are used to validate the HDT-based algorithms on small random trees.
#pragma once

#include <cstdint>

#include "count_type.hpp"
#include "tree.hpp"

namespace soda13 {

struct CaseCounts {
    count_t A = 0, B = 0, C = 0, D = 0, E = 0;
    count_t distance() const { return B + C + D; }
    // Parameterized distance of Bansal et al. (Sec. 2): B + p (C + D).
    // `p` is taken as long double so that a decimal p read with strtold is
    // reproduced exactly by an LDBL_DIG-digit print (trial 9, Finding 2 of
    // trial 8); a double argument converts implicitly.
    long double parameterized(long double p) const {
        return static_cast<long double>(B) + p * static_cast<long double>(C + D);
    }
};

CaseCounts naive_triplet_counts(const Tree& t1, const Tree& t2);
CaseCounts naive_quartet_counts(const Tree& t1, const Tree& t2);

// Rooted triplet topology of {a,b,c}: 0 = ab|c, 1 = ac|b, 2 = bc|a, 3 = star.
// Unrooted quartet topology of {a,b,c,d}: 0 = ab|cd, 1 = ac|bd, 2 = ad|bc, 3 = star.
struct NaiveTopology {
    explicit NaiveTopology(const Tree& t);
    int triplet(int a, int b, int c) const;   // leaf ids
    int quartet(int a, int b, int c, int d) const;
    int lca(int u, int v) const;              // node ids
    int dist(int u, int v) const;             // edge count
    std::vector<int> depth, parent, leaf_node;
};

}  // namespace soda13
