// single_tree.cpp — see single_tree.hpp (Sec. 3 of soda13.pdf).
#include "single_tree.hpp"

namespace soda13 {

count_t binom2(count_t n) { return n < 2 ? 0 : n * (n - 1) / 2; }
count_t binom3(count_t n) { return n < 3 ? 0 : n * (n - 1) / 2 * (n - 2) / 3; }
count_t binom4(count_t n) { return n < 4 ? 0 : n * (n - 1) / 2 * (n - 2) / 3 * (n - 3) / 4; }

SingleTreeCounts count_single_tree(const Tree& t) {
    SingleTreeCounts r;
    const std::vector<int> cnt = t.subtree_leaf_counts();
    const count_t n = cnt[t.root];
    r.n = n;
    count_t unres_trip = 0, unres_quart = 0;
    for (int v = 0; v < t.num_nodes(); ++v) {
        const auto& ch = t.nodes[v].children;
        if (ch.empty()) continue;
        // Sec. 3 recurrences over the children in order.
        count_t s = 0, p = 0, tt = 0, q = 0;
        for (int c : ch) {
            const count_t ni = cnt[c];
            q += ni * tt;
            tt += ni * p;
            p += ni * s;
            s += ni;
        }
        unres_trip += tt;
        unres_quart += q + tt * (n - cnt[v]);
    }
    r.unresolved_triplets = unres_trip;
    r.resolved_triplets = binom3(n) - unres_trip;
    r.unresolved_quartets = unres_quart;
    r.resolved_quartets = binom4(n) - unres_quart;
    return r;
}

}  // namespace soda13
