// tests/test_single_tree.cpp — Sec. 3 counts vs. the naive oracle.
#include <cstdio>
#include <cstdlib>

#include "../naive_distance.hpp"
#include "../single_tree.hpp"
#include "test_util.hpp"

using namespace soda13;

int main() {
    std::mt19937 rng(12345);
    int fails = 0;
    for (int iter = 0; iter < 300; ++iter) {
        int n = 3 + iter % 14;
        Tree t1 = test::random_tree(n, rng, (iter % 5) * 0.25);
        Tree t2 = test::random_tree(n, rng, (iter % 3) * 0.5);
        assign_leaf_ids(t1, t2);
        t1.contract_unary();
        t2.contract_unary();
        SingleTreeCounts c1 = count_single_tree(t1), c2 = count_single_tree(t2);
        CaseCounts trip = naive_triplet_counts(t1, t2), quart = naive_quartet_counts(t1, t2);
        bool ok = true;
        ok &= c1.resolved_triplets == trip.A + trip.B + trip.C;
        ok &= c2.resolved_triplets == trip.A + trip.B + trip.D;
        ok &= c1.unresolved_triplets == trip.D + trip.E;
        ok &= c2.unresolved_triplets == trip.C + trip.E;
        ok &= c1.resolved_quartets == quart.A + quart.B + quart.C;
        ok &= c2.resolved_quartets == quart.A + quart.B + quart.D;
        ok &= c1.unresolved_quartets == quart.D + quart.E;
        ok &= c2.unresolved_quartets == quart.C + quart.E;
        if (!ok) {
            ++fails;
            std::printf("FAIL iter %d n=%d\n  T1=%s\n  T2=%s\n", iter, n, t1.to_newick().c_str(), t2.to_newick().c_str());
        }
    }
    std::printf("test_single_tree: %s (%d failures)\n", fails ? "FAILED" : "ok", fails);
    return fails ? 1 : 0;
}
