// tests/test_distance.cpp — HDT algorithms vs. the naive oracle on random
// small trees (CLAUDE.md: "the naive oracle is the source of truth").
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../count.hpp"
#include "../distance.hpp"
#include "../naive_distance.hpp"
#include "test_util.hpp"

using namespace soda13;

static bool same(const CaseCounts& a, const CaseCounts& b) {
    return a.A == b.A && a.B == b.B && a.C == b.C && a.D == b.D && a.E == b.E;
}

int main(int argc, char** argv) {
    int iters = argc > 1 ? std::atoi(argv[1]) : 400;
    int nmax = argc > 2 ? std::atoi(argv[2]) : 14;
    bool do_trip = true, do_quart = true;
    if (argc > 3) {
        do_trip = std::strstr(argv[3], "t") != nullptr;
        do_quart = std::strstr(argv[3], "q") != nullptr;
        if (std::strstr(argv[3], "n")) count_no_compress() = true;  // "n": no contract/extract
    }
    std::mt19937 rng(2024);
    int fails = 0;
    for (int iter = 0; iter < iters; ++iter) {
        int n = 3 + iter % (nmax - 2);
        Tree t1 = test::random_tree(n, rng, (iter % 5) * 0.25);
        Tree t2 = test::random_tree(n, rng, ((iter / 5) % 4) * 0.33);
        if (iter % 7 == 0) t2 = t1;  // identical topologies
        test::shuffle_labels(t2, rng);
        assign_leaf_ids(t1, t2);
        if (do_trip) {
            CaseCounts got = triplet_case_counts(t1, t2), exp = naive_triplet_counts(t1, t2);
            if (!same(got, exp)) {
                ++fails;
                std::printf("TRIPLET FAIL iter %d n=%d\n  T1=%s\n  T2=%s\n  got A=%s B=%s C=%s D=%s E=%s\n  exp A=%s B=%s C=%s D=%s E=%s\n",
                            iter, n, t1.to_newick().c_str(), t2.to_newick().c_str(), count_to_string(got.A).c_str(), count_to_string(got.B).c_str(),
                            count_to_string(got.C).c_str(), count_to_string(got.D).c_str(), count_to_string(got.E).c_str(), count_to_string(exp.A).c_str(), count_to_string(exp.B).c_str(),
                            count_to_string(exp.C).c_str(), count_to_string(exp.D).c_str(), count_to_string(exp.E).c_str());
            }
        }
        if (do_quart) {
            CaseCounts got = quartet_case_counts(t1, t2), exp = naive_quartet_counts(t1, t2);
            if (!same(got, exp)) {
                ++fails;
                std::printf("QUARTET FAIL iter %d n=%d\n  T1=%s\n  T2=%s\n  got A=%s B=%s C=%s D=%s E=%s\n  exp A=%s B=%s C=%s D=%s E=%s\n",
                            iter, n, t1.to_newick().c_str(), t2.to_newick().c_str(), count_to_string(got.A).c_str(), count_to_string(got.B).c_str(),
                            count_to_string(got.C).c_str(), count_to_string(got.D).c_str(), count_to_string(got.E).c_str(), count_to_string(exp.A).c_str(), count_to_string(exp.B).c_str(),
                            count_to_string(exp.C).c_str(), count_to_string(exp.D).c_str(), count_to_string(exp.E).c_str());
            }
        }
        if (fails >= 5) break;
    }
    // Second phase: very high-degree T1 (every new leaf becomes an extra
    // child), so that the pair-table (★) rows see many colours per node.
    for (int iter = 0; iter < iters / 2 && fails < 5; ++iter) {
        int n = 5 + iter % (nmax - 4);
        Tree t1 = test::random_tree(n, rng, 0.0);
        Tree t2 = test::random_tree(n, rng, (iter % 3) * 0.35);
        test::shuffle_labels(t2, rng);
        assign_leaf_ids(t1, t2);
        if (do_trip) {
            CaseCounts got = triplet_case_counts(t1, t2), exp = naive_triplet_counts(t1, t2);
            if (!same(got, exp)) {
                ++fails;
                std::printf("TRIPLET FAIL (high degree) iter %d n=%d\n  T1=%s\n  T2=%s\n", iter, n, t1.to_newick().c_str(),
                            t2.to_newick().c_str());
            }
        }
        if (do_quart) {
            CaseCounts got = quartet_case_counts(t1, t2), exp = naive_quartet_counts(t1, t2);
            if (!same(got, exp)) {
                ++fails;
                std::printf("QUARTET FAIL (high degree) iter %d n=%d\n  T1=%s\n  T2=%s\n", iter, n, t1.to_newick().c_str(),
                            t2.to_newick().c_str());
            }
        }
    }
    std::printf("test_distance: %s (%d failures)\n", fails ? "FAILED" : "ok", fails);
    return fails ? 1 : 0;
}
