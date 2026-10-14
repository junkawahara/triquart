// End-to-end checks of Sections 2–5 / Figure 10, with the independent naive
// oracle linked ONLY into this executable. Exhaustive rooted multifurcating
// trees include every assignment of up to five labelled leaves to topologies.
#include "naive_distance.hpp"
#include "test_support.hpp"
#include "triplet_distance.hpp"

#include <cmath>
#include "all_rooted_trees.hpp"

using namespace soda2013;
using namespace test_support;


void compare(const Tree& first, const Tree& second) {
    TripletDistanceStats stats;
    const auto actual = triplet_counts(first, second, &stats);
    const auto expected = oracle::triplet_counts(first, second);
    CHECK(actual.A == expected.A);
    CHECK(actual.B == expected.B);
    CHECK(actual.C == expected.C);
    CHECK(actual.D == expected.D);
    CHECK(actual.E == expected.E);
    CHECK(actual.total() == choose(first.leaf_count(), 3));
    CHECK(stats.queries == first.size() - first.leaf_count());
    CHECK(stats.hdt_builds == first.leaf_count() + stats.heavy_contractions);
    CHECK(stats.maximum_components_per_active_leaf <= 32);
    CHECK(stats.paired_leaf_updates <= stats.hdt_input_components);
    for (long double p : {0.0L, 0.25L, 0.5L, 1.0L})
        CHECK(actual.parameterized(p) == expected.parameterized(p));
}

void check_work(const Tree& first, const Tree& second, Count expected) {
    TripletDistanceStats stats;
    CHECK(triplet_counts(first, second, &stats).distance() == expected);
    const auto n = first.leaf_count();
    const auto work = stats.hdt_round_visits + stats.paired_leaf_updates + stats.compression.extraction_nodes +
                      stats.compression.queue_pops + stats.updates.recomputed_components +
                      stats.updates.color_rows;
    const double budget = 500.0 * static_cast<double>(n) * (1 + std::log2(static_cast<double>(n)));
    CHECK(static_cast<double>(work) <= budget);
    CHECK(stats.maximum_components_per_active_leaf <= 32);
    CHECK(stats.hdt_builds == n + stats.heavy_contractions);
    std::cout << "n=" << n << " structural/counter work=" << work
              << " heavy contractions=" << stats.heavy_contractions << '\n';
}

int main() try {
    std::size_t pairs = 0;
    const std::vector<std::size_t> expected_sizes{0, 1, 1, 4, 26, 236};
    for (std::size_t n = 1; n <= 5; ++n) {
        std::vector<std::size_t> labels;
        for (std::size_t i = 0; i < n; ++i) labels.push_back(i);
        const auto newicks = all_subtrees(labels);
        CHECK(newicks.size() == expected_sizes[n]);
        std::vector<Tree> trees;
        for (const auto& newick : newicks) trees.push_back(Tree::from_newick(newick + ";"));
        for (const auto& first : trees)
            for (const auto& second : trees) { compare(first, second); ++pairs; }
    }
    std::mt19937_64 rng(5102013);
    for (std::size_t n = 1; n <= 12; ++n)
        for (int trial = 0; trial < 100; ++trial) {
            const auto first = random_tree(n, rng, trial % 2 ? 2 : 12);
            const auto second = random_tree(n, rng);
            compare(first, second);
            compare(second, first);
        }
    const auto a = Tree::from_newick("((a,b),c);"), b = Tree::from_newick("(a,b,c);");
    CHECK(triplet_distance(a, b) == 1);
    CHECK(parameterized_triplet_distance(a, b, 0.5L) == 0.5L);
    expect_exception<std::invalid_argument>([&] { (void)parameterized_triplet_distance(a, b, -1); });
    expect_exception<std::invalid_argument>([&] { (void)triplet_distance(a, Tree::from_newick("(x,y,z);")); });
    for (std::size_t n : {100, 1000, 10000}) {
        const auto s = star(n), c = caterpillar(n), t = balanced(n);
        check_work(s, s, 0);
        check_work(c, c, 0);
        check_work(t, t, 0);
        check_work(s, c, choose(n, 3));
        check_work(c, s, choose(n, 3));
    }
    std::cout << "triplet distance: OK (" << pairs << " exhaustive pairs, 2,400 random directions)\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
