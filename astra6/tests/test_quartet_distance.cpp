// Sections 2–7 end to end. The naive oracle is linked only into this test.
#include "quartet_distance.hpp"
#include "quartet_topology_oracle.hpp"
#include "all_rooted_trees.hpp"

#include <cmath>

using namespace soda2013;
using namespace test_support;

std::array<bool, 9> shared_seen{}, disagreement_seen{};

void traversal_stats(const Tree& tree, const CountTraversalStats& stats) {
    CHECK(stats.queries == tree.size() - tree.leaf_count());
    CHECK(stats.hdt_builds == tree.leaf_count() + stats.heavy_contractions);
    CHECK(stats.maximum_components_per_active_leaf <= 32);
    CHECK(stats.paired_leaf_updates <= stats.hdt_input_components);
}

void compare(const Tree& first, const Tree& second) {
    QuartetDistanceStats stats;
    const auto actual = quartet_counts(first, second, &stats);
    const auto expected = oracle::quartet_counts(first, second);
    CHECK(actual.A == expected.A);
    CHECK(actual.B == expected.B);
    CHECK(actual.C == expected.C);
    CHECK(actual.D == expected.D);
    CHECK(actual.E == expected.E);
    CHECK(actual.total() == choose(first.leaf_count(), 4));
    CHECK(stats.agreeing_cases.values == quartet_oracle::shared_cases(first, second).values);
    CHECK(stats.disagreeing_cases.values == quartet_oracle::disagreement_cases(first, second).values);
    for (std::size_t i = 0; i < 9; ++i) {
        shared_seen[i] = shared_seen[i] || stats.agreeing_cases.values[i];
        disagreement_seen[i] = disagreement_seen[i] || stats.disagreeing_cases.values[i];
    }
    traversal_stats(first, stats.forward);
    traversal_stats(second, stats.reverse);
    for (long double p : {0.0L, 0.25L, 0.5L, 1.0L})
        CHECK(actual.parameterized(p) == expected.parameterized(p));
}

void check_work(const Tree& first, const Tree& second, Count expected) {
    QuartetDistanceStats stats;
    CHECK(quartet_counts(first, second, &stats).distance() == expected);
    std::size_t degree = 0;
    for (const Tree* tree : {&first, &second})
        for (NodeId v : tree->preorder()) degree = std::max(degree, tree->node(v).children.size());
    std::size_t work = 0;
    for (const auto* pass : {&stats.forward, &stats.reverse}) {
        work += pass->hdt_round_visits + pass->paired_leaf_updates + pass->compression.extraction_nodes + pass->compression.queue_pops +
                pass->updates.recomputed_components + pass->updates.color_rows + pass->updates.color_pairs;
        CHECK(pass->maximum_components_per_active_leaf <= 32);
    }
    const auto n = first.leaf_count();
    const double budget = 500.0 * static_cast<double>(degree * n) * (1 + std::log2(static_cast<double>(n)));
    CHECK(static_cast<double>(work) <= budget);
    std::cout << "n=" << n << " degree=" << degree << " quartet work=" << work
              << " contractions=" << stats.forward.heavy_contractions + stats.reverse.heavy_contractions << '\n';
}

int main() try {
    std::size_t pairs = 0;
    for (std::size_t n = 1; n <= 5; ++n) {
        std::vector<std::size_t> labels;
        for (std::size_t i = 0; i < n; ++i) labels.push_back(i);
        std::vector<Tree> trees;
        for (const auto& text : all_subtrees(labels)) trees.push_back(Tree::from_newick(text + ";"));
        for (const auto& first : trees)
            for (const auto& second : trees) { compare(first, second); ++pairs; }
    }
    std::mt19937_64 rng(7162013);
    for (std::size_t n = 1; n <= 12; ++n)
        for (int trial = 0; trial < 100; ++trial) {
            const auto first = random_tree(n, rng, trial % 2 ? 2 : 12);
            const auto second = random_tree(n, rng);
            compare(first, second);
            compare(second, first);
        }
    for (bool seen : shared_seen) CHECK(seen);
    for (bool seen : disagreement_seen) CHECK(seen);
    const auto a = Tree::from_newick("((a,b),(c,d));"), b = Tree::from_newick("((a,c),(b,d));");
    const auto unresolved = Tree::from_newick("(a,b,c,d);");
    CHECK(quartet_distance(a, b) == 1);
    CHECK(parameterized_quartet_distance(a, b, 0) == 1);
    CHECK(parameterized_quartet_distance(a, unresolved, 0.5L) == 0.5L);
    expect_exception<std::invalid_argument>([&] { (void)quartet_distance(a, Tree::from_newick("(a,b,c,x);")); });
    for (long double p : {-1.0L, 1.1L, std::numeric_limits<long double>::infinity(),
                          std::numeric_limits<long double>::quiet_NaN()})
        expect_exception<std::invalid_argument>([&] { (void)parameterized_quartet_distance(a, b, p); });
    for (std::size_t n : {100, 1000, 10000}) {
        const auto c = caterpillar(n), t = balanced(n);
        check_work(c, c, 0);
        check_work(t, t, 0);
        if (n <= 1000) {
            const auto s = star(n);
            check_work(s, s, 0);
            check_work(s, c, choose(n, 4));
            check_work(c, s, choose(n, 4));
        }
    }
    std::cout << "quartet distance: OK (" << pairs << " exhaustive pairs; 2,400 random directions; all nine A/B cases)\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
