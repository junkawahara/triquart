// Independent checks for Sections 2–3 of soda13.pdf. Later HDT tests use the
// oracle here, so first check its topology classification against hand examples.
#include "distance_statistics.hpp"
#include "naive_distance.hpp"
#include "test_support.hpp"

#include <limits>

using namespace soda2013;
using namespace test_support;

int main() try {
    CHECK(choose(2, 3) == 0);
    CHECK(choose(12, 4) == 495);
    CHECK(choose(67, 33) == 14226520737620288370ULL);
    expect_exception<std::overflow_error>([] { (void)choose(68, 34); });
    expect_exception<std::overflow_error>([] { (void)add(~Count(0), 1); });
    expect_exception<std::overflow_error>([] { (void)multiply(~Count(0), 2); });
    CHECK(multiply(~Count(0), 0) == 0);
    expect_exception<std::logic_error>([] { (void)subtract(0, 1); });

    for (const std::string text : {"", ";", "();", "(a);", "(a,b)", "(a,a);",
                                    "(a,,b);", "(,a);", "(a,b));", "a;b;", "a:;",
                                    "a:nan;", "a:1x;", "a:1e9999;", "'a;", "a[;"})
        expect_exception<std::invalid_argument>([&] { (void)Tree::from_newick(text); });
    expect_exception<std::invalid_argument>([] { (void)Tree({}, 0); });
    expect_exception<std::invalid_argument>([] { (void)Tree({{"a", {}}}, 1); });
    expect_exception<std::invalid_argument>([] { (void)Tree({{"", {1, 1}}, {"a", {}}}, 0); });
    expect_exception<std::invalid_argument>([] { (void)Tree({{"", {0, 1}}, {"a", {}}}, 0); });
    expect_exception<std::invalid_argument>([] { (void)Tree({{"a", {}}, {"b", {}}}, 0); });
    const auto quoted = Tree::from_newick(" [hi] ('a b':1.2e-3,'x''y':2,(c,d)[[note]])root:0; ");
    CHECK(quoted.leaf_count() == 4);
    CHECK(quoted.node(quoted.leaf("x'y")).label == "x'y");
    CHECK(quoted.subtree_leaves(quoted.root()) == 4);

    const auto ab_c = Tree::from_newick("((a,b),c);");
    const auto ac_b = Tree::from_newick("((a,c),b);");
    const auto abc = Tree::from_newick("(a,b,c);");
    CHECK(oracle::triplet_counts(ab_c, ab_c).A == 1);
    CHECK(oracle::triplet_counts(ab_c, ac_b).B == 1);
    CHECK(oracle::triplet_counts(ab_c, abc).C == 1);
    CHECK(oracle::triplet_counts(abc, ab_c).D == 1);
    CHECK(oracle::triplet_counts(abc, abc).E == 1);
    CHECK(oracle::triplet_topology(ab_c, ab_c.leaf("c"), ab_c.leaf("a"), ab_c.leaf("b")) == 3);
    expect_exception<std::invalid_argument>([&] {
        (void)oracle::triplet_counts(abc, Tree::from_newick("(a,b,z);"));
    });
    expect_exception<std::invalid_argument>([&] { validate_leaf_sets(abc, quoted); });

    const auto ab_cd = Tree::from_newick("((a,b),(c,d));");
    const auto ac_bd = Tree::from_newick("((a,c),(b,d));");
    const auto ad_bc = Tree::from_newick("((a,d),(b,c));");
    const auto rerooted = Tree::from_newick("(a,b,(c,d));");
    const auto abcd = Tree::from_newick("(a,b,c,d);");
    CHECK(oracle::quartet_counts(ab_cd, rerooted).A == 1);
    CHECK(oracle::quartet_counts(ab_cd, ac_bd).B == 1);
    CHECK(oracle::quartet_counts(ac_bd, ad_bc).B == 1);
    CHECK(oracle::quartet_counts(ab_cd, abcd).C == 1);
    CHECK(oracle::quartet_counts(abcd, ab_cd).D == 1);
    CHECK(oracle::quartet_counts(abcd, abcd).E == 1);
    CHECK(oracle::quartet_counts(ad_bc, ad_bc).A == 1);

    const DistanceCounts categories{1, 2, 3, 4, 5};
    CHECK(categories.total() == 15);
    CHECK(categories.distance() == 9);
    CHECK(categories.parameterized(0) == 2);
    CHECK(categories.parameterized(0.5L) == 5.5L);
    CHECK(categories.parameterized(1) == 9);
    for (long double p : {-0.1L, 1.1L, std::numeric_limits<long double>::quiet_NaN(),
                           std::numeric_limits<long double>::infinity()})
        expect_exception<std::invalid_argument>([&] { (void)categories.parameterized(p); });
    CHECK(from_agreeing(15, 6, 7, 1, 5).B == 2);
    CHECK(from_resolved(15, 6, 7, 1, 2).E == 5);

    std::mt19937_64 rng(2013);
    for (std::size_t n = 1; n <= 12; ++n) {
        for (int trial = 0; trial < 100; ++trial) {
            const auto tree = random_tree(n, rng);
            const auto dp = count_single_tree(tree);
            const auto triplets = oracle::triplet_counts(tree, tree);
            const auto quartets = oracle::quartet_counts(tree, tree);
            CHECK(dp.triplets.resolved == triplets.A);
            CHECK(dp.triplets.unresolved == triplets.E);
            CHECK(dp.quartets.resolved == quartets.A);
            CHECK(dp.quartets.unresolved == quartets.E);
            CHECK(triplets.total() == choose(n, 3));
            CHECK(quartets.total() == choose(n, 4));
        }
        CHECK(count_single_tree(star(n)).triplets.unresolved == choose(n, 3));
        CHECK(count_single_tree(star(n)).quartets.unresolved == choose(n, 4));
        CHECK(count_single_tree(caterpillar(n)).triplets.resolved == choose(n, 3));
        CHECK(count_single_tree(caterpillar(n)).quartets.resolved == choose(n, 4));
    }
    const auto deep = caterpillar(20000);
    CHECK(count_single_tree(deep).triplets.resolved == choose(20000, 3));
    CHECK(count_single_tree(deep).quartets.resolved == choose(20000, 4));
    std::string newick(9999, '(');
    newick += "0";
    for (int i = 1; i < 10000; ++i) newick += "," + std::to_string(i) + ")";
    newick += ";";
    CHECK(Tree::from_newick(newick).leaf_count() == 10000);
    // Quartet counts overflow here, but triplet-only counting must still work.
    CHECK(count_triplet_topologies(star(200000)).unresolved == choose(200000, 3));
    std::cout << "foundations: OK (1,200 random trees, n <= 12; deep trees; input validation)\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
