// Sections 5 and 7: batch-update path union and compression must preserve
// every Figure 13/16 case under arbitrary recolorings of the retained leaves.
#include "quartet_topology_oracle.hpp"
#include "quartet_hdt.hpp"

using namespace soda2013;
using namespace test_support;

ResolvedQuartets colored_oracle(const Tree& tree, const std::vector<Color>& colors) {
    ResolvedQuartets result;
    const auto& leaves = tree.leaves();
    for (std::size_t i = 0; i < leaves.size(); ++i)
        for (std::size_t j = i + 1; j < leaves.size(); ++j)
            for (std::size_t k = j + 1; k < leaves.size(); ++k)
                for (std::size_t l = k + 1; l < leaves.size(); ++l) {
                    const std::array<NodeId, 4> selected{leaves[i], leaves[j], leaves[k], leaves[l]};
                    const auto first = quartet_oracle::coloring(colors, selected);
                    const auto second = quartet_oracle::topology(tree, selected);
                    if (!first.split || !second.split || first.type < second.type) continue;
                    auto& cases = first.split == second.split ? result.agreeing : result.disagreeing;
                    ++cases.at(first.type, second.type);
                }
    return result;
}

void equal(const ResolvedQuartets& first, const ResolvedQuartets& second) {
    CHECK(first.agreeing.values == second.agreeing.values);
    CHECK(first.disagreeing.values == second.disagreeing.values);
}

void check_updates(const Tree& tree, std::mt19937_64& rng) {
    ColoredQuartetHdt hdt(tree);
    std::vector<Color> colors(tree.size(), 1);
    expect_exception<std::invalid_argument>([&] { hdt.recolor({{tree.leaves()[0], 2}, {no_node, 3}}); });
    CHECK(hdt.color(tree.leaves()[0]) == 1);
    for (int batch = 0; batch < 10; ++batch) {
        std::vector<Recoloring> changes;
        std::set<NodeId> paths;
        for (NodeId leaf : tree.leaves()) {
            if (rng() % 3 == 0) continue;
            const Color color = rng() % 5 * (batch % 2 ? 100000003 : 1);
            changes.push_back({leaf, color});
            if (colors[leaf] != color) {
                NodeId uid = hdt.hdt().leaf_for_vertex(leaf);
                while (uid != no_node) { paths.insert(uid); uid = hdt.hdt().node(uid).parent; }
                colors[leaf] = color;
            }
        }
        const auto before = hdt.stats().recomputed_components;
        hdt.recolor(changes);
        CHECK(hdt.stats().recomputed_components - before == paths.size());
        equal(hdt.query(), colored_oracle(tree, colors));
    }
    const auto before = hdt.stats().recomputed_components;
    hdt.recolor({});
    CHECK(hdt.stats().recomputed_components == before);
    const NodeId leaf = tree.leaves()[0];
    hdt.recolor({{leaf, 19}, {leaf, 0}, {leaf, 7}});
    CHECK(hdt.color(leaf) == 7);
    colors[leaf] = 7;
    equal(hdt.query(), colored_oracle(tree, colors));
}

void check_compressed(const Tree& tree, const std::vector<NodeId>& kept,
                      ComponentTree compressed, std::mt19937_64& rng) {
    CHECK(compressed.size() <= 8 * kept.size());
    ColoredQuartetHdt full(tree, 0), small(std::move(compressed));
    for (int trial = 0; trial < 5; ++trial) {
        std::vector<Recoloring> changes;
        for (NodeId leaf : kept) changes.push_back({leaf, rng() % 5});
        full.recolor(changes);
        small.recolor(changes);
        equal(full.query(), small.query());
        CHECK(full.counters(full.hdt().root()).shared[SharedScalar::zero] ==
              small.counters(small.hdt().root()).shared[SharedScalar::zero]);
    }
    if (kept.size() > 1) {
        const std::vector<NodeId> nested(kept.begin(), kept.begin() + static_cast<std::ptrdiff_t>(kept.size() / 2));
        ColoredQuartetHdt nested_hdt(small.extract(nested)), reference(tree, 0);
        std::vector<Recoloring> changes;
        for (NodeId leaf : nested) changes.push_back({leaf, 1 + rng() % 4});
        nested_hdt.recolor(changes);
        reference.recolor(changes);
        equal(nested_hdt.query(), reference.query());
    }
}

int main() try {
    std::mt19937_64 rng(5172013);
    std::size_t subsets = 0;
    for (std::size_t n = 1; n <= 12; ++n) {
        for (int trial = 0; trial < 25; ++trial) {
            const auto tree = trial == 0 ? star(n) : trial == 1 ? caterpillar(n) : random_tree(n, rng);
            check_updates(tree, rng);
            ColoredQuartetHdt full(tree, 4);
            const auto original = full.query();
            const std::size_t count = n <= 6 ? (std::size_t(1) << n) - 1 : 12;
            for (std::size_t subset = 0; subset < count; ++subset) {
                const std::size_t mask = n <= 6 ? subset + 1 : rng();
                std::vector<NodeId> kept;
                for (std::size_t i = 0; i < n; ++i)
                    if (mask & (std::size_t(1) << i)) kept.push_back(tree.leaves()[i]);
                if (kept.empty()) kept.push_back(tree.leaves()[0]);
                check_compressed(tree, kept, full.extract(kept), rng);
                check_compressed(tree, kept, full.contract(kept), rng);
                equal(full.query(), original);
                ++subsets;
            }
        }
    }
    std::cout << "quartet HDT: OK (3,000 update batches; " << subsets
              << " subsets, both compression methods, nested extracts, every case)\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
