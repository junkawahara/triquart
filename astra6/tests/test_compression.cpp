// Section 5: compression must preserve queries under ANY coloring of the kept
// leaves, reset discarded leaves to zero, preserve their multiplicity, and use
// linear contraction work / only the selected HDT path union for extraction.
#include "test_support.hpp"
#include "triplet_counters.hpp"

using namespace soda2013;
using namespace test_support;

void compare(const Tree& original, const std::vector<NodeId>& kept,
             ComponentTree compressed, std::mt19937_64& rng) {
    CHECK(compressed.size() <= 8 * kept.size());
    std::size_t active = 0;
    for (const auto& atom : compressed.atoms()) if (atom.active_leaf != no_node) ++active;
    CHECK(active == kept.size());
    ColoredTripletHdt full(original, 0), small(std::move(compressed));
    CHECK(small.hdt().node(small.hdt().root()).leaf_mass == original.leaf_count());
    for (int coloring = 0; coloring < 5; ++coloring) {
        std::vector<Recoloring> changes;
        for (NodeId leaf : kept) changes.push_back({leaf, rng() % (kept.size() + 1)});
        full.recolor(changes);
        small.recolor(changes);
        CHECK(full.query().resolved == small.query().resolved);
        CHECK(full.query().unresolved == small.query().unresolved);
        CHECK(full.counters(full.hdt().root()).zero == small.counters(small.hdt().root()).zero);
        CHECK(full.counters(full.hdt().root()).positive == small.counters(small.hdt().root()).positive);
    }
    if (kept.size() > 1) {
        const std::vector<NodeId> nested(kept.begin(), kept.begin() + static_cast<std::ptrdiff_t>(kept.size() / 2));
        auto extracted = small.extract(nested);
        CHECK(extracted.size() <= 8 * nested.size());
        ColoredTripletHdt nested_hdt(std::move(extracted)), reference(original, 0);
        std::vector<Recoloring> changes;
        for (NodeId leaf : nested) changes.push_back({leaf, 1 + rng() % 4});
        reference.recolor(changes);
        nested_hdt.recolor(changes);
        CHECK(reference.query().resolved == nested_hdt.query().resolved);
        CHECK(reference.query().unresolved == nested_hdt.query().unresolved);
        CHECK(reference.counters(reference.hdt().root()).zero == nested_hdt.counters(nested_hdt.hdt().root()).zero);
    }
}

void check_subsets(const Tree& tree, std::mt19937_64& rng, std::size_t trials) {
    ColoredTripletHdt full(tree, 4);
    expect_exception<std::invalid_argument>([&] { (void)full.extract({}); });
    expect_exception<std::invalid_argument>([&] { (void)full.contract({}); });
    expect_exception<std::invalid_argument>([&] { (void)full.extract({no_node}); });
    for (std::size_t trial = 0; trial < trials; ++trial) {
        std::vector<NodeId> kept;
        const std::size_t mask = trials == (std::size_t(1) << tree.leaf_count()) - 1 ? trial + 1 : rng();
        for (std::size_t i = 0; i < tree.leaf_count(); ++i)
            if (mask & (std::size_t(1) << i)) kept.push_back(tree.leaves()[i]);
        if (kept.empty()) kept.push_back(tree.leaves()[0]);
        CompressionStats extraction, contraction;
        auto extracted = full.extract(kept, &extraction);
        auto contracted = full.contract(kept, &contraction);
        std::size_t path_sum = 0;
        for (NodeId leaf : kept) {
            NodeId uid = full.hdt().leaf_for_vertex(leaf);
            while (uid != no_node) { ++path_sum; uid = full.hdt().node(uid).parent; }
        }
        CHECK(extraction.extraction_nodes <= 2 * path_sum);
        CHECK(extraction.contract_input_nodes <= extraction.extraction_nodes);
        CHECK(contraction.contract_input_nodes == tree.size());
        CHECK(extraction.queue_pops <= 5 * extraction.contract_input_nodes);
        CHECK(contraction.queue_pops <= 5 * contraction.contract_input_nodes);
        CHECK(extraction.merges + extracted.size() == extraction.contract_input_nodes);
        CHECK(contraction.merges + contracted.size() == contraction.contract_input_nodes);
        compare(tree, kept, std::move(extracted), rng);
        compare(tree, kept, std::move(contracted), rng);
        // Extraction must leave the source coloring and query unchanged.
        CHECK(full.query().resolved == 0 && full.query().unresolved == 0);
        for (NodeId leaf : tree.leaves()) CHECK(full.color(leaf) == 4);
    }
}

int main() try {
    std::mt19937_64 rng(5012013);
    // Original IDs may be huge and sparse. The optional external-ID lookup
    // must not allocate up to the largest ID. Exercise both validation paths.
    for (std::size_t n : {8, 300}) {
        const auto tree = star(n);
        auto atoms = ComponentTree::from_tree(tree).atoms();
        for (NodeId leaf : tree.leaves()) atoms[leaf].active_leaf = no_node - 1 - 1000003 * leaf;
        const auto valid = ComponentTree(atoms, tree.root());
        ColoredTripletHdt by_id(valid), by_atom(valid);
        const NodeId a = tree.leaves().front(), b = tree.leaves().back();
        by_id.recolor({{atoms[a].active_leaf, 7}, {atoms[b].active_leaf, 9}});
        by_atom.recolor_atoms({{a, 7}, {b, 9}});
        CHECK(by_id.query().resolved == by_atom.query().resolved);
        CHECK(by_id.query().unresolved == by_atom.query().unresolved);
        CHECK(by_id.color(atoms[b].active_leaf) == 9);
        const auto selected_id = by_id.extract({atoms[a].active_leaf, atoms[b].active_leaf});
        const auto selected_atom = by_atom.extract_atoms({a, b});
        CHECK(selected_id.size() == selected_atom.size());
        expect_exception<std::invalid_argument>([&] { by_atom.recolor_atoms({{a, 0}, {tree.root(), 2}}); });
        CHECK(by_atom.color(atoms[a].active_leaf) == 7);
        expect_exception<std::invalid_argument>([&] { (void)by_atom.extract_atoms({no_node}); });
        expect_exception<std::invalid_argument>([&] { (void)by_atom.contract_atoms({tree.root()}); });
        atoms[b].active_leaf = atoms[a].active_leaf;
        expect_exception<std::invalid_argument>([&] { (void)ComponentTree(atoms, tree.root()); });
    }
    for (std::size_t n = 1; n <= 12; ++n) {
        const std::size_t trials = n <= 6 ? (std::size_t(1) << n) - 1 : 20;
        check_subsets(star(n), rng, trials);
        check_subsets(caterpillar(n), rng, trials);
        for (int tree = 0; tree < 15; ++tree) check_subsets(random_tree(n, rng), rng, trials);
    }
    for (std::size_t n : {1000, 10000, 100000}) {
        const auto tree = caterpillar(n);
        ColoredTripletHdt full(tree);
        for (const auto& kept : {std::vector<NodeId>{tree.leaves().front()},
                                 std::vector<NodeId>{tree.leaves().front(), tree.leaves().back()},
                                 std::vector<NodeId>{tree.leaves()[n / 4], tree.leaves()[n / 2], tree.leaves().back()}}) {
            CompressionStats stats;
            const auto compressed = full.extract(kept, &stats);
            CHECK(compressed.size() <= 8 * kept.size());
            CHECK(stats.extraction_nodes <= 2 * kept.size() * (1 + full.hdt().node(full.hdt().root()).height));
            CHECK(stats.queue_pops <= 5 * stats.contract_input_nodes);
        }
    }
    std::cout << "compression: OK (subsets, recolorings, nested extracts; up to 100,000 leaves)\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
