// Semantic oracle for every Figure 10–12 triplet counter, not just the root
// answer. Enumerate original leaves, LCAs and branch points independently of
// the HDT composition formulas, for arbitrary colorings including color zero.
#include "naive_distance.hpp"
#include "test_support.hpp"
#include "triplet_counters.hpp"

#include <array>
#include <map>
#include <set>

using namespace soda2013;
using namespace test_support;

NodeId lca(const Tree& tree, NodeId a, NodeId b) {
    while (tree.depth(a) > tree.depth(b)) a = tree.parent(a);
    while (tree.depth(b) > tree.depth(a)) b = tree.parent(b);
    while (a != b) { a = tree.parent(a); b = tree.parent(b); }
    return a;
}

void check_counters(const Tree& tree, const ColoredTripletHdt& colored) {
    const auto& hdt = colored.hdt();
    std::vector<std::vector<NodeId>> leaves(hdt.nodes().size());
    std::vector<std::set<NodeId>> vertices(hdt.nodes().size()), paths(hdt.nodes().size());
    for (NodeId uid = 0; uid < hdt.nodes().size(); ++uid) {
        const auto& node = hdt.node(uid);
        if (node.left == no_node) {
            vertices[uid].insert(node.source);
            paths[uid].insert(node.source);
            if (tree.is_leaf(node.source)) leaves[uid].push_back(node.source);
        } else {
            for (NodeId child : {node.left, node.right}) {
                leaves[uid].insert(leaves[uid].end(), leaves[child].begin(), leaves[child].end());
                vertices[uid].insert(vertices[child].begin(), vertices[child].end());
            }
            if (node.composition == Composition::CC) {
                paths[uid] = paths[node.left];
                paths[uid].insert(paths[node.right].begin(), paths[node.right].end());
            } else if (node.composition == Composition::IG) {
                paths[uid] = vertices[node.left];
            }
        }
        const auto& actual = colored.counters(uid);
        // Each array holds n_i, n_(ii), n_i•, n_i↑• as independently defined.
        std::map<Color, std::array<Count, 4>> expected;
        Count zero = 0;
        for (NodeId leaf : leaves[uid]) {
            if (colored.color(leaf)) ++expected[colored.color(leaf)][0];
            else ++zero;
        }
        NodeId super_root = no_node;
        if (node.type == ComponentType::G) {
            for (NodeId v : vertices[uid]) {
                if (!vertices[uid].count(tree.parent(v))) {
                    super_root = tree.parent(v);
                    break;
                }
            }
            CHECK(super_root != no_node);
        }
        const auto branch = [&](NodeId leaf) {
            while (!paths[uid].count(leaf)) {
                leaf = tree.parent(leaf);
                CHECK(leaf != no_node);
            }
            return tree.depth(leaf);
        };
        for (std::size_t i = 0; i < leaves[uid].size(); ++i) {
            for (std::size_t j = i + 1; j < leaves[uid].size(); ++j) {
                const NodeId a = leaves[uid][i], b = leaves[uid][j];
                const Color ca = colored.color(a), cb = colored.color(b);
                if (!ca || !cb) continue;
                const NodeId ancestor = lca(tree, a, b);
                const bool at = node.type == ComponentType::G ? ancestor == super_root :
                                paths[uid].count(ancestor) && branch(a) == branch(b);
                const bool inside = node.type == ComponentType::G ? ancestor != super_root :
                                    !paths[uid].count(ancestor);
                if (ca == cb && inside) ++expected[ca][1];
                if (ca != cb && at) { ++expected[ca][2]; ++expected[cb][2]; }
                if (ca != cb && node.type == ComponentType::C) {
                    if (branch(a) > branch(b)) ++expected[ca][3];
                    if (branch(b) > branch(a)) ++expected[cb][3];
                }
            }
        }
        Count twice_at = 0;
        CHECK(actual.zero == zero);
        CHECK(actual.positive == leaves[uid].size() - zero);
        CHECK(actual.colors.size() == expected.size());
        auto it = expected.begin();
        for (const auto& row : actual.colors) {
            CHECK(row.color == it->first);
            CHECK(row.leaves == it->second[0]);
            CHECK(row.same_inside == it->second[1]);
            CHECK(row.different_at == it->second[2]);
            CHECK(row.i_up_other == it->second[3]);
            twice_at += it->second[2];
            ++it;
        }
        CHECK(actual.different_at == twice_at / 2);
        SharedTriplets shared;
        for (std::size_t i = 0; i < leaves[uid].size(); ++i)
            for (std::size_t j = i + 1; j < leaves[uid].size(); ++j)
                for (std::size_t k = j + 1; k < leaves[uid].size(); ++k) {
                    const NodeId a = leaves[uid][i], b = leaves[uid][j], c = leaves[uid][k];
                    const Color ca = colored.color(a), cb = colored.color(b), cc = colored.color(c);
                    if (!ca || !cb || !cc) continue;
                    const int topology = oracle::triplet_topology(tree, a, b, c);
                    if (topology == 0 && ca != cb && ca != cc && cb != cc) ++shared.unresolved;
                    if ((topology == 1 && ca == cb && ca != cc) ||
                        (topology == 2 && ca == cc && ca != cb) ||
                        (topology == 3 && cb == cc && cb != ca)) ++shared.resolved;
                }
        CHECK(actual.shared.resolved == shared.resolved);
        CHECK(actual.shared.unresolved == shared.unresolved);
    }
}

int main() try {
    std::mt19937_64 rng(7012013);
    for (std::size_t n = 1; n <= 12; ++n) {
        for (int trial = 0; trial < 50; ++trial) {
            const auto tree = random_tree(n, rng, trial % 2 ? 2 : 12);
            ColoredTripletHdt colored(tree);
            check_counters(tree, colored);
            for (int batch = 0; batch < 10; ++batch) {
                std::vector<Recoloring> changes;
                std::set<NodeId> dirty;
                for (NodeId leaf : tree.leaves()) {
                    if (rng() % 2) continue;
                    const Color color = rng() % (n + 1);
                    changes.push_back({leaf, color});
                    if (color == colored.color(leaf)) continue;
                    NodeId uid = colored.hdt().leaf_for_vertex(leaf);
                    while (uid != no_node) { dirty.insert(uid); uid = colored.hdt().node(uid).parent; }
                }
                const auto before = colored.stats().recomputed_components;
                colored.recolor(changes);
                CHECK(colored.stats().recomputed_components - before == dirty.size());
                check_counters(tree, colored);
            }
            const NodeId leaf = tree.leaves()[0];
            colored.recolor({{leaf, 0}, {leaf, 1000000}, {leaf, 2000000}});
            CHECK(colored.color(leaf) == 2000000);
            check_counters(tree, colored);
            const auto before = colored.stats().recomputed_components;
            colored.recolor({});
            colored.recolor({{leaf, 2000000}});
            CHECK(colored.stats().recomputed_components == before);
            expect_exception<std::invalid_argument>([&] { colored.recolor({{no_node, 5}}); });
            if (!tree.is_leaf(tree.root())) {
                expect_exception<std::invalid_argument>([&] {
                    colored.recolor({{leaf, 0}, {tree.root(), 5}});
                });
                CHECK(colored.color(leaf) == 2000000);
            }
        }
    }
    std::cout << "triplet counters: OK (every component; 600 trees, 6,000 recoloring batches)\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
