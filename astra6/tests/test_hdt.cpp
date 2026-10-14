// Structural rather than numerical validation of Section 4 / Lemma 4.1.
// On small inputs explicitly expand EVERY component's vertex set and inspect
// its original-tree boundary. On large inputs check work and balance bounds.
#include "hdt.hpp"
#include "test_support.hpp"

#include <cmath>
#include <set>

using namespace soda2013;
using namespace test_support;

void check_hdt(const Tree& tree, bool expand, bool* covered = nullptr) {
    const Hdt hdt(tree);
    const auto& nodes = hdt.nodes();
    CHECK(nodes.size() == tree.size() * 2 - 1);
    CHECK(nodes[hdt.root()].type == ComponentType::C);
    CHECK(nodes[hdt.root()].mass == tree.size());
    CHECK(nodes[hdt.root()].parent == no_node);
    const auto& stats = hdt.stats();
    CHECK(stats.leaf_to_c == tree.leaf_count());
    CHECK(stats.cc + stats.gg + stats.ig == tree.size() - 1);
    CHECK(stats.visited_components <= 9 * tree.size());
    for (std::size_t i = 1; i < stats.round_sizes.size(); ++i) {
        const auto old_size = stats.round_sizes[i - 1], new_size = stats.round_sizes[i];
        CHECK(new_size < old_size);
        CHECK(9 * (old_size - new_size) >= old_size);
    }

    std::vector<std::set<NodeId>> vertices(expand ? nodes.size() : 0);
    for (NodeId uid = 0; uid < nodes.size(); ++uid) {
        const auto& node = nodes[uid];
        if (node.left == no_node) {
            CHECK(node.right == no_node);
            CHECK(node.source < tree.size());
            CHECK(hdt.leaf_for_vertex(node.source) == uid);
            CHECK(node.mass == 1 && node.height == 0);
            if (expand) vertices[uid].insert(node.source);
        } else {
            CHECK(node.left < uid && node.right < uid);
            CHECK(nodes[node.left].parent == uid && nodes[node.right].parent == uid);
            CHECK(node.mass == nodes[node.left].mass + nodes[node.right].mass);
            CHECK(node.height == 1 + std::max(nodes[node.left].height, nodes[node.right].height));
            const auto left_type = nodes[node.left].type, right_type = nodes[node.right].type;
            if (node.composition == Composition::CC)
                CHECK(left_type == ComponentType::C && right_type == ComponentType::C);
            else if (node.composition == Composition::GG)
                CHECK(left_type == ComponentType::G && right_type == ComponentType::G);
            else {
                CHECK(node.composition == Composition::IG);
                CHECK(left_type == ComponentType::I && right_type == ComponentType::G);
            }
            if (expand) {
                vertices[uid] = vertices[node.left];
                for (NodeId v : vertices[node.right]) CHECK(vertices[uid].insert(v).second);
            }
        }
        const auto height_bound = std::ceil(std::log(static_cast<double>(node.mass)) /
                                             std::log(10.0 / 9.0));
        CHECK(static_cast<double>(node.height) <= height_bound);
        if (node.c_to_g) CHECK(node.type == ComponentType::G);
        if (covered) {
            covered[static_cast<int>(node.composition)] = true;
            covered[5] = covered[5] || node.c_to_g;
        }
        if (!expand) continue;
        const auto& component = vertices[uid];
        CHECK(component.size() == node.mass);
        std::vector<NodeId> tops;
        std::size_t downward = 0;
        for (NodeId v : component) {
            if (tree.parent(v) == no_node || !component.count(tree.parent(v))) tops.push_back(v);
            for (NodeId child : tree.node(v).children)
                if (!component.count(child)) ++downward;
        }
        if (node.type == ComponentType::I) {
            CHECK(component.size() == 1 && !tree.is_leaf(*component.begin()));
        } else if (node.type == ComponentType::C) {
            CHECK(tops.size() == 1 && downward <= 1);
        } else {
            CHECK(node.type == ComponentType::G && downward == 0 && !tops.empty());
            const NodeId super_root = tree.parent(tops[0]);
            CHECK(super_root != no_node);
            for (NodeId top : tops) CHECK(tree.parent(top) == super_root);
        }
        if (node.c_to_g) {
            // The component was a connected, closed C before the transform.
            CHECK(tops.size() == 1 && downward == 0);
        }
        if (node.composition == Composition::CC) {
            bool correct_orientation = false;
            for (NodeId v : vertices[node.left])
                if (vertices[node.right].count(tree.parent(v))) correct_orientation = true;
            CHECK(correct_orientation);
        }
    }
}

int main() try {
    std::mt19937_64 rng(4012013);
    bool covered[6]{};
    for (std::size_t n = 1; n <= 100; ++n) {
        check_hdt(star(n), true, covered);
        check_hdt(caterpillar(n), true, covered);
        for (int trial = 0; trial < 20; ++trial)
            check_hdt(random_tree(n, rng, trial % 2 ? 2 : 12), true, covered);
    }
    for (bool action : covered) CHECK(action);
    for (std::size_t n : {1000, 10000, 100000}) {
        check_hdt(star(n), false);
        check_hdt(caterpillar(n), false);
    }
    std::cout << "HDT: OK (all five actions, 2,200 structural cases; up to 100,000 leaves)\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
