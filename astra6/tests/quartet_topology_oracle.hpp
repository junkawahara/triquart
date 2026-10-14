#pragma once

// Test-only interpretation of Figures 4, 8 and 9, from original-tree paths.
#include "naive_distance.hpp"
#include "quartet_contributions.hpp"
#include "quartet_counter_oracle.hpp"

namespace quartet_oracle {
using namespace soda2013;

struct Topology {
    int split = 0;  // 0 means unresolved (or an incompatible color pattern).
    QuartetType type = QuartetType::alpha;
    NodeId anchor = no_node;
};

inline Topology topology(const Tree& tree, const std::array<NodeId, 4>& leaves) {
    const int split = oracle::quartet_topology(tree, leaves[0], leaves[1], leaves[2], leaves[3]);
    if (!split) return {};
    NodeId root = leaves[0];
    for (NodeId leaf : leaves) root = counter_oracle::lca(tree, root, leaf);
    std::array<NodeId, 4> children{no_node, no_node, no_node, no_node};
    std::array<std::size_t, 4> sizes{};
    std::size_t groups = 0;
    for (NodeId leaf : leaves) {
        NodeId child = leaf;
        while (tree.parent(child) != root) child = tree.parent(child);
        std::size_t i = 0;
        while (i < groups && children[i] != child) ++i;
        if (i == groups) children[groups++] = child;
        ++sizes[i];
    }
    if (groups == 3) return {split, QuartetType::beta, root};
    CHECK(groups == 2);
    if (sizes[0] == 2) return {split, QuartetType::alpha, root};
    const NodeId triple_child = children[sizes[0] == 3 ? 0 : 1];
    NodeId anchor = no_node;
    for (NodeId leaf : leaves) {
        NodeId child = leaf;
        while (tree.parent(child) != root) child = tree.parent(child);
        if (child == triple_child) anchor = anchor == no_node ? leaf : counter_oracle::lca(tree, anchor, leaf);
    }
    return {split, QuartetType::gamma, anchor};
}

inline int split_for_pair(std::size_t first, std::size_t second) {
    if ((first == 0 && second == 1) || (first == 2 && second == 3)) return 1;
    if ((first == 0 && second == 2) || (first == 1 && second == 3)) return 2;
    return 3;
}

inline Topology coloring(const std::vector<Color>& colors, const std::array<NodeId, 4>& leaves) {
    std::array<Color, 4> used{};
    std::array<std::size_t, 4> sizes{};
    std::size_t count = 0, zeros = 0;
    for (NodeId leaf : leaves) {
        if (colors[leaf] == 0) { ++zeros; continue; }
        std::size_t i = 0;
        while (i < count && used[i] != colors[leaf]) ++i;
        if (i == count) used[count++] = colors[leaf];
        ++sizes[i];
    }
    QuartetType type;
    if (!zeros && count == 2 && sizes[0] == 2) type = QuartetType::alpha;
    else if (!zeros && count == 3) type = QuartetType::beta;
    else if (zeros == 1 && count == 2) type = QuartetType::gamma;
    else return {};
    for (std::size_t i = 0; i < count; ++i) {
        if (sizes[i] != 2) continue;
        for (std::size_t a = 0; a < 4; ++a)
            for (std::size_t b = a + 1; b < 4; ++b)
                if (colors[leaves[a]] == used[i] && colors[leaves[b]] == used[i])
                    return {split_for_pair(a, b), type, no_node};
    }
    throw std::logic_error("color pattern without repeated pair");
}

inline QuartetCases resolved_cases(const Tree& first, const Tree& second, bool disagreement) {
    validate_leaf_sets(first, second);
    const auto& a = first.leaves();
    std::vector<NodeId> b;
    for (NodeId leaf : a) b.push_back(second.leaf(first.node(leaf).label));
    QuartetCases result;
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = i + 1; j < a.size(); ++j)
            for (std::size_t k = j + 1; k < a.size(); ++k)
                for (std::size_t l = k + 1; l < a.size(); ++l) {
                    const auto x = topology(first, {a[i], a[j], a[k], a[l]});
                    const auto y = topology(second, {b[i], b[j], b[k], b[l]});
                    if (x.split && y.split && ((x.split != y.split) == disagreement)) ++result.at(x.type, y.type);
                }
    return result;
}

inline QuartetCases shared_cases(const Tree& first, const Tree& second) {
    return resolved_cases(first, second, false);
}
inline QuartetCases disagreement_cases(const Tree& first, const Tree& second) {
    return resolved_cases(first, second, true);
}

}  // namespace quartet_oracle
