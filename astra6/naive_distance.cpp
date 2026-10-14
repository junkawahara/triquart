// Independent test oracle for soda13.pdf, Section 2. LCA depths classify rooted
// triplets, and the four-point condition on unweighted path lengths classifies
// unrooted quartets. No HDT, decoration, coloring or single-tree DP is used.
#include "naive_distance.hpp"

#include <array>

namespace soda2013::oracle {
namespace {

NodeId lca(const Tree& tree, NodeId a, NodeId b) {
    while (tree.depth(a) > tree.depth(b)) a = tree.parent(a);
    while (tree.depth(b) > tree.depth(a)) b = tree.parent(b);
    while (a != b) {
        a = tree.parent(a);
        b = tree.parent(b);
    }
    return a;
}

std::size_t path_length(const Tree& tree, NodeId a, NodeId b) {
    return tree.depth(a) + tree.depth(b) - 2 * tree.depth(lca(tree, a, b));
}

void tally(DistanceCounts& result, int first, int second) {
    if (first == 0 && second == 0) result.E = add(result.E, 1);
    else if (first == 0) result.D = add(result.D, 1);
    else if (second == 0) result.C = add(result.C, 1);
    else if (first == second) result.A = add(result.A, 1);
    else result.B = add(result.B, 1);
}

std::vector<NodeId> matching_leaves(const Tree& first, const Tree& second) {
    validate_leaf_sets(first, second);
    std::vector<NodeId> result;
    for (NodeId v : first.leaves()) result.push_back(second.leaf(first.node(v).label));
    return result;
}

}  // namespace

int triplet_topology(const Tree& tree, NodeId a, NodeId b, NodeId c) {
    const std::array<std::size_t, 3> depths{
        tree.depth(lca(tree, a, b)), tree.depth(lca(tree, a, c)), tree.depth(lca(tree, b, c))};
    for (std::size_t i = 0; i < 3; ++i)
        if (depths[i] > depths[(i + 1) % 3] && depths[i] > depths[(i + 2) % 3])
            return static_cast<int>(i + 1);
    return 0;
}

int quartet_topology(const Tree& tree, NodeId a, NodeId b, NodeId c, NodeId d) {
    const std::array<std::size_t, 3> sums{
        path_length(tree, a, b) + path_length(tree, c, d),
        path_length(tree, a, c) + path_length(tree, b, d),
        path_length(tree, a, d) + path_length(tree, b, c)};
    for (std::size_t i = 0; i < 3; ++i)
        if (sums[i] < sums[(i + 1) % 3] && sums[i] < sums[(i + 2) % 3])
            return static_cast<int>(i + 1);
    return 0;
}

DistanceCounts triplet_counts(const Tree& first, const Tree& second) {
    const auto& a = first.leaves();
    const auto b = matching_leaves(first, second);
    DistanceCounts result;
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = i + 1; j < a.size(); ++j)
            for (std::size_t k = j + 1; k < a.size(); ++k)
                tally(result, triplet_topology(first, a[i], a[j], a[k]),
                      triplet_topology(second, b[i], b[j], b[k]));
    return result;
}

DistanceCounts quartet_counts(const Tree& first, const Tree& second) {
    const auto& a = first.leaves();
    const auto b = matching_leaves(first, second);
    DistanceCounts result;
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = i + 1; j < a.size(); ++j)
            for (std::size_t k = j + 1; k < a.size(); ++k)
                for (std::size_t l = k + 1; l < a.size(); ++l)
                    tally(result, quartet_topology(first, a[i], a[j], a[k], a[l]),
                          quartet_topology(second, b[i], b[j], b[k], b[l]));
    return result;
}

}  // namespace soda2013::oracle
