#pragma once

// Rooted input trees for Sections 2–5 of soda13.pdf. For quartet distances a
// root of degree two is harmless: induced unrooted topologies suppress it.
// Every other internal vertex has at least two children (Section 4, p. 6).
#include "counts.hpp"

#include <cstddef>
#include <limits>
#include <string>
#include <map>
#include <vector>

namespace soda2013 {

using NodeId = std::size_t;
inline constexpr NodeId no_node = std::numeric_limits<NodeId>::max();

struct Node {
    std::string label;
    std::vector<NodeId> children;
};

class Tree {
public:
    Tree(std::vector<Node> nodes, NodeId root);
    // Labels may be single-quoted (escape a quote by doubling it). Optional
    // internal labels, comments and finite branch lengths do not affect topology.
    static Tree from_newick(const std::string& text);

    NodeId root() const { return root_; }
    std::size_t size() const { return nodes_.size(); }
    std::size_t leaf_count() const { return leaves_.size(); }
    const Node& node(NodeId v) const { return nodes_.at(v); }
    bool is_leaf(NodeId v) const { return node(v).children.empty(); }
    NodeId parent(NodeId v) const { return parents_.at(v); }
    std::size_t depth(NodeId v) const { return depths_.at(v); }
    const std::vector<NodeId>& preorder() const { return preorder_; }
    const std::vector<NodeId>& leaves() const { return leaves_; }
    // Subtrees occupy [leaf_begin(v), leaf_end(v)) in leaves(), enabling
    // Section 5 recolorings without repeatedly traversing internal vertices.
    std::size_t leaf_begin(NodeId v) const { return begin_.at(v); }
    std::size_t leaf_end(NodeId v) const { return end_.at(v); }
    Count subtree_leaves(NodeId v) const { return end_.at(v) - begin_.at(v); }
    NodeId leaf(const std::string& label) const;

private:
    std::vector<Node> nodes_;
    NodeId root_;
    std::vector<NodeId> parents_, preorder_, leaves_;
    std::vector<std::size_t> depths_, begin_, end_;
    // Label matching is preprocessing only: O(n log n) comparisons, with no
    // hash worst case. The Count traversal then uses paired integer indices.
    std::map<std::string, NodeId> labels_;
};

// Section 2: both inputs must have the same nonempty set of leaf labels.
void validate_leaf_sets(const Tree& first, const Tree& second);

}  // namespace soda2013
