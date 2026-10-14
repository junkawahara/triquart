// tree.hpp — rooted input trees of arbitrary degree (Sec. 2 of soda13.pdf).
//
// The two input trees of every distance routine must share an identical,
// non-empty leaf-label set; `assign_leaf_ids` enforces this and gives every
// leaf an integer id 0..n-1 that is the same in both trees ("leaves in T1
// and T2 with the same labels are kept paired together", Sec. 5).
//
// Unary internal nodes do not change any induced triplet/quartet topology,
// so `contract_unary` removes them; all algorithms assume they are absent.
#pragma once

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace soda13 {

struct Tree {
    struct Node {
        int parent = -1;
        std::vector<int> children;
        std::string label;  // non-empty exactly for leaves
        int leaf_id = -1;   // shared id, set by assign_leaf_ids()
    };
    std::vector<Node> nodes;
    int root = -1;

    bool is_leaf(int v) const { return nodes[v].children.empty(); }
    int num_nodes() const { return static_cast<int>(nodes.size()); }
    int num_leaves() const;
    std::vector<int> leaves() const;  // in DFS order

    // Parse a Newick string, e.g. "((a,b),(c,d,e));".  Branch lengths and
    // internal labels are accepted and ignored.  Throws std::invalid_argument.
    static Tree from_newick(const std::string& s);
    std::string to_newick() const;

    // Remove internal nodes with exactly one child (and a root with one
    // child).  Keeps leaf labels and ids.
    void contract_unary();

    // Depth-first helpers.
    std::vector<int> subtree_leaf_counts() const;  // per node
    int degree_max() const;                        // max #children
};

// Verifies both trees have the same non-empty label set and assigns
// consistent leaf ids.  Throws std::invalid_argument otherwise.
void assign_leaf_ids(Tree& t1, Tree& t2);

}  // namespace soda13
