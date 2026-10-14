// ctree.hpp — the "working" copy of T2 on which the HDT is built.
//
// Sec. 5 of soda13.pdf keeps T2 (and its HDT) compressed during the
// recursion: `contract` / `extract` replace whole subtrees that only contain
// leaves of colour 0 by single nodes "with the counting information reset to
// represent the situation where all leaves in these subtrees have color 0".
// A CTree therefore has two kinds of leaves:
//
//   * designated (real) leaves: `leaf_id >= 0`, weight 1; they carry the
//     colour of the T1 leaf with the same id, and
//   * zero-leaves: `leaf_id == -1`, `weight >= 1`; they stand for `weight`
//     original leaves that all have colour 0 for the rest of the recursion.
//
// Internal nodes always have >= 2 children (unary nodes are contracted;
// they would not change any induced topology).
#pragma once

#include <cstdint>
#include <vector>

#include "tree.hpp"

namespace soda13 {

struct CTree {
    struct Node {
        int parent = -1;
        std::vector<int> children;
        int leaf_id = -1;         // >= 0 for designated leaves
        std::int64_t weight = 0;  // zero-leaves: number of merged colour-0 leaves; real leaves: 1
    };
    std::vector<Node> nodes;
    int root = -1;
    int num_designated = 0;  // number of real leaves

    bool is_leaf(int v) const { return nodes[v].children.empty(); }
    int num_nodes() const { return static_cast<int>(nodes.size()); }

    int add_node() {
        nodes.emplace_back();
        return static_cast<int>(nodes.size()) - 1;
    }
    int add_real_leaf(int parent, int leaf_id) {
        int v = add_node();
        nodes[v].leaf_id = leaf_id;
        nodes[v].weight = 1;
        attach(parent, v);
        ++num_designated;
        return v;
    }
    int add_zero_leaf(int parent, std::int64_t w) {
        int v = add_node();
        nodes[v].weight = w;
        attach(parent, v);
        return v;
    }
    void attach(int parent, int v) {
        nodes[v].parent = parent;
        if (parent >= 0) nodes[parent].children.push_back(v);
    }

    // Build from an input tree (all leaves designated; requires leaf ids).
    static CTree from_tree(const Tree& t) {
        CTree c;
        c.nodes.resize(t.nodes.size());
        c.root = t.root;
        for (int v = 0; v < t.num_nodes(); ++v) {
            c.nodes[v].parent = t.nodes[v].parent;
            c.nodes[v].children = t.nodes[v].children;
            if (t.is_leaf(v)) {
                c.nodes[v].leaf_id = t.nodes[v].leaf_id;
                c.nodes[v].weight = 1;
                ++c.num_designated;
            }
        }
        return c;
    }
};

}  // namespace soda13
