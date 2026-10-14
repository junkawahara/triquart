#pragma once

// Section 5 of soda13.pdf: a contracted T2 is a tree of Section 4 components.
// A frozen all-zero C/G component needs only its number of leaves: every other
// Figure 11/14 counter contains at least one positive color and is zero.
#include "tree.hpp"

namespace soda2013 {

enum class ComponentType { L, I, C, G };

struct ComponentAtom {
    ComponentType type = ComponentType::I;
    Count leaves = 0;
    NodeId active_leaf = no_node;  // Stable original T2 leaf ID, if incontractible.
    std::vector<NodeId> children;
};

class ComponentTree {
public:
    ComponentTree(std::vector<ComponentAtom> atoms, NodeId root);
    static ComponentTree from_tree(const Tree& tree);
    std::size_t size() const { return atoms_.size(); }
    NodeId root() const { return root_; }
    const ComponentAtom& atom(NodeId v) const { return atoms_.at(v); }
    const std::vector<ComponentAtom>& atoms() const { return atoms_; }
    NodeId parent(NodeId v) const { return parents_.at(v); }

private:
    std::vector<ComponentAtom> atoms_;
    NodeId root_;
    std::vector<NodeId> parents_;
};

struct CompressionStats {
    std::size_t extraction_nodes = 0;
    std::size_t contract_input_nodes = 0;
    std::size_t contract_output_nodes = 0;
    std::size_t queue_pops = 0;
    std::size_t merges = 0;
};

// Section 5: greedily apply the Section 4 actions using a queue. A merge may
// never consume a designated active leaf. All other atoms represent color zero.
ComponentTree contract_components(const ComponentTree& tree, CompressionStats* stats = nullptr);

}  // namespace soda2013
