#pragma once

// soda13.pdf, Section 4 / Figures 5–6. Transformations do not add HDT nodes;
// each binary composition adds exactly one. C1 is below C2 in CC -> C.
#include "component_tree.hpp"

namespace soda2013 {

enum class Composition { LeafToC, IdentityI, CC, GG, IG, Atom };

struct HdtNode {
    ComponentType type = ComponentType::I;
    Composition composition = Composition::IdentityI;
    // A final C -> G transformation, after the composition (or L -> C).
    bool c_to_g = false;
    bool downwards_open = false;
    NodeId left = no_node;
    NodeId right = no_node;
    NodeId parent = no_node;
    NodeId source = no_node;  // Original vertex at an HDT leaf.
    std::size_t mass = 1;     // Original vertices / HDT leaves, not colored leaves.
    std::size_t height = 0;
    Count leaf_mass = 0;      // Original leaves, including frozen color-zero leaves.
};

struct HdtBuildStats {
    std::vector<std::size_t> round_sizes;
    std::size_t visited_components = 0;
    std::size_t leaf_to_c = 0;
    std::size_t c_to_g = 0;
    std::size_t cc = 0;
    std::size_t gg = 0;
    std::size_t ig = 0;
};

class Hdt {
public:
    explicit Hdt(const Tree& tree);
    explicit Hdt(const ComponentTree& tree);
    NodeId root() const { return root_; }
    const std::vector<HdtNode>& nodes() const { return nodes_; }
    const HdtNode& node(NodeId uid) const { return nodes_.at(uid); }
    NodeId leaf_for_vertex(NodeId v) const { return vertex_to_leaf_.at(v); }
    const HdtBuildStats& stats() const { return stats_; }

private:
    std::vector<HdtNode> nodes_;
    std::vector<NodeId> vertex_to_leaf_;
    NodeId root_ = no_node;
    HdtBuildStats stats_;
};

}  // namespace soda2013
