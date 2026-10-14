#pragma once

// soda13.pdf, Section 5: the common coloring, path-union and compression
// machinery. Decorations supply only the bottom-up recomputation of a node.
#include "hdt.hpp"

#include <utility>

namespace soda2013 {

using Color = std::size_t;
struct Recoloring { NodeId leaf; Color color; };

struct ColorUpdateStats {
    std::size_t batches = 0;
    std::size_t recomputed_components = 0;
    std::size_t color_rows = 0;
    std::size_t color_pairs = 0;
};

class ColoredHdt {
public:
    virtual ~ColoredHdt() = default;
    // Mark the UNION of changed leaf paths, then recompute each node once.
    // Validate all leaf IDs before mutation; repeated entries use the last color.
    void recolor(const std::vector<Recoloring>& changes);
    // Section 5's paired-leaf pointers are component indices in the production
    // traversal. These variants accept CURRENT component atom IDs and perform
    // O(1) address lookup, with no dictionary construction or hash assumptions.
    void recolor_atoms(const std::vector<Recoloring>& changes);
    const Hdt& hdt() const { return hdt_; }
    Color color(NodeId leaf) const { return colors_.at(original_atom(leaf)); }
    const ColorUpdateStats& stats() const { return stats_; }
    const ComponentTree& components() const { return tree_; }
    // Side components become all-zero summaries; selected leaves retain their
    // original IDs. extract visits paths and frontiers, contract scans the tree.
    ComponentTree extract(const std::vector<NodeId>& leaves, CompressionStats* stats = nullptr);
    ComponentTree contract(const std::vector<NodeId>& leaves, CompressionStats* stats = nullptr) const;
    ComponentTree extract_atoms(const std::vector<NodeId>& atoms, CompressionStats* stats = nullptr);
    ComponentTree contract_atoms(const std::vector<NodeId>& atoms, CompressionStats* stats = nullptr) const;

protected:
    explicit ColoredHdt(ComponentTree tree, Color initial_color);
    ComponentTree tree_;
    Hdt hdt_;
    // Optional original-ID lookup, built in O(m + k log k) for m atoms/k leaves.
    // It is not used by Count, which maintains direct paired-leaf indices.
    mutable std::vector<std::pair<NodeId, NodeId>> leaf_atoms_;
    mutable bool leaf_index_built_ = false;
    std::vector<Color> colors_;
    std::vector<bool> marked_;
    std::vector<std::size_t> scratch_index_;
    ColorUpdateStats stats_;
    NodeId original_atom(NodeId leaf) const;
    void validate_active_atom(NodeId atom) const;
    virtual void recompute(NodeId uid) = 0;
};

}  // namespace soda2013
