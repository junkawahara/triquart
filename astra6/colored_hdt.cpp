// soda13.pdf, Section 5: batched recoloring uses a marked-only traversal of
// the HDT. Its cost depends on the union of paths, not on the full tree size.
#include "colored_hdt.hpp"

#include <algorithm>

namespace soda2013 {

ColoredHdt::ColoredHdt(ComponentTree tree, Color initial_color)
    : tree_(std::move(tree)), hdt_(tree_), colors_(tree_.size(), initial_color),
      marked_(hdt_.nodes().size(), false), scratch_index_(hdt_.nodes().size()) {}

NodeId ColoredHdt::original_atom(NodeId leaf) const {
    if (!leaf_index_built_) {
        for (NodeId v = 0; v < tree_.size(); ++v)
            if (tree_.atom(v).active_leaf != no_node) leaf_atoms_.push_back({tree_.atom(v).active_leaf, v});
        std::sort(leaf_atoms_.begin(), leaf_atoms_.end());
        leaf_index_built_ = true;
    }
    const auto found = std::lower_bound(leaf_atoms_.begin(), leaf_atoms_.end(), std::make_pair(leaf, NodeId(0)));
    if (found == leaf_atoms_.end() || found->first != leaf)
        throw std::invalid_argument("unknown original leaf ID");
    return found->second;
}

void ColoredHdt::validate_active_atom(NodeId atom) const {
    if (atom >= tree_.size() || tree_.atom(atom).active_leaf == no_node)
        throw std::invalid_argument("only active leaf atoms can be selected");
}

void ColoredHdt::recolor(const std::vector<Recoloring>& changes) {
    std::vector<Recoloring> atoms;
    atoms.reserve(changes.size());
    for (const auto& change : changes) atoms.push_back({original_atom(change.leaf), change.color});
    recolor_atoms(atoms);
}

void ColoredHdt::recolor_atoms(const std::vector<Recoloring>& changes) {
    for (const auto& change : changes)
        validate_active_atom(change.leaf);
    ++stats_.batches;
    for (const auto& change : changes) {
        const NodeId atom = change.leaf;
        if (colors_[atom] == change.color) continue;
        colors_[atom] = change.color;
        NodeId uid = hdt_.leaf_for_vertex(atom);
        while (uid != no_node && !marked_[uid]) {
            marked_[uid] = true;
            uid = hdt_.node(uid).parent;
        }
    }
    if (!marked_[hdt_.root()]) return;
    std::vector<NodeId> order{hdt_.root()};
    for (std::size_t i = 0; i < order.size(); ++i) {
        const auto& node = hdt_.node(order[i]);
        if (node.left != no_node && marked_[node.left]) order.push_back(node.left);
        if (node.right != no_node && marked_[node.right]) order.push_back(node.right);
    }
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        recompute(*it);
        marked_[*it] = false;
    }
}

}  // namespace soda2013
