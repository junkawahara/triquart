// soda13.pdf, Section 5: extract the union of selected leaf-to-root HDT paths,
// replace every side subtree with an all-zero component, and contract greedily.
// The fragments below invert Figure 6 only on the visited part of the HDT.
// A G fragment is a sibling list; a C fragment has one top and a downward port.
#include "colored_hdt.hpp"

#include <utility>

namespace soda2013 {
namespace {

struct Fragment {
    NodeId first = no_node;
    NodeId last = no_node;
    NodeId bottom = no_node;
};

}  // namespace

ComponentTree ColoredHdt::extract(const std::vector<NodeId>& leaves, CompressionStats* stats) {
    std::vector<NodeId> atoms;
    atoms.reserve(leaves.size());
    for (NodeId leaf : leaves) atoms.push_back(original_atom(leaf));
    return extract_atoms(atoms, stats);
}

ComponentTree ColoredHdt::extract_atoms(const std::vector<NodeId>& leaves, CompressionStats* stats) {
    if (leaves.empty()) throw std::invalid_argument("extract requires a nonempty leaf set");
    for (NodeId leaf : leaves) validate_active_atom(leaf);
    for (NodeId leaf : leaves) {
        NodeId uid = hdt_.leaf_for_vertex(leaf);
        while (uid != no_node && !marked_[uid]) {
            marked_[uid] = true;
            uid = hdt_.node(uid).parent;
        }
    }

    // Includes one boundary atom for each side subtree, but does not descend it.
    std::vector<NodeId> order{hdt_.root()};
    for (std::size_t i = 0; i < order.size(); ++i) {
        const NodeId uid = order[i];
        scratch_index_[uid] = i;
        const auto& node = hdt_.node(uid);
        if (marked_[uid] && node.left != no_node) {
            order.push_back(node.left);
            order.push_back(node.right);
        }
    }
    if (stats) stats->extraction_nodes += order.size();
    std::vector<Fragment> fragments(order.size());
    std::vector<ComponentAtom> atoms;
    std::vector<NodeId> next, first_child, last_child;
    const auto attach_forest = [&](NodeId parent, const Fragment& forest) {
        if (last_child[parent] == no_node) first_child[parent] = forest.first;
        else next[last_child[parent]] = forest.first;
        last_child[parent] = forest.last;
    };
    for (std::size_t i = order.size(); i-- > 0;) {
        const NodeId uid = order[i];
        const auto& node = hdt_.node(uid);
        auto& fragment = fragments[i];
        if (!marked_[uid] || node.left == no_node) {
            const NodeId active = marked_[uid] ? tree_.atom(node.source).active_leaf : no_node;
            const NodeId v = atoms.size();
            atoms.push_back({node.type, node.leaf_mass, active, {}});
            next.push_back(no_node);
            first_child.push_back(no_node);
            last_child.push_back(no_node);
            fragment = {v, v, node.downwards_open ? v : no_node};
        } else {
            const auto& left = fragments[scratch_index_[node.left]];
            const auto& right = fragments[scratch_index_[node.right]];
            switch (node.composition) {
            case Composition::CC:
                if (right.bottom == no_node) throw std::logic_error("closed upper CC fragment");
                attach_forest(right.bottom, left);
                fragment = {right.first, right.last, left.bottom};
                break;
            case Composition::GG:
                next[left.last] = right.first;
                fragment = {left.first, right.last, no_node};
                break;
            case Composition::IG:
                attach_forest(left.first, right);
                fragment = {left.first, left.last, node.downwards_open ? left.first : no_node};
                break;
            default:
                throw std::logic_error("noncomposition inside extract path union");
            }
            if (node.c_to_g) fragment.bottom = no_node;
        }
        marked_[uid] = false;
    }
    for (NodeId v = 0; v < atoms.size(); ++v)
        for (NodeId child = first_child[v]; child != no_node; child = next[child])
            atoms[v].children.push_back(child);
    const NodeId root = fragments[0].first;
    if (root != fragments[0].last) throw std::logic_error("extraction produced a root forest");
    return contract_components(ComponentTree(std::move(atoms), root), stats);
}

ComponentTree ColoredHdt::contract(const std::vector<NodeId>& leaves,
                                        CompressionStats* stats) const {
    std::vector<NodeId> atoms;
    atoms.reserve(leaves.size());
    for (NodeId leaf : leaves) atoms.push_back(original_atom(leaf));
    return contract_atoms(atoms, stats);
}

ComponentTree ColoredHdt::contract_atoms(const std::vector<NodeId>& leaves, CompressionStats* stats) const {
    if (leaves.empty()) throw std::invalid_argument("contract requires a nonempty leaf set");
    std::vector<bool> keep(tree_.size(), false);
    for (NodeId leaf : leaves) {
        validate_active_atom(leaf);
        keep[leaf] = true;
    }
    auto atoms = tree_.atoms();
    for (NodeId v = 0; v < atoms.size(); ++v)
        if (!keep[v]) atoms[v].active_leaf = no_node;
    return contract_components(ComponentTree(std::move(atoms), tree_.root()), stats);
}

}  // namespace soda2013
