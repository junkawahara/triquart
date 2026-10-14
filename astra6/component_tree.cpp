// soda13.pdf, Section 5 (contract) using precisely the five Section 4 actions.
// Intrusive child and zero-G sibling lists let every queue visit/merge take O(1).
#include "component_tree.hpp"

#include <algorithm>
#include <array>
#include <deque>
#include <utility>

namespace soda2013 {
namespace {

// Linear validation of fixed-width leaf IDs. ComponentTree is reconstructed
// inside Section 5, so a comparison sort here would add an unwanted log factor.
// At most ceil(bit-width(NodeId)/8) radix passes suffice. For <=256 keys, sort
// is also linear with a fixed bound of eight comparisons per key up to constants.
void validate_distinct_ids(std::vector<NodeId> ids) {
    if (ids.size() <= 256) std::sort(ids.begin(), ids.end());
    else {
        const NodeId maximum = *std::max_element(ids.begin(), ids.end());
        std::vector<NodeId> scratch(ids.size());
        for (unsigned shift = 0; shift < std::numeric_limits<NodeId>::digits; shift += 8) {
            if ((maximum >> shift) == 0) break;
            std::array<std::size_t, 256> positions{};
            for (NodeId id : ids) ++positions[(id >> shift) & 255];
            std::size_t start = 0;
            for (auto& position : positions) {
                const auto frequency = position;
                position = start;
                start += frequency;
            }
            for (NodeId id : ids) scratch[positions[(id >> shift) & 255]++] = id;
            ids.swap(scratch);
        }
    }
    if (std::adjacent_find(ids.begin(), ids.end()) != ids.end())
        throw std::invalid_argument("duplicate designated component leaf");
}

}  // namespace

ComponentTree::ComponentTree(std::vector<ComponentAtom> atoms, NodeId root)
    : atoms_(std::move(atoms)), root_(root), parents_(atoms_.size(), no_node) {
    if (atoms_.empty() || root_ >= size()) throw std::invalid_argument("invalid component root");
    std::vector<bool> seen(size(), false);
    std::vector<NodeId> stack{root_};
    std::vector<NodeId> labels;
    seen[root_] = true;
    std::size_t visited = 0;
    while (!stack.empty()) {
        const NodeId v = stack.back();
        stack.pop_back();
        ++visited;
        const auto& atom = atoms_[v];
        if (atom.type == ComponentType::L) throw std::invalid_argument("apply L -> C to component atoms");
        if (atom.type == ComponentType::I && (atom.leaves || atom.children.empty()))
            throw std::invalid_argument("I atoms are internal vertices without leaves");
        if (atom.type == ComponentType::C && atom.children.size() > 1)
            throw std::invalid_argument("C atom has more than one downward edge");
        if (atom.type == ComponentType::G &&
            (!atom.children.empty() || parents_[v] == no_node || atoms_[parents_[v]].type != ComponentType::I))
            throw std::invalid_argument("G atom must be a closed child of I");
        if (atom.active_leaf != no_node) {
            if (atom.leaves != 1 || !atom.children.empty())
                throw std::invalid_argument("invalid designated component leaf");
            labels.push_back(atom.active_leaf);
        }
        for (NodeId child : atom.children) {
            if (child >= size() || seen[child]) throw std::invalid_argument("invalid component edges");
            seen[child] = true;
            parents_[child] = v;
            stack.push_back(child);
        }
    }
    if (visited != size()) throw std::invalid_argument("unreachable component atom");
    validate_distinct_ids(std::move(labels));
}

ComponentTree ComponentTree::from_tree(const Tree& tree) {
    std::vector<ComponentAtom> atoms;
    atoms.reserve(tree.size());
    for (NodeId v = 0; v < tree.size(); ++v) {
        if (tree.is_leaf(v)) atoms.push_back({ComponentType::C, 1, v, {}});
        else atoms.push_back({ComponentType::I, 0, no_node, tree.node(v).children});
    }
    return ComponentTree(std::move(atoms), tree.root());
}

namespace {

struct MutableAtom {
    ComponentType type;
    Count leaves;
    NodeId active_leaf;
    bool alive = true;
    bool queued = false;
    NodeId parent = no_node, first = no_node, last = no_node;
    NodeId previous = no_node, next = no_node;
    NodeId zero_previous = no_node, zero_next = no_node, first_zero = no_node;
    std::size_t children = 0, groups = 0, zero_groups = 0;
};

class Contractor {
public:
    Contractor(const ComponentTree& tree, CompressionStats* stats) : root_(tree.root()), stats_(stats) {
        nodes_.reserve(tree.size());
        for (const auto& atom : tree.atoms()) nodes_.push_back({atom.type, atom.leaves, atom.active_leaf});
        for (NodeId v = 0; v < tree.size(); ++v)
            for (NodeId child : tree.atom(v).children) attach(v, child);
        // Normalize every closed C below I before considering any composition.
        for (NodeId v = 0; v < tree.size(); ++v) normalize(v);
        for (NodeId v = 0; v < tree.size(); ++v) enqueue(v);
        if (stats_) stats_->contract_input_nodes += tree.size();
    }

    ComponentTree run() {
        while (!queue_.empty()) {
            const NodeId v = queue_.front();
            queue_.pop_front();
            nodes_[v].queued = false;
            if (stats_) ++stats_->queue_pops;
            if (!nodes_[v].alive) continue;
            normalize(v);
            auto& node = nodes_[v];
            if (node.type == ComponentType::I && node.zero_groups >= 2) {
                const NodeId a = node.first_zero, b = nodes_[a].zero_next;
                nodes_[a].leaves = add(nodes_[a].leaves, nodes_[b].leaves);
                detach(b);
                nodes_[b].alive = false;
                changed(v);  // GG -> G; I remains available for another pair.
            } else if (node.type == ComponentType::I && node.groups == 1 &&
                       node.zero_groups == 1 && node.children - node.groups <= 1) {
                const NodeId g = node.first_zero;
                node.leaves = add(node.leaves, nodes_[g].leaves);
                detach(g);
                nodes_[g].alive = false;
                node.type = ComponentType::C;
                normalize(v);
                changed(v);  // IG -> C, possibly followed by C -> G.
            } else if (node.type == ComponentType::C && node.active_leaf == no_node &&
                       node.first != no_node && nodes_[node.first].type == ComponentType::C &&
                       nodes_[node.first].active_leaf == no_node) {
                const NodeId lower = node.first, child = nodes_[lower].first;
                node.leaves = add(node.leaves, nodes_[lower].leaves);
                detach(lower);
                if (child != no_node) { detach(child); attach(v, child); }
                nodes_[lower].alive = false;
                normalize(v);
                changed(v);  // CC -> C, lower is C1.
            }
        }
        std::vector<NodeId> remap(nodes_.size(), no_node);
        std::vector<ComponentAtom> result;
        for (NodeId v = 0; v < nodes_.size(); ++v) {
            if (!nodes_[v].alive) continue;
            remap[v] = result.size();
            const auto& node = nodes_[v];
            result.push_back({node.type, node.leaves, node.active_leaf, {}});
        }
        for (NodeId v = 0; v < nodes_.size(); ++v) {
            if (!nodes_[v].alive) continue;
            for (NodeId c = nodes_[v].first; c != no_node; c = nodes_[c].next)
                result[remap[v]].children.push_back(remap[c]);
        }
        if (stats_) stats_->contract_output_nodes += result.size();
        return ComponentTree(std::move(result), remap[root_]);
    }

private:
    NodeId root_;
    CompressionStats* stats_;
    std::vector<MutableAtom> nodes_;
    std::deque<NodeId> queue_;

    bool zero_group(NodeId v) const {
        return nodes_[v].type == ComponentType::G && nodes_[v].active_leaf == no_node;
    }
    void enqueue(NodeId v) {
        if (v == no_node || !nodes_[v].alive || nodes_[v].queued) return;
        nodes_[v].queued = true;
        queue_.push_back(v);
    }
    void changed(NodeId v) {
        if (stats_) ++stats_->merges;
        enqueue(v);
        enqueue(nodes_[v].parent);
    }
    void attach(NodeId parent, NodeId child) {
        auto& p = nodes_[parent];
        auto& c = nodes_[child];
        c.parent = parent;
        c.previous = p.last;
        c.next = no_node;
        if (p.last == no_node) p.first = child;
        else nodes_[p.last].next = child;
        p.last = child;
        ++p.children;
        if (c.type == ComponentType::G) ++p.groups;
        if (zero_group(child)) {
            c.zero_previous = no_node;
            c.zero_next = p.first_zero;
            if (p.first_zero != no_node) nodes_[p.first_zero].zero_previous = child;
            p.first_zero = child;
            ++p.zero_groups;
        }
    }
    void detach(NodeId child) {
        auto& c = nodes_[child];
        auto& p = nodes_[c.parent];
        if (c.previous == no_node) p.first = c.next;
        else nodes_[c.previous].next = c.next;
        if (c.next == no_node) p.last = c.previous;
        else nodes_[c.next].previous = c.previous;
        --p.children;
        if (c.type == ComponentType::G) --p.groups;
        if (zero_group(child)) {
            if (c.zero_previous == no_node) p.first_zero = c.zero_next;
            else nodes_[c.zero_previous].zero_next = c.zero_next;
            if (c.zero_next != no_node) nodes_[c.zero_next].zero_previous = c.zero_previous;
            --p.zero_groups;
        }
        c.parent = c.previous = c.next = c.zero_previous = c.zero_next = no_node;
    }
    void normalize(NodeId v) {
        auto& node = nodes_[v];
        const NodeId parent = node.parent;
        if (node.type == ComponentType::C && node.children == 0 && parent != no_node &&
            nodes_[parent].type == ComponentType::I) {
            // Reattach so the parent's G / zero-G lists reflect the new type.
            detach(v);
            node.type = ComponentType::G;
            attach(parent, v);
            enqueue(parent);
        }
    }
};

}  // namespace

ComponentTree contract_components(const ComponentTree& tree, CompressionStats* stats) {
    // Size bound after greedy termination, for k > 0 designated leaves:
    // - There is at most one zero G per I (otherwise GG still applies).
    // - After deleting zero branches, an I either branches toward >=2 active
    //   leaves or is adjacent to its sole active G. Thus at most 2k-1 I atoms.
    // - A zero C must point to I or an active C (otherwise CC still applies).
    //   There are at most 3k-1 such C atoms, and at most 2k-1 zero G atoms.
    // Together with k designated atoms this is at most 8k-3 <= 8k.
    return Contractor(tree, stats).run();
}

}  // namespace soda2013
