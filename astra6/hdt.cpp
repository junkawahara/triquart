// The non-overlapping rounds in soda13.pdf, Section 4, p. 5, rules (i)–(iii).
// We rebuild only the CURRENT component tree after each round. Scanning the
// original tree in every round would lose the O(n) construction bound.
#include "hdt.hpp"

#include <utility>

namespace soda2013 {
namespace {

struct WorkingComponent {
    NodeId hdt;
    NodeId parent = no_node;
    std::vector<NodeId> children;
};

struct Pair {
    NodeId first, second;
    Composition composition;
};

}  // namespace

Hdt::Hdt(const Tree& tree) : Hdt(ComponentTree::from_tree(tree)) {}

Hdt::Hdt(const ComponentTree& tree) {
    nodes_.reserve(tree.size() * 2 - 1);
    std::vector<WorkingComponent> current;
    current.reserve(tree.size());
    vertex_to_leaf_.reserve(tree.size());
    for (NodeId v = 0; v < tree.size(); ++v) {
        HdtNode node;
        node.source = v;
        const auto& atom = tree.atom(v);
        node.type = atom.type;
        node.leaf_mass = atom.leaves;
        node.downwards_open = !atom.children.empty();
        node.composition = Composition::Atom;
        if (atom.type == ComponentType::I) node.composition = Composition::IdentityI;
        if (atom.type == ComponentType::C && atom.active_leaf != no_node) {
            node.composition = Composition::LeafToC;
            ++stats_.leaf_to_c;
        }
        nodes_.push_back(node);
        vertex_to_leaf_.push_back(v);
        current.push_back({v, tree.parent(v), atom.children});
    }

    stats_.round_sizes.push_back(current.size());
    while (current.size() > 1) {
        const std::size_t size = current.size();
        stats_.visited_components += size;
        const auto type = [&](NodeId v) { return nodes_[current[v].hdt].type; };

        // Section 4: convert all closed C children of I to G before pairing.
        for (NodeId v = 0; v < size; ++v) {
            const NodeId p = current[v].parent;
            if (type(v) == ComponentType::C && current[v].children.empty() &&
                p != no_node && type(p) == ComponentType::I) {
                auto& node = nodes_[current[v].hdt];
                node.type = ComponentType::G;
                node.c_to_g = true;
                ++stats_.c_to_g;
            }
        }

        std::vector<Pair> pairs;
        std::vector<bool> paired(size, false);
        const auto pair = [&](NodeId a, NodeId b, Composition composition) {
            if (paired[a] || paired[b]) throw std::logic_error("overlapping HDT compositions");
            paired[a] = paired[b] = true;
            pairs.push_back({a, b, composition});
        };

        for (NodeId v = 0; v < size; ++v) {
            if (type(v) != ComponentType::I) continue;
            std::vector<NodeId> groups;
            std::size_t open = 0;
            for (NodeId child : current[v].children) {
                if (type(child) == ComponentType::G) groups.push_back(child);
                else ++open;
            }
            if (groups.size() >= 2) {
                // (i): pair G siblings. I itself is not merged this round.
                for (std::size_t i = 0; i + 1 < groups.size(); i += 2)
                    pair(groups[i], groups[i + 1], Composition::GG);
            } else if (open <= 1 && groups.size() == 1) {
                // (ii): the I may also have one downwards-open child. "Only
                // one child being a G" does NOT mean "only one child in all".
                pair(v, groups[0], Composition::IG);
            }
        }
        for (NodeId v = 0; v < size; ++v) {
            const NodeId parent = current[v].parent;
            if (type(v) != ComponentType::C ||
                (parent != no_node && type(parent) == ComponentType::C)) continue;
            // (iii): pair from the TOP of every maximal consecutive C path.
            NodeId upper = v;
            while (!current[upper].children.empty()) {
                if (current[upper].children.size() != 1)
                    throw std::logic_error("C component has multiple downward edges");
                const NodeId lower = current[upper].children[0];
                if (type(lower) != ComponentType::C) break;
                pair(lower, upper, Composition::CC);  // C1 below C2, Fig. 6.
                if (current[lower].children.empty()) break;
                upper = current[lower].children[0];
                if (type(upper) != ComponentType::C) break;
            }
        }
        if (pairs.empty()) throw std::logic_error("HDT round made no progress");

        std::vector<NodeId> remap(size, no_node);
        std::vector<WorkingComponent> next;
        next.reserve(size - pairs.size());
        for (const auto& p : pairs) {
            const NodeId left = current[p.first].hdt;
            const NodeId right = current[p.second].hdt;
            HdtNode node;
            node.left = left;
            node.right = right;
            node.composition = p.composition;
            node.type = p.composition == Composition::GG ? ComponentType::G : ComponentType::C;
            node.downwards_open = p.composition == Composition::CC ? nodes_[left].downwards_open :
                                    p.composition == Composition::IG && current[p.first].children.size() > 1;
            node.mass = nodes_[left].mass + nodes_[right].mass;
            node.leaf_mass = add(nodes_[left].leaf_mass, nodes_[right].leaf_mass);
            node.height = 1 + std::max(nodes_[left].height, nodes_[right].height);
            nodes_[left].parent = nodes_[right].parent = nodes_.size();
            remap[p.first] = remap[p.second] = next.size();
            next.push_back({nodes_.size(), no_node, {}});
            nodes_.push_back(node);
            if (p.composition == Composition::CC) ++stats_.cc;
            if (p.composition == Composition::GG) ++stats_.gg;
            if (p.composition == Composition::IG) ++stats_.ig;
        }
        for (NodeId v = 0; v < size; ++v) {
            if (remap[v] != no_node) continue;
            remap[v] = next.size();
            next.push_back({current[v].hdt, no_node, {}});
        }
        for (NodeId v = 0; v < size; ++v) {
            const NodeId p = current[v].parent;
            if (p == no_node || remap[v] == remap[p]) continue;
            const NodeId child = remap[v], parent = remap[p];
            if (next[child].parent == no_node) {
                next[child].parent = parent;
                next[parent].children.push_back(child);
            } else if (next[child].parent != parent) {
                throw std::logic_error("component has two parents");
            }
            // Two G siblings merged above have the same new parent edge.
        }
        current = std::move(next);
        stats_.round_sizes.push_back(current.size());
    }
    root_ = current[0].hdt;
}

}  // namespace soda2013
