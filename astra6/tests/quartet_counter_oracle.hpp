#pragma once

// Semantic interpreter of Section 7 / Figure 11 subscripts, for TESTS ONLY.
// It enumerates distinct leaves, checks colors, LCAs and external-path branch
// points. It neither knows nor uses any Figure 12 composition equation.
#include "quartet_shared_counters.hpp"
#include "test_support.hpp"

#include <array>
#include <set>
#include <sstream>

namespace counter_oracle {
using namespace soda2013;

inline NodeId lca(const Tree& tree, NodeId a, NodeId b) {
    while (tree.depth(a) > tree.depth(b)) a = tree.parent(a);
    while (tree.depth(b) > tree.depth(a)) b = tree.parent(b);
    while (a != b) { a = tree.parent(a); b = tree.parent(b); }
    return a;
}

struct Geometry {
    const Tree& tree;
    ComponentType type;
    std::vector<NodeId> leaves;
    std::vector<bool> path;
    std::vector<NodeId> branch;
    NodeId super_root = no_node;
    const std::vector<NodeId>& lcas;
    NodeId ancestor(NodeId a, NodeId b) const { return lcas[a * tree.size() + b]; }
};

class ComponentGeometry {
public:
    ComponentGeometry(const Tree& tree, const Hdt& hdt) : tree_(tree), lcas_(tree.size() * tree.size()) {
        for (NodeId a = 0; a < tree.size(); ++a)
            for (NodeId b = 0; b < tree.size(); ++b) lcas_[a * tree.size() + b] = lca(tree, a, b);
        std::vector<std::set<NodeId>> vertices(hdt.nodes().size()), paths(hdt.nodes().size());
        for (NodeId uid = 0; uid < hdt.nodes().size(); ++uid) {
            const auto& node = hdt.node(uid);
            if (node.left == no_node) {
                vertices[uid].insert(node.source);
                paths[uid].insert(node.source);
            } else {
                for (NodeId child : {node.left, node.right})
                    vertices[uid].insert(vertices[child].begin(), vertices[child].end());
                if (node.composition == Composition::CC) {
                    paths[uid] = paths[node.left];
                    paths[uid].insert(paths[node.right].begin(), paths[node.right].end());
                } else if (node.composition == Composition::IG) paths[uid] = vertices[node.left];
            }
            Geometry geometry{tree_, node.type, {}, std::vector<bool>(tree.size()),
                              std::vector<NodeId>(tree.size(), no_node), no_node, lcas_};
            for (NodeId v : paths[uid]) geometry.path[v] = true;
            for (NodeId v : vertices[uid]) {
                if (tree.is_leaf(v)) geometry.leaves.push_back(v);
                if (node.type == ComponentType::G && !vertices[uid].count(tree.parent(v)))
                    geometry.super_root = tree.parent(v);
            }
            if (node.type == ComponentType::C)
                for (NodeId leaf : geometry.leaves) {
                    NodeId v = leaf;
                    while (!geometry.path[v]) { v = tree.parent(v); CHECK(v != no_node); }
                    geometry.branch[leaf] = v;
                }
            geometries_.push_back(std::move(geometry));
        }
    }
    const Geometry& operator[](NodeId uid) const { return geometries_.at(uid); }
private:
    const Tree& tree_;
    std::vector<NodeId> lcas_;
    std::vector<Geometry> geometries_;
};

enum class Placement { at, inside, any };
struct Term {
    int leaf = -1;
    std::vector<std::size_t> children;
};
struct Segment {
    Placement placement;
    std::size_t root, begin, end;
};

class Pattern {
public:
    explicit Pattern(std::string_view expression) : expression_(expression) {
        has_i_ = expression.find('i') != std::string_view::npos;
        has_j_ = expression.find('j') != std::string_view::npos;
        std::size_t start = 0;
        int depth = 0;
        for (std::size_t pos = 0; pos <= expression.size(); ++pos) {
            if (pos == expression.size() || (expression[pos] == '^' && depth == 0)) {
                auto body = expression.substr(start, pos - start);
                Placement placement = Placement::at;
                if (body.front() == '[') { placement = Placement::any; body = body.substr(1, body.size() - 2); }
                else if (body.front() == '(') { placement = Placement::inside; body = body.substr(1, body.size() - 2); }
                const std::size_t begin = roles_.size();
                std::size_t p = 0;
                const auto root = parse(body, p);
                CHECK(p == body.size());
                segments_.push_back({placement, root, begin, roles_.size()});
                start = pos + 1;
            } else if (expression[pos] == '(' || expression[pos] == '[') ++depth;
            else if (expression[pos] == ')' || expression[pos] == ']') --depth;
        }
        CHECK(roles_.size() >= 1 && roles_.size() <= 3 && terms_.size() <= 8);
        // Leaf permutations within an unordered sibling pair are duplicates
        // when the colors coincide, or when both colors are anonymous b/c.
        for (const auto& term : terms_) {
            if (term.children.size() != 2) continue;
            const auto& a = terms_[term.children[0]];
            const auto& b = terms_[term.children[1]];
            if (a.leaf < 0 || b.leaf < 0) continue;
            const char x = roles_[static_cast<std::size_t>(a.leaf)], y = roles_[static_cast<std::size_t>(b.leaf)];
            const bool anonymous_pair = (x == 'b' && y == 'c') || (x == 'c' && y == 'b');
            // If b also occurs outside this pair (e.g. b(bc)), its color is
            // fixed by that occurrence and b/c can no longer be exchanged.
            const bool anonymous_exchange = anonymous_pair &&
                std::count(roles_.begin(), roles_.end(), 'b') == 1 &&
                std::count(roles_.begin(), roles_.end(), 'c') == 1;
            if (x == y || anonymous_exchange) normalizer_ *= 2;
        }
    }

    Count count(const Geometry& geometry, const std::vector<Color>& colors, Color index = 0, Color second_index = 0) const {
        if (segments_.size() > 1 && geometry.type != ComponentType::C) return 0;  // C-only counters.
        Count result = 0;
        std::array<NodeId, 3> selected{};
        for (NodeId a : geometry.leaves) {
            selected[0] = a;
            if (roles_.size() == 1) { result += matches(geometry, colors, index, second_index, selected); continue; }
            for (NodeId b : geometry.leaves) {
                if (a == b) continue;
                selected[1] = b;
                if (roles_.size() == 2) { result += matches(geometry, colors, index, second_index, selected); continue; }
                for (NodeId c : geometry.leaves) {
                    if (c == a || c == b) continue;
                    selected[2] = c;
                    result += matches(geometry, colors, index, second_index, selected);
                }
            }
        }
        if (result % normalizer_ != 0)
            throw std::runtime_error("oracle symmetry mismatch for " + std::string(expression_) +
                                     ": assignments=" + std::to_string(result) +
                                     " divisor=" + std::to_string(normalizer_));
        return result / normalizer_;
    }
    std::string_view expression() const { return expression_; }

private:
    std::string_view expression_;
    std::vector<char> roles_;
    std::vector<Term> terms_;
    std::vector<Segment> segments_;
    bool has_i_ = false;
    bool has_j_ = false;
    Count normalizer_ = 1;

    std::size_t parse(std::string_view body, std::size_t& pos) {
        std::vector<std::size_t> children;
        while (pos < body.size() && body[pos] != ')') {
            if (body[pos] == '(') {
                ++pos;
                children.push_back(parse(body, pos));
                CHECK(pos < body.size() && body[pos] == ')');
                ++pos;
            } else {
                CHECK(body[pos] == '0' || body[pos] == 'i' || body[pos] == 'j' || body[pos] == 'b' || body[pos] == 'c');
                children.push_back(terms_.size());
                terms_.push_back({static_cast<int>(roles_.size()), {}});
                roles_.push_back(body[pos++]);
            }
        }
        CHECK(!children.empty());
        if (children.size() == 1) return children[0];
        const auto root = terms_.size();
        terms_.push_back({-1, std::move(children)});
        return root;
    }

    bool matches(const Geometry& g, const std::vector<Color>& colors, Color index, Color second_index,
                 const std::array<NodeId, 3>& leaves) const {
        Color b = 0, c = 0;
        for (std::size_t i = 0; i < roles_.size(); ++i) {
            const Color color = colors[leaves[i]];
            if (roles_[i] == '0') { if (color != 0) return false; }
            else if (roles_[i] == 'i') { if (color != index) return false; }
            else if (roles_[i] == 'j') { if (color != second_index) return false; }
            else {
                if (color == 0 || (has_i_ && color == index) || (has_j_ && color == second_index)) return false;
                Color& variable = roles_[i] == 'b' ? b : c;
                if (variable && color != variable) return false;
                variable = color;
                if (b && b == c) return false;
            }
        }
        std::array<NodeId, 8> ancestors{};
        for (std::size_t i = 0; i < terms_.size(); ++i) {
            const auto& term = terms_[i];
            if (term.leaf >= 0) { ancestors[i] = leaves[static_cast<std::size_t>(term.leaf)]; continue; }
            ancestors[i] = ancestors[term.children[0]];
            for (std::size_t child : term.children) ancestors[i] = g.ancestor(ancestors[i], ancestors[child]);
            // Parentheses inside a topology require a strictly lower LCA.
            for (std::size_t child : term.children)
                if (terms_[child].leaf < 0 && g.tree.depth(ancestors[child]) <= g.tree.depth(ancestors[i])) return false;
        }
        std::size_t previous_depth = g.tree.size();
        for (const auto& segment : segments_) {
            const NodeId ancestor = ancestors[segment.root];
            if (segment.end - segment.begin > 1 && segment.placement != Placement::any) {
                const bool on_boundary = g.type == ComponentType::G ? ancestor == g.super_root : g.path[ancestor];
                if (segment.placement == Placement::inside) { if (on_boundary) return false; }
                else {
                    if (!on_boundary) return false;
                    if (g.type == ComponentType::C)
                        for (std::size_t i = segment.begin; i < segment.end; ++i)
                            if (g.branch[leaves[i]] != ancestor) return false;
                }
            }
            if (segments_.size() > 1) {
                const NodeId branch = g.branch[leaves[segment.begin]];
                for (std::size_t i = segment.begin; i < segment.end; ++i)
                    if (g.branch[leaves[i]] != branch) return false;
                const auto depth = g.tree.depth(branch);
                if (depth >= previous_depth) return false;
                previous_depth = depth;
            }
        }
        return true;
    }
};

inline std::string diagnostic(const Tree& tree, const std::vector<Color>& colors, NodeId uid,
                              std::string_view subscript, Color index, Count actual, Count expected) {
    std::vector<std::string> subtrees(tree.size());
    for (auto it = tree.preorder().rbegin(); it != tree.preorder().rend(); ++it) {
        if (tree.is_leaf(*it)) subtrees[*it] = tree.node(*it).label + ":" + std::to_string(colors[*it]);
        else {
            auto& text = subtrees[*it];
            text = "(";
            for (NodeId child : tree.node(*it).children) text += subtrees[child] + ",";
            text.back() = ')';
        }
    }
    std::ostringstream text;
    text << "n_" << subscript << " i=" << index << " uid=" << uid
         << " actual=" << actual << " expected=" << expected << " tree(colors)=" << subtrees[tree.root()];
    return text.str();
}

}  // namespace counter_oracle
