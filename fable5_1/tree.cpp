// tree.cpp — see tree.hpp.
#include "tree.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <functional>
#include <set>

namespace soda13 {

int Tree::num_leaves() const {
    int c = 0;
    for (const auto& nd : nodes) c += nd.children.empty() ? 1 : 0;
    return c;
}

std::vector<int> Tree::leaves() const {
    std::vector<int> out;
    std::vector<int> st;
    if (root >= 0) st.push_back(root);
    while (!st.empty()) {
        int v = st.back();
        st.pop_back();
        if (nodes[v].children.empty()) out.push_back(v);
        for (auto it = nodes[v].children.rbegin(); it != nodes[v].children.rend(); ++it) st.push_back(*it);
    }
    return out;
}

Tree Tree::from_newick(const std::string& s) {
    Tree t;
    size_t pos = 0;
    auto skip_ws = [&]() {
        while (pos < s.size() && std::isspace(static_cast<unsigned char>(s[pos]))) ++pos;
    };
    auto read_label = [&]() {
        std::string lab;
        skip_ws();
        if (pos < s.size() && (s[pos] == '\'' || s[pos] == '"')) {
            char q = s[pos++];
            while (pos < s.size() && s[pos] != q) lab += s[pos++];
            if (pos < s.size()) ++pos;
        } else {
            while (pos < s.size() && !std::strchr("(),:;", s[pos]) &&
                   !std::isspace(static_cast<unsigned char>(s[pos])))
                lab += s[pos++];
        }
        skip_ws();
        if (pos < s.size() && s[pos] == ':') {  // branch length, ignored
            ++pos;
            while (pos < s.size() && !std::strchr("(),;", s[pos])) ++pos;
        }
        return lab;
    };
    // Iterative descent: `open` holds the ids of the currently open internal nodes.
    std::vector<int> open;
    auto new_node = [&]() {
        int id = static_cast<int>(t.nodes.size());
        t.nodes.emplace_back();
        if (!open.empty()) {
            t.nodes[id].parent = open.back();
            t.nodes[open.back()].children.push_back(id);
        }
        return id;
    };
    skip_ws();
    while (true) {
        skip_ws();
        if (pos < s.size() && s[pos] == '(') {
            ++pos;
            open.push_back(new_node());
            continue;
        }
        if (pos < s.size() && s[pos] == ')') {
            if (open.empty()) throw std::invalid_argument("Newick: unbalanced ')'");
            ++pos;
            int id = open.back();
            open.pop_back();
            if (t.nodes[id].children.empty()) throw std::invalid_argument("Newick: empty group");
            read_label();  // internal label ignored
            if (open.empty()) {
                t.root = id;
                break;
            }
            skip_ws();
            if (pos < s.size() && s[pos] == ',') {
                ++pos;
                continue;
            }
            if (pos < s.size() && s[pos] == ')') continue;
            throw std::invalid_argument("Newick: expected ',' or ')'");
        }
        // leaf
        int id = new_node();
        std::string lab = read_label();
        if (lab.empty()) throw std::invalid_argument("Newick: leaf without label");
        t.nodes[id].label = lab;
        if (open.empty()) {
            t.root = id;
            break;
        }
        skip_ws();
        if (pos < s.size() && s[pos] == ',') {
            ++pos;
            continue;
        }
        if (pos < s.size() && s[pos] == ')') continue;
        throw std::invalid_argument("Newick: expected ',' or ')'");
    }
    skip_ws();
    if (pos < s.size() && s[pos] == ';') ++pos;
    skip_ws();
    if (pos != s.size()) throw std::invalid_argument("Newick: trailing characters");
    return t;
}

std::string Tree::to_newick() const {
    std::string out;
    // explicit stack: (node, next child index)
    std::vector<std::pair<int, size_t>> st{{root, 0}};
    while (!st.empty()) {
        auto& [v, k] = st.back();
        const Node& nd = nodes[v];
        if (nd.children.empty()) {
            out += nd.label;
            st.pop_back();
            continue;
        }
        if (k == 0) out += '(';
        if (k < nd.children.size()) {
            if (k) out += ',';
            int c = nd.children[k++];
            st.emplace_back(c, 0);
        } else {
            out += ')';
            st.pop_back();
        }
    }
    out += ';';
    return out;
}

void Tree::contract_unary() {
    // Rebuild the node array keeping only leaves and internal nodes with >= 2
    // children (iterative to survive very deep trees).
    Tree t;
    if (root < 0) return;
    auto skip = [&](int v) {
        while (nodes[v].children.size() == 1) v = nodes[v].children[0];
        return v;
    };
    struct Item {
        int old, parent_new;
    };
    std::vector<Item> st{{skip(root), -1}};
    while (!st.empty()) {
        Item it = st.back();
        st.pop_back();
        int id = static_cast<int>(t.nodes.size());
        t.nodes.emplace_back();
        t.nodes[id].label = nodes[it.old].label;
        t.nodes[id].leaf_id = nodes[it.old].leaf_id;
        t.nodes[id].parent = it.parent_new;
        if (it.parent_new >= 0) t.nodes[it.parent_new].children.push_back(id);
        else t.root = id;
        const auto& ch = nodes[it.old].children;
        for (auto c = ch.rbegin(); c != ch.rend(); ++c) st.push_back({skip(*c), id});
    }
    *this = std::move(t);
}

std::vector<int> Tree::subtree_leaf_counts() const {
    std::vector<int> cnt(nodes.size(), 0);
    // post-order via explicit stack
    std::vector<int> order;
    std::vector<int> st{root};
    while (!st.empty()) {
        int v = st.back();
        st.pop_back();
        order.push_back(v);
        for (int c : nodes[v].children) st.push_back(c);
    }
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        int v = *it;
        if (nodes[v].children.empty()) {
            cnt[v] = 1;
        } else {
            for (int c : nodes[v].children) cnt[v] += cnt[c];
        }
    }
    return cnt;
}

int Tree::degree_max() const {
    int d = 0;
    for (const auto& nd : nodes) d = std::max<int>(d, static_cast<int>(nd.children.size()));
    return d;
}

void assign_leaf_ids(Tree& t1, Tree& t2) {
    std::map<std::string, int> ids;
    for (auto& nd : t1.nodes) {
        if (!nd.children.empty()) continue;
        if (nd.label.empty()) throw std::invalid_argument("leaf without label in T1");
        if (ids.count(nd.label)) throw std::invalid_argument("duplicate leaf label in T1: " + nd.label);
        int id = static_cast<int>(ids.size());
        ids[nd.label] = id;
        nd.leaf_id = id;
    }
    if (ids.empty()) throw std::invalid_argument("T1 has no leaves");
    std::set<int> seen;
    for (auto& nd : t2.nodes) {
        if (!nd.children.empty()) continue;
        auto it = ids.find(nd.label);
        if (it == ids.end()) throw std::invalid_argument("leaf label only in T2: " + nd.label);
        if (!seen.insert(it->second).second) throw std::invalid_argument("duplicate leaf label in T2: " + nd.label);
        nd.leaf_id = it->second;
    }
    if (seen.size() != ids.size()) throw std::invalid_argument("leaf label sets differ");
}

}  // namespace soda13
