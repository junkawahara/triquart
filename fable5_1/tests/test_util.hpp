// tests/test_util.hpp — random tree generation for oracle-based tests.
#pragma once

#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "../tree.hpp"

namespace soda13 {
namespace test {

// Random rooted tree on leaves "L0".."L{n-1}" with arbitrary degree.
// `max_deg` bounds the number of children; `binary_bias` in [0,1] is the
// probability of attaching a new leaf by subdividing an edge (binary-ish)
// instead of adding it as an extra child of an existing internal node
// (higher degree, i.e. unresolved topologies).
inline Tree random_tree(int n, std::mt19937& rng, double binary_bias = 0.5) {
    Tree t;
    if (n == 1) {
        t.nodes.emplace_back();
        t.nodes[0].label = "L0";
        t.root = 0;
        return t;
    }
    // Start from a cherry (L0,L1).
    t.nodes.resize(3);
    t.root = 0;
    t.nodes[0].children = {1, 2};
    t.nodes[1].parent = t.nodes[2].parent = 0;
    t.nodes[1].label = "L0";
    t.nodes[2].label = "L1";
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (int k = 2; k < n; ++k) {
        int leaf = static_cast<int>(t.nodes.size());
        t.nodes.emplace_back();
        t.nodes[leaf].label = "L" + std::to_string(k);
        std::uniform_int_distribution<int> pick(0, leaf - 1);
        int v = pick(rng);
        if (U(rng) < binary_bias || t.is_leaf(v)) {
            // subdivide the edge above v (or create a new root above it)
            int mid = static_cast<int>(t.nodes.size());
            t.nodes.emplace_back();
            int p = t.nodes[v].parent;
            t.nodes[mid].parent = p;
            if (p >= 0) {
                for (int& c : t.nodes[p].children)
                    if (c == v) c = mid;
            } else {
                t.root = mid;
            }
            t.nodes[v].parent = mid;
            t.nodes[mid].children = {v, leaf};
            t.nodes[leaf].parent = mid;
        } else {
            t.nodes[v].children.push_back(leaf);
            t.nodes[leaf].parent = v;
        }
    }
    return t;
}

// Shuffle leaf labels of a tree (same label set, different placement).
inline void shuffle_labels(Tree& t, std::mt19937& rng) {
    std::vector<int> lv = t.leaves();
    std::vector<std::string> labs;
    for (int v : lv) labs.push_back(t.nodes[v].label);
    std::shuffle(labs.begin(), labs.end(), rng);
    for (size_t k = 0; k < lv.size(); ++k) t.nodes[lv[k]].label = labs[k];
}

// Root with ~sqrt(n) children, each a random binary tree (max degree ~ sqrt(n)).
inline Tree star_of_binaries(int n, std::mt19937& rng) {
    int k = std::max(2, (int)std::sqrt((double)n));
    Tree t;
    t.nodes.emplace_back();
    t.root = 0;
    int next = 0;
    for (int g = 0; g < k; ++g) {
        int lo = (int)((long long)n * g / k), hi = (int)((long long)n * (g + 1) / k);
        int m = hi - lo;
        if (m <= 0) continue;
        Tree sub = random_tree(m, rng, 1.0);
        int base = (int)t.nodes.size();
        for (auto& nd : sub.nodes) {
            Tree::Node c = nd;
            c.parent = nd.parent < 0 ? 0 : nd.parent + base;
            for (int& ch : c.children) ch += base;
            if (nd.children.empty()) c.label = "L" + std::to_string(next++);
            t.nodes.push_back(c);
        }
        t.nodes[0].children.push_back(sub.root + base);
    }
    return t;
}

// Caterpillar ((((L0,L1),L2),L3),...).
inline Tree caterpillar(int n) {
    Tree t;
    t.nodes.resize(2);
    t.nodes[0].children = {1};
    t.nodes[1].parent = 0;
    t.nodes[1].label = "L0";
    t.root = 0;
    int spine = 0;
    for (int k = 1; k < n; ++k) {
        int leaf = (int)t.nodes.size();
        t.nodes.emplace_back();
        t.nodes[leaf].label = "L" + std::to_string(k);
        t.nodes[leaf].parent = spine;
        t.nodes[spine].children.push_back(leaf);
        if (k + 1 < n) {
            int ns = (int)t.nodes.size();
            t.nodes.emplace_back();
            t.nodes[ns].parent = -1;
            t.nodes[ns].children = {spine};
            t.nodes[spine].parent = ns;
            spine = ns;
            t.root = ns;
        }
    }
    return t;
}

}  // namespace test
}  // namespace soda13
