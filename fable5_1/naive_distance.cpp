// naive_distance.cpp — see naive_distance.hpp.
#include "naive_distance.hpp"

#include <algorithm>

namespace soda13 {

NaiveTopology::NaiveTopology(const Tree& t) {
    const int m = t.num_nodes();
    depth.assign(m, 0);
    parent.assign(m, -1);
    int nl = 0;
    for (int v = 0; v < m; ++v)
        if (t.is_leaf(v)) nl = std::max(nl, t.nodes[v].leaf_id + 1);
    leaf_node.assign(nl, -1);
    std::vector<int> st{t.root};
    while (!st.empty()) {
        int v = st.back();
        st.pop_back();
        parent[v] = t.nodes[v].parent;
        depth[v] = parent[v] < 0 ? 0 : depth[parent[v]] + 1;
        if (t.is_leaf(v)) leaf_node[t.nodes[v].leaf_id] = v;
        for (int c : t.nodes[v].children) st.push_back(c);
    }
}

int NaiveTopology::lca(int u, int v) const {
    while (u != v) {
        if (depth[u] < depth[v]) std::swap(u, v);
        u = parent[u];
    }
    return u;
}

int NaiveTopology::dist(int u, int v) const {
    int w = lca(u, v);
    return depth[u] + depth[v] - 2 * depth[w];
}

int NaiveTopology::triplet(int a, int b, int c) const {
    // Rooted topology by LCA depths (Sec. 2, Fig. 1): the pair whose LCA is
    // strictly deepest is the resolved pair; all equal => star.
    const int na = leaf_node[a], nb = leaf_node[b], nc = leaf_node[c];
    const int dab = depth[lca(na, nb)], dac = depth[lca(na, nc)], dbc = depth[lca(nb, nc)];
    if (dab > dac && dab > dbc) return 0;
    if (dac > dab && dac > dbc) return 1;
    if (dbc > dab && dbc > dac) return 2;
    return 3;
}

int NaiveTopology::quartet(int a, int b, int c, int d) const {
    // Four-point condition with unit edge lengths (exact for tree metrics):
    // ab|cd iff d(a,b)+d(c,d) < d(a,c)+d(b,d) (= d(a,d)+d(b,c)); star if all
    // three sums coincide (Sec. 2, Fig. 2).
    const int na = leaf_node[a], nb = leaf_node[b], nc = leaf_node[c], nd = leaf_node[d];
    const int s0 = dist(na, nb) + dist(nc, nd);
    const int s1 = dist(na, nc) + dist(nb, nd);
    const int s2 = dist(na, nd) + dist(nb, nc);
    const int mn = std::min({s0, s1, s2});
    int cnt = (s0 == mn) + (s1 == mn) + (s2 == mn);
    if (cnt != 1) return 3;
    if (s0 == mn) return 0;
    if (s1 == mn) return 1;
    return 2;
}

static void classify(CaseCounts& cc, int top1, int top2) {
    const bool r1 = top1 != 3, r2 = top2 != 3;
    if (r1 && r2) {
        if (top1 == top2) ++cc.A;
        else ++cc.B;
    } else if (r1 && !r2) {
        ++cc.C;
    } else if (!r1 && r2) {
        ++cc.D;
    } else {
        ++cc.E;
    }
}

CaseCounts naive_triplet_counts(const Tree& t1, const Tree& t2) {
    NaiveTopology p1(t1), p2(t2);
    const int n = static_cast<int>(p1.leaf_node.size());
    CaseCounts cc;
    for (int a = 0; a < n; ++a)
        for (int b = a + 1; b < n; ++b)
            for (int c = b + 1; c < n; ++c) classify(cc, p1.triplet(a, b, c), p2.triplet(a, b, c));
    return cc;
}

CaseCounts naive_quartet_counts(const Tree& t1, const Tree& t2) {
    NaiveTopology p1(t1), p2(t2);
    const int n = static_cast<int>(p1.leaf_node.size());
    CaseCounts cc;
    for (int a = 0; a < n; ++a)
        for (int b = a + 1; b < n; ++b)
            for (int c = b + 1; c < n; ++c)
                for (int d = c + 1; d < n; ++d)
                    classify(cc, p1.quartet(a, b, c, d), p2.quartet(a, b, c, d));
    return cc;
}

}  // namespace soda13
