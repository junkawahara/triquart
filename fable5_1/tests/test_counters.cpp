// tests/test_counters.cpp — brute-force oracle for the Sec. 7 counters and
// the Fig. 13 / Fig. 16 contributions.
//
// Part 1: for a random T2 and a random colouring, every HDT node's counters
// (as composed by counters.hpp) are compared with a direct enumeration over
// the leaves of the node's component, using the definitions of Sec. 7 /
// Fig. 11 / Fig. 14 (see "Notation" in paper_figures.md).
//
// Part 2: the root contribution (A, B) is compared, case by case
// (αα, βα, ββ, γα, γβ, γγ), with a direct enumeration of the coloured
// quartets of T2 using the T1-type implied by the colours (Sec. 6) and the
// T2-type implied by the rooted T2 (Fig. 4).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <set>
#include <string>

#include "../hdt_deco.hpp"
#include "../naive_distance.hpp"
#include "test_util.hpp"

using namespace soda13;

namespace {

struct CompInfo {
    bool isC = true;
    std::vector<int> leaves;        // CTree leaf nodes in the component
    std::set<int> members;          // CTree nodes in the component
    std::vector<int> path;          // C: external path top..bottom
    int super_root = -1;            // G
};

struct Ctx {
    const CTree* t;
    const std::vector<int>* colour;
    std::vector<int> depth, parent;
    int lca(int u, int v) const {
        while (u != v) {
            if (depth[u] < depth[v]) std::swap(u, v);
            u = parent[u];
        }
        return u;
    }
    int col(int leaf) const {
        const auto& nd = t->nodes[leaf];
        return nd.leaf_id >= 0 ? (*colour)[nd.leaf_id] : 0;
    }
};

// Structural description of a leaf relative to the component.
struct LeafPos {
    int att;  // path node (C) / super root (G) where the leaf branches off
    int sc;   // side child through which it hangs (the leaf itself if att==leaf)
    int attdepth;
};

CompInfo component(const HdtDeco<QuartetFields>& d, int h, const Ctx& cx) {
    const Hdt& H = d.hdt();
    CompInfo ci;
    const HNode& n = H.nodes[h];
    ci.isC = n.is_C();
    std::vector<int> st{h};
    while (!st.empty()) {
        int x = st.back();
        st.pop_back();
        const HNode& m = H.nodes[x];
        if (m.op == HOp::LEAF) {
            ci.leaves.push_back(m.tnode);
            ci.members.insert(m.tnode);
        } else if (m.op == HOp::IG2C) {
            ci.members.insert(m.tnode);
        }
        if (m.c1 >= 0) st.push_back(m.c1);
        if (m.c2 >= 0) st.push_back(m.c2);
    }
    if (ci.isC) {
        int cur = n.bottom;
        std::vector<int> p;
        while (true) {
            p.push_back(cur);
            if (cur == n.top) break;
            cur = cx.parent[cur];
        }
        ci.path.assign(p.rbegin(), p.rend());
    } else {
        for (int v : ci.members)
            if (!ci.members.count(cx.parent[v])) {
                ci.super_root = cx.parent[v];
                break;
            }
    }
    return ci;
}

LeafPos leaf_pos(const CompInfo& ci, const Ctx& cx, int leaf) {
    LeafPos lp;
    if (ci.isC) {
        std::set<int> onpath(ci.path.begin(), ci.path.end());
        int cur = leaf, prev = leaf;
        while (!onpath.count(cur)) {
            prev = cur;
            cur = cx.parent[cur];
        }
        lp.att = cur;
        lp.sc = prev;
    } else {
        int cur = leaf, prev = leaf;
        while (cur != ci.super_root) {
            prev = cur;
            cur = cx.parent[cur];
        }
        lp.att = cur;
        lp.sc = prev;
    }
    lp.attdepth = cx.depth[lp.att];
    return lp;
}

// Token matching: '0' zero, 'i' colour i, 'j' colour j (two-colour counters),
// 'b' = •, 's' = □.
struct Pattern {
    const char* name;
    std::string kind;    // "pairNode","pairOff","node3","paren3","br3","up2","upNode","upPair","pairUp","up3","any2","br2"
    std::string tokens;  // one char per leaf role in order (x,y,z) / (x,y)
};

bool colours_match(const std::string& tok, const std::vector<int>& cols, int i, int j) {
    int bcol = -1;
    for (size_t k = 0; k < tok.size(); ++k) {
        const int c = cols[k];
        switch (tok[k]) {
            case '0':
                if (c != 0) return false;
                break;
            case 'i':
                if (c != i) return false;
                break;
            case 'j':
                if (c != j) return false;
                break;
            case 'b':
                if (c == 0) return false;
                if (tok.find('i') != std::string::npos && c == i) return false;
                if (bcol < 0) bcol = c;
                else if (bcol != c) return false;
                break;
            case 's':
                if (c == 0) return false;
                if (tok.find('i') != std::string::npos && c == i) return false;
                break;
        }
    }
    // □ must differ from •
    for (size_t k = 0; k < tok.size(); ++k)
        if (tok[k] == 's' && bcol >= 0 && cols[k] == bcol) return false;
    // two □ tokens in the same pattern denote the same colour (n_{•↑(□□)})
    {
        int scol = -1;
        for (size_t k = 0; k < tok.size(); ++k)
            if (tok[k] == 's') {
                if (scol < 0) scol = cols[k];
                else if (scol != cols[k]) return false;
            }
    }
    return true;
}

// Brute force a counter over the component.
i64 brute(const Pattern& P, const CompInfo& ci, const Ctx& cx, const std::vector<LeafPos>& lp, int i, int j) {
    const int L = static_cast<int>(ci.leaves.size());
    auto pairOff = [&](int a, int b) { return lp[a].att == lp[b].att && lp[a].sc == lp[b].sc; };
    auto pairNode = [&](int a, int b) { return lp[a].att == lp[b].att && lp[a].sc != lp[b].sc; };
    auto lcad = [&](int a, int b) { return cx.depth[cx.lca(ci.leaves[a], ci.leaves[b])]; };
    std::set<std::vector<int>> counted;
    i64 total = 0;
    const bool two = P.tokens.size() == 2;
    std::vector<int> idx(two ? 2 : 3);
    auto consider = [&](const std::vector<int>& perm) {
        std::vector<int> cols;
        for (int x : perm) cols.push_back(cx.col(ci.leaves[x]));
        if (!colours_match(P.tokens, cols, i, j)) return;
        bool ok = false;
        const int x = perm[0], y = perm[1], z = two ? -1 : perm[2];
        if (P.kind == "pairNode") ok = pairNode(x, y);
        else if (P.kind == "pairOff") ok = pairOff(x, y);
        else if (P.kind == "any2") ok = true;
        else if (P.kind == "up2") ok = lp[x].att != lp[y].att && lp[x].attdepth > lp[y].attdepth;
        else if (P.kind == "node3") ok = pairOff(y, z) && lp[x].att == lp[y].att && lp[x].sc != lp[y].sc;
        else if (P.kind == "paren3")
            ok = pairOff(x, y) && pairOff(y, z) && lcad(y, z) > lcad(x, y) && lcad(y, z) > lcad(x, z);
        else if (P.kind == "br3") ok = lcad(y, z) > lcad(x, y) && lcad(y, z) > lcad(x, z);
        else if (P.kind == "upPair") ok = pairOff(y, z) && lp[x].att != lp[y].att && lp[x].attdepth > lp[y].attdepth;
        else if (P.kind == "pairUp") ok = pairOff(y, z) && lp[x].att != lp[y].att && lp[y].attdepth > lp[x].attdepth;
        else if (P.kind == "upNode") ok = pairNode(y, z) && lp[x].att != lp[y].att && lp[x].attdepth > lp[y].attdepth;
        else if (P.kind == "up3")
            ok = lp[x].att != lp[y].att && lp[y].att != lp[z].att && lp[x].att != lp[z].att &&
                 lp[x].attdepth > lp[y].attdepth && lp[y].attdepth > lp[z].attdepth;
        if (!ok) return;
        std::vector<int> key(perm);
        std::sort(key.begin(), key.end());
        if (counted.insert(key).second) ++total;
    };
    if (two) {
        for (int a = 0; a < L; ++a)
            for (int b = 0; b < L; ++b)
                if (a != b) consider({a, b});
    } else {
        for (int a = 0; a < L; ++a)
            for (int b = 0; b < L; ++b)
                for (int c = 0; c < L; ++c)
                    if (a != b && b != c && a != c) consider({a, b, c});
    }
    return total;
}

// per-colour counters: name, kind, tokens, accessor
struct PerColour {
    Pattern p;
    std::function<i64(const QuartetFields&)> get;
};
#define PC(field, kind, tok) PerColour{{#field, kind, tok}, [](const QuartetFields& f) { return f.field; }}

std::vector<PerColour> per_colour_counters() {
    return {
        PC(n_0i, "pairNode", "0i"), PC(n_ii, "pairNode", "ii"), PC(n_ib, "pairNode", "ib"),
        PC(p_0i, "pairOff", "0i"), PC(p_ii, "pairOff", "ii"), PC(p_ib, "pairOff", "ib"),
        PC(n_0_pii, "node3", "0ii"), PC(n_b_pii, "node3", "bii"), PC(n_i_p0b, "node3", "i0b"),
        PC(n_i_pbb, "node3", "ibb"), PC(n_i_pbs, "node3", "ibs"), PC(n_i_p0i, "node3", "i0i"),
        PC(n_b_p0i, "node3", "b0i"), PC(n_0_pib, "node3", "0ib"), PC(n_b_pis, "node3", "bis"),
        PC(n_i_pib, "node3", "iib"), PC(n_b_pib, "node3", "bib"),
        PC(p_0_pii, "paren3", "0ii"), PC(p_b_pii, "paren3", "bii"), PC(p_i_p0b, "paren3", "i0b"),
        PC(p_i_pib, "paren3", "iib"), PC(p_i_p0i, "paren3", "i0i"), PC(p_0_pib, "paren3", "0ib"),
        PC(p_b_p0i, "paren3", "b0i"),
        PC(br_0_pii, "br3", "0ii"), PC(br_b_pii, "br3", "bii"), PC(br_i_p0b, "br3", "i0b"),
        PC(br_i_pib, "br3", "iib"), PC(br_0_pib, "br3", "0ib"), PC(br_i_p0i, "br3", "i0i"),
        PC(br_b_p0i, "br3", "b0i"),
        PC(n_0Ui, "up2", "0i"), PC(n_iU0, "up2", "i0"), PC(n_iUi, "up2", "ii"), PC(n_iUb, "up2", "ib"),
        PC(n_bUi, "up2", "bi"),
        PC(n_iU0b, "upNode", "i0b"), PC(n_bUii, "upNode", "bii"), PC(n_iUbb, "upNode", "ibb"),
        PC(n_iUbs, "upNode", "ibs"), PC(n_0Uib, "upNode", "0ib"), PC(n_bU0i, "upNode", "b0i"),
        PC(n_iUib, "upNode", "iib"), PC(n_0Uii, "upNode", "0ii"), PC(n_iU0i, "upNode", "i0i"),
        PC(n_bUis, "upNode", "bis"), PC(n_bUib, "upNode", "bib"),
        PC(n_0Upii, "upPair", "0ii"), PC(n_iUp0b, "upPair", "i0b"), PC(n_iUpbs, "upPair", "ibs"),
        PC(n_bUpii, "upPair", "bii"), PC(n_iUpbb, "upPair", "ibb"), PC(n_0Upib, "upPair", "0ib"),
        PC(n_iUpib, "upPair", "iib"), PC(n_bUp0i, "upPair", "b0i"), PC(n_iUp0i, "upPair", "i0i"),
        PC(n_bUpis, "upPair", "bis"), PC(n_bUpib, "upPair", "bib"),
        PC(n_piiU0, "pairUp", "0ii"), PC(n_piiUb, "pairUp", "bii"), PC(n_pbbUi, "pairUp", "ibb"),
        PC(n_p0bUi, "pairUp", "i0b"), PC(n_pibU0, "pairUp", "0ib"), PC(n_pibUi, "pairUp", "iib"),
        PC(n_p0iUi, "pairUp", "i0i"), PC(n_p0iUb, "pairUp", "b0i"),
        PC(n_0UiUi, "up3", "0ii"), PC(n_bUiUi, "up3", "bii"), PC(n_iUbU0, "up3", "ib0"),
        PC(n_iU0Ub, "up3", "i0b"), PC(n_bUiU0, "up3", "bi0"), PC(n_iUiU0, "up3", "ii0"),
        PC(n_0UbUi, "up3", "0bi"), PC(n_bU0Ui, "up3", "b0i"), PC(n_iU0Ui, "up3", "i0i"),
        PC(n_0UiUb, "up3", "0ib"), PC(n_iUiUb, "up3", "iib"), PC(n_iUbUi, "up3", "ibi"),
    };
}

struct Total {
    Pattern p;
    std::function<i64(const QuartetCounters::Totals&)> get;
};
#define TT(field, kind, tok) Total{{#field, kind, tok}, [](const QuartetCounters::Totals& t) { return t.field; }}

std::vector<Total> total_counters() {
    return {
        TT(n_0b, "pairNode", "0b"), TT(n_bb, "pairNode", "bb"), TT(n_bs, "pairNode", "bs"),
        TT(p_0b, "pairOff", "0b"), TT(p_bb, "pairOff", "bb"), TT(p_bs, "pairOff", "bs"),
        TT(n_s_pbb, "node3", "sbb"), TT(br_bs, "any2", "bs"), TT(n_0Ub, "up2", "0b"), TT(n_bU0, "up2", "b0"),
        TT(n_0Ubb, "upNode", "0bb"), TT(n_bUpss, "upPair", "bss"), TT(n_b_pbs, "node3", "bbs"),
        TT(n_bUbs, "upNode", "bbs"), TT(n_pbsUb, "pairUp", "bbs"), TT(n_bUbUs, "up3", "bbs"),
        TT(n_bUsUb, "up3", "bsb"), TT(p_b_pbs, "paren3", "bbs"), TT(br_b_pbs, "br3", "bbs"),
        TT(n_sUbb, "upNode", "sbb"), TT(n_bUpbs, "upPair", "bbs"),
    };
}

int check_counters(const CTree& tree, const std::vector<int>& colour, const HdtDeco<QuartetFields>& d, int maxreport) {
    Ctx cx;
    cx.t = &tree;
    cx.colour = &colour;
    const int m = tree.num_nodes();
    cx.depth.assign(m, 0);
    cx.parent.assign(m, -1);
    std::vector<int> st{tree.root};
    while (!st.empty()) {
        int v = st.back();
        st.pop_back();
        cx.parent[v] = tree.nodes[v].parent;
        cx.depth[v] = cx.parent[v] < 0 ? 0 : cx.depth[cx.parent[v]] + 1;
        for (int c : tree.nodes[v].children) st.push_back(c);
    }
    static const auto PCS = per_colour_counters();
    static const auto TTS = total_counters();
    int bad = 0;
    const Hdt& H = d.hdt();
    for (int h = 0; h < static_cast<int>(H.nodes.size()); ++h) {
        if (H.nodes[h].op == HOp::LEAF) continue;
        CompInfo ci = component(d, h, cx);
        std::vector<LeafPos> lp;
        for (int leaf : ci.leaves) lp.push_back(leaf_pos(ci, cx, leaf));
        const QuartetCounters& C = d.counters(h);
        std::set<int> cols;
        for (int leaf : ci.leaves)
            if (cx.col(leaf)) cols.insert(cx.col(leaf));
        auto report = [&](const char* name, int i, int j, i64 got, i64 exp) {
            ++bad;
            if (bad <= maxreport)
                std::printf("  hnode %d op %d (%s) counter %s i=%d j=%d: got %lld expected %lld\n", h, (int)H.nodes[h].op,
                            ci.isC ? "C" : "G", name, i, j, (long long)got, (long long)exp);
        };
        // colour list
        std::vector<int> expect_cols(cols.begin(), cols.end());
        if (expect_cols != C.colours) {
            report("colours", -1, -1, (i64)C.colours.size(), (i64)expect_cols.size());
            continue;
        }
        // n0 / nb
        i64 n0 = 0;
        for (int leaf : ci.leaves)
            if (cx.col(leaf) == 0) n0 += tree.nodes[leaf].weight;
        if (C.t.n0 != n0) report("n0", -1, -1, C.t.n0, n0);
        for (size_t k = 0; k < C.colours.size(); ++k) {
            const int i = C.colours[k];
            i64 ni = 0;
            for (int leaf : ci.leaves)
                if (cx.col(leaf) == i) ++ni;
            if (C.ent[k].n_i != ni) report("n_i", i, -1, C.ent[k].n_i, ni);
            for (const auto& pc : PCS) {
                if (!ci.isC && (pc.p.kind == "up2" || pc.p.kind == "upNode" || pc.p.kind == "upPair" ||
                                pc.p.kind == "pairUp" || pc.p.kind == "up3"))
                    continue;  // not defined for G (must be 0 in our representation)
                const i64 exp = brute(pc.p, ci, cx, lp, i, -1);
                const i64 got = pc.get(C.ent[k]);
                if (got != exp) report(pc.p.name, i, -1, got, exp);
            }
            // pair table row
            for (int r = C.row_begin[k]; r < C.row_begin[k + 1]; ++r) {
                const int j = C.pairs[r].j;
                i64 e1 = brute({"n_ij", "pairNode", "ij"}, ci, cx, lp, i, j);
                i64 e2 = brute({"p_ij", "pairOff", "ij"}, ci, cx, lp, i, j);
                if (C.pairs[r].n_ij != e1) report("n_ij", i, j, C.pairs[r].n_ij, e1);
                if (C.pairs[r].p_ij != e2) report("p_ij", i, j, C.pairs[r].p_ij, e2);
            }
            // pairs missing from the table must be zero
            for (size_t l = k + 1; l < C.colours.size(); ++l) {
                const int j = C.colours[l];
                bool present = false;
                for (int r = C.row_begin[k]; r < C.row_begin[k + 1]; ++r)
                    if (C.pairs[r].j == j) present = true;
                if (present) continue;
                i64 e1 = brute({"n_ij", "pairNode", "ij"}, ci, cx, lp, i, j);
                i64 e2 = brute({"p_ij", "pairOff", "ij"}, ci, cx, lp, i, j);
                if (e1 || e2) report("missing pair", i, j, 0, e1 + e2);
            }
        }
        for (const auto& tc : TTS) {
            if (!ci.isC && (tc.p.kind == "up2" || tc.p.kind == "upNode" || tc.p.kind == "upPair" ||
                            tc.p.kind == "pairUp" || tc.p.kind == "up3"))
                continue;
            const i64 exp = brute(tc.p, ci, cx, lp, -1, -1);
            const i64 got = tc.get(C.t);
            if (got != exp) report(tc.p.name, -1, -1, got, exp);
        }
    }
    return bad;
}

// ---------------------------------------------------------------------------
// Part 2: root contributions per case.
// ---------------------------------------------------------------------------
// T1 type from colours: 0 = α, 1 = β, 2 = γ, -1 = not anchored here.
int t1_type(const std::vector<int>& c, std::vector<int>& pairidx) {
    std::map<int, std::vector<int>> by;
    for (int k = 0; k < 4; ++k) by[c[k]].push_back(k);
    int dbl = -1;
    for (auto& [col, v] : by)
        if (col != 0 && v.size() == 2) {
            if (dbl >= 0) {  // {i,i,j,j}
                pairidx = by[dbl];
                return 0;
            }
            dbl = col;
        }
    if (dbl < 0) return -1;
    pairidx = by[dbl];
    std::vector<int> rest;
    for (int k = 0; k < 4; ++k)
        if (c[k] != dbl) rest.push_back(c[k]);
    if (rest.size() != 2) return -1;
    if (rest[0] != 0 && rest[1] != 0 && rest[0] != rest[1]) return 1;
    if ((rest[0] == 0) != (rest[1] == 0)) return 2;
    return -1;
}

struct Brute2 {
    i64 a[6] = {0, 0, 0, 0, 0, 0}, b[6] = {0, 0, 0, 0, 0, 0};
};

Brute2 brute_root(const Tree& t2, const std::vector<int>& colour) {
    NaiveTopology nt(t2);
    const int n = static_cast<int>(nt.leaf_node.size());
    Brute2 out;
    auto t2_type = [&](int p0, int p1, int q0, int q1) {
        int a = nt.leaf_node[p0], b = nt.leaf_node[p1], c = nt.leaf_node[q0], d = nt.leaf_node[q1];
        int z = nt.lca(nt.lca(a, b), nt.lca(c, d));
        int x = nt.lca(a, b), y = nt.lca(c, d);
        if (x != z && y != z) return 0;
        if (x == z) {  // swap so that y == z
            std::swap(a, c);
            std::swap(b, d);
        }
        // pair (a,b) below z; c,d in different child subtrees of z
        if (nt.lca(nt.lca(a, b), c) == z && nt.lca(nt.lca(a, b), d) == z) return 1;
        return 2;
    };
    for (int a = 0; a < n; ++a)
        for (int b = a + 1; b < n; ++b)
            for (int c = b + 1; c < n; ++c)
                for (int d = c + 1; d < n; ++d) {
                    int ids[4] = {a, b, c, d};
                    std::vector<int> cols = {colour[a], colour[b], colour[c], colour[d]};
                    std::vector<int> pairidx;
                    int ty1 = t1_type(cols, pairidx);
                    if (ty1 < 0) continue;
                    int q = nt.quartet(a, b, c, d);
                    if (q == 3) continue;
                    // T2 split pairs
                    int P0, P1, Q0, Q1;
                    if (q == 0) { P0 = 0; P1 = 1; Q0 = 2; Q1 = 3; }
                    else if (q == 1) { P0 = 0; P1 = 2; Q0 = 1; Q1 = 3; }
                    else { P0 = 0; P1 = 3; Q0 = 1; Q1 = 2; }
                    int ty2 = t2_type(ids[P0], ids[P1], ids[Q0], ids[Q1]);
                    // shared?
                    std::set<int> t1pair(pairidx.begin(), pairidx.end());
                    bool shared = (t1pair == std::set<int>{P0, P1}) || (t1pair == std::set<int>{Q0, Q1});
                    // case index among αα,βα,ββ,γα,γβ,γγ (T1 type >= T2 type), else swapped
                    int cidx = -1;
                    if (ty1 == 0 && ty2 == 0) cidx = 0;
                    else if (ty1 == 1 && ty2 == 0) cidx = 1;
                    else if (ty1 == 1 && ty2 == 1) cidx = 2;
                    else if (ty1 == 2 && ty2 == 0) cidx = 3;
                    else if (ty1 == 2 && ty2 == 1) cidx = 4;
                    else if (ty1 == 2 && ty2 == 2) cidx = 5;
                    if (cidx < 0) continue;  // belongs to the swapped run
                    if (shared) ++out.a[cidx];
                    else ++out.b[cidx];
                }
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    int iters = argc > 1 ? std::atoi(argv[1]) : 200;
    int nmax = argc > 2 ? std::atoi(argv[2]) : 9;
    std::mt19937 rng(777);
    int bad_total = 0;
    static const char* names[6] = {"aa", "ba", "bb", "ga", "gb", "gg"};
    for (int iter = 0; iter < iters; ++iter) {
        int n = 3 + iter % (nmax - 2);
        Tree t = test::random_tree(n, rng, (iter % 4) * 0.33);
        Tree tt = t;
        assign_leaf_ids(t, tt);
        t.contract_unary();
        CTree ct = CTree::from_tree(t);
        int d = 1 + iter % 8;  // up to 8 colours so that the pair-table (★) paths see many colours
        std::vector<int> colour(n);
        std::uniform_int_distribution<int> U(0, d);
        for (int& c : colour) c = U(rng);
        HdtDeco<QuartetFields> deco;
        deco.build(ct, colour, false);
        int bad = check_counters(ct, colour, deco, 10);
        if (bad) {
            std::printf("iter %d n=%d T2=%s colours:", iter, n, t.to_newick().c_str());
            for (int c : colour) std::printf(" %d", c);
            std::printf("\n  %d counter mismatches\n", bad);
            bad_total += bad;
            if (bad_total > 30) break;
            continue;
        }
        // Part 2
        Brute2 br = brute_root(t, colour);
        for (int cidx = 0; cidx < 6; ++cidx) {
            contrib_case_mask() = 1u << cidx;
            HdtDeco<QuartetFields> dd;
            dd.build(ct, colour, false);
            Contrib got = dd.root_sum();
            if (got.a != br.a[cidx] || got.b != br.b[cidx]) {
                ++bad_total;
                std::printf("iter %d n=%d case %s: got A=%s B=%s expected A=%s B=%s\n  T2=%s colours:", iter, n,
                            names[cidx], count_to_string(got.a).c_str(), count_to_string(got.b).c_str(), count_to_string(br.a[cidx]).c_str(), count_to_string(br.b[cidx]).c_str(),
                            t.to_newick().c_str());
                for (int c : colour) std::printf(" %d", c);
                std::printf("\n");
            }
        }
        contrib_case_mask() = CASE_ALL;
        if (bad_total > 30) break;
    }
    std::printf("test_counters: %s (%d problems)\n", bad_total ? "FAILED" : "ok", bad_total);
    return bad_total ? 1 : 0;
}
