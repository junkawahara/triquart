// hdt_deco.hpp — the HDT decorated with counters (Sec. 5–7 of soda13.pdf).
//
// Every HDT node stores the counters of its component (counters.hpp) and,
// for CC→C and GG→G nodes, its contribution (contributions.hpp); each node
// also stores the sum of the contributions in its subtree, so that "the
// value at the root then contains A(τ_v ∩ τ_u)" summed over all u (Sec. 2).
//
// Recolouring (Sec. 5): "The total updating work in the first for-loop in
// COUNT can be performed by first marking the paths in the HDT from all
// leaves recolored, and then propagating the change of information in the
// nodes of the HDT in a bottom-up manner, passing each marked node only
// once."  `mark()` marks a root path dirty, `refresh()` recomputes every
// dirty node bottom-up exactly once.
#pragma once

#include <vector>

#include "contributions.hpp"
#include "counters.hpp"
#include "hdt.hpp"

namespace soda13 {

template <class F>
class HdtDeco {
public:
    using Counters = CountersT<F>;

    // `colour` is indexed by leaf id and must outlive this object.
    void build(const CTree& tree, const std::vector<int>& colour, bool asym_only) {
        tree_ = &tree;
        colour_ = &colour;
        asym_only_ = asym_only;
        hdt_ = build_hdt(tree);
        const int m = static_cast<int>(hdt_.nodes.size());
        ctr_.assign(m, Counters());
        own_.assign(m, Contrib());
        sum_.assign(m, Contrib());
        dirty_.assign(m, 0);
        for (int h = 0; h < m; ++h) recompute(h);  // children precede parents
    }

    const Hdt& hdt() const { return hdt_; }
    const CTree& tree() const { return *tree_; }
    // L→C and IG→C are the identity on the counters (Fig. 12 / Fig. 15
    // columns "n^C_-" and "n^G_-"), so those nodes share their child's
    // counter set instead of storing a copy.
    const Counters& counters(int h) const {
        while (hdt_.nodes[h].op == HOp::L2C || hdt_.nodes[h].op == HOp::IG2C) h = hdt_.nodes[h].c1;
        return ctr_[h];
    }
    const Contrib& own(int h) const { return own_[h]; }
    Contrib root_sum() const { return sum_[hdt_.root]; }
    int leaf_hnode(int tnode) const { return hdt_.leaf_hnode[tnode]; }

    // Mark the root path of a CTree leaf dirty (after its colour changed).
    void mark_leaf(int tnode) {
        int h = hdt_.leaf_hnode[tnode];
        while (h >= 0 && !dirty_[h]) {
            dirty_[h] = 1;
            h = hdt_.nodes[h].parent;
        }
    }

    void refresh() {
        if (dirty_[hdt_.root]) refresh_rec(hdt_.root);
    }

    // Recompute one node from its children (or from its leaf colour).
    void recompute(int h) {
        const HNode& n = hdt_.nodes[h];
        Counters& out = ctr_[h];
        switch (n.op) {
            case HOp::LEAF: {
                const CTree::Node& ln = tree_->nodes[n.tnode];
                const int col = ln.leaf_id >= 0 ? (*colour_)[ln.leaf_id] : 0;
                compose_leaf(out, col, ln.weight);
                break;
            }
            case HOp::L2C:
            case HOp::IG2C:
                break;  // identity: shared with the child, see counters()
            case HOp::C2G:
                compose_C2G(out, counters(n.c1));
                break;
            case HOp::CC2C:
                compose_CC2C(out, counters(n.c1), counters(n.c2), sc_);
                own_[h] = contrib_CC2C(counters(n.c1), counters(n.c2), asym_only_, sc_);
                break;
            case HOp::GG2G:
                compose_GG2G(out, counters(n.c1), counters(n.c2), sc_);
                own_[h] = contrib_GG2G(counters(n.c1), counters(n.c2), asym_only_, sc_);
                break;
        }
        Contrib s = own_[h];
        if (n.c1 >= 0) s += sum_[n.c1];
        if (n.c2 >= 0) s += sum_[n.c2];
        sum_[h] = s;
        dirty_[h] = 0;
    }

private:
    void refresh_rec(int h) {
        const HNode& n = hdt_.nodes[h];
        if (n.c1 >= 0 && dirty_[n.c1]) refresh_rec(n.c1);
        if (n.c2 >= 0 && dirty_[n.c2]) refresh_rec(n.c2);
        recompute(h);
    }

    const CTree* tree_ = nullptr;
    const std::vector<int>* colour_ = nullptr;
    bool asym_only_ = false;
    Hdt hdt_;
    std::vector<Counters> ctr_;
    std::vector<Contrib> own_, sum_;
    std::vector<char> dirty_;
    Scratch sc_;
};

}  // namespace soda13
