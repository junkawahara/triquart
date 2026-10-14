// count.hpp — Sec. 5 of soda13.pdf: the main algorithm COUNT(v) (Fig. 7)
// with the `contract` and `extract` compressions of T2.
//
//   COUNT(v):
//     if v is a leaf: colour v by 0
//     else
//       let c_1 be the child of v with the largest subtree, c_2..c_k the rest
//       for i = 2..k: colour the leaves of c_i by i        (leaves are now coloured according to v)
//       query the HDT for the number of triplets/quartets in T2 compatible with the colouring
//       for i = 2..k: colour the leaves of c_i by 0
//       COUNT(c_1)
//       for i = 2..k: colour the leaves of c_i by 1; COUNT(c_i)
//
// Invariant (Sec. 5): when entering v all leaves below v have colour 1 and
// all others colour 0; on exit everything is 0.
//
// Compression (Sec. 5, "We now outline how to perform the compression"):
//   * contract — "performed just before the recursive call COUNT(c_1), if
//     the size of c_1 is less than some fixed fraction of the size of (the
//     current version of) T2"; builds a tree of size O(|c_1|) containing the
//     leaves of c_1 with the same induced topologies, the other leaves being
//     merged into colour-0 nodes.  Implemented as the Steiner-tree
//     compression `contract_tree`, O(|T2|) time.
//   * extract — for every c_i, i ≥ 2, "just after the query to the HDT":
//     "finding the union of k_i root-to-leaf paths in the HDT from the leaves
//     of c_i to the root ... all subtrees hanging to sides of this union of
//     paths are converted to single nodes, with the counting information
//     reset to represent the situation where all leaves in these subtrees
//     have color 0", followed by contract.  O(k_i + k_i log(|v|/k_i)).
//
// Zero-leaves are only kept in quartet mode: the triplet counters of Fig. 10
// never involve colour 0 (they count {i,i,j} and {i,j,k}).
#pragma once

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <vector>

#include "ctree.hpp"
#include "hdt_deco.hpp"
#include "tree.hpp"

namespace soda13 {

// Debug aid: when true, COUNT never contracts / extracts and follows Fig. 7
// literally on the full T2 (O(n^2)-ish).  Used by the tests to separate
// errors in the counters from errors in the compression.
inline bool& count_no_compress() {
    static bool flag = false;
    return flag;
}

namespace detail {

// Steiner-tree compression of `in` with respect to the leaves accepted by
// `designated(leaf_id)`.  Non-designated leaves become zero-leaf weights
// attached at the Steiner node / edge they hang from (dropped entirely if
// !keep_zero).  Output nodes always have >= 2 children.
template <class Pred>
inline CTree contract_tree(const CTree& in, Pred designated, bool keep_zero) {
    const int m = in.num_nodes();
    std::vector<int> order;
    order.reserve(m);
    std::vector<int> st{in.root};
    while (!st.empty()) {
        int v = st.back();
        st.pop_back();
        order.push_back(v);
        for (int c : in.nodes[v].children) st.push_back(c);
    }
    std::vector<int> desig(m, 0);
    std::vector<i64> wt(m, 0);
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        const int v = *it;
        const CTree::Node& nd = in.nodes[v];
        if (nd.children.empty()) {
            if (nd.leaf_id >= 0 && designated(nd.leaf_id))
                desig[v] = 1;
            else
                wt[v] = nd.weight;
        } else {
            for (int c : nd.children) {
                desig[v] += desig[c];
                wt[v] += wt[c];
            }
        }
    }
    CTree out;
    // Descend from `c` along chain nodes; returns the first kept node and the
    // zero weight hanging off the chain strictly above it.
    auto descend = [&](int c, i64& ew) {
        ew = 0;
        while (!in.nodes[c].children.empty()) {
            int next = -1, cnt = 0;
            for (int k : in.nodes[c].children)
                if (desig[k] > 0) {
                    ++cnt;
                    next = k;
                }
            if (cnt != 1) break;  // branching node: kept
            ew += wt[c] - wt[next];
            c = next;
        }
        return c;
    };
    // emit(v, parent): v is a kept node (designated leaf or branching node).
    struct Frame {
        int v, parent;
    };
    std::vector<Frame> work;
    auto attach_chain = [&](int c, int parent) {
        i64 ew = 0;
        const int kept = descend(c, ew);
        int par = parent;
        if (keep_zero && ew > 0) {
            par = out.add_node();
            out.attach(parent, par);
            out.add_zero_leaf(par, ew);
        }
        work.push_back({kept, par});
        return par;
    };
    // root handling
    {
        i64 ew = 0;
        const int kept = descend(in.root, ew);
        if (keep_zero && ew > 0) {
            out.root = out.add_node();
            out.add_zero_leaf(out.root, ew);
            work.push_back({kept, out.root});
        } else {
            work.push_back({kept, -1});
        }
    }
    while (!work.empty()) {
        Frame f = work.back();
        work.pop_back();
        const CTree::Node& nd = in.nodes[f.v];
        if (nd.children.empty()) {
            const int leaf = out.add_real_leaf(f.parent, nd.leaf_id);
            if (f.parent < 0) out.root = leaf;
            continue;
        }
        const int p = out.add_node();
        out.attach(f.parent, p);
        if (f.parent < 0) out.root = p;
        i64 here = wt[f.v];
        for (int c : nd.children)
            if (desig[c] > 0) here -= wt[c];
        if (keep_zero && here > 0) out.add_zero_leaf(p, here);
        for (int c : nd.children)
            if (desig[c] > 0) attach_chain(c, p);
    }
    return out;
}

}  // namespace detail

template <class F>
class Engine {
public:
    // Counts over (T1, T2).  `asym_only` restricts the quartet contributions
    // to the cases βα, γα, γβ (used for the swapped run, Fig. 8 / Fig. 9).
    Engine(const Tree& t1, const Tree& t2, bool asym_only)
        : t1_(t1), t2_(t2), asym_only_(asym_only), keep_zero_(F::kQuartet) {}

    Contrib run() {
        prepare_t1();
        const int n = static_cast<int>(order_.size());
        if (n < 3) return Contrib();
        colour_.assign(n, 1);
        cur_tnode_.assign(n, -1);
        auto inst = std::make_unique<Instance>();
        inst->tree = CTree::from_tree(t2_);
        for (int v = 0; v < inst->tree.num_nodes(); ++v)
            if (inst->tree.is_leaf(v)) inst->leaf_map.emplace_back(inst->tree.nodes[v].leaf_id, v);
        finish_instance(*inst);
        activate(*inst);
        count(t1_.root, *inst);
        // Sec. 5 invariant, "when exiting v, all leaves in T1 have the colour 0",
        // applied to the root call (guard in the spirit of the Lemma 4.1 check).
        for (int c : colour_)
            if (c != 0) throw std::logic_error("COUNT: exit invariant violated (a leaf is not coloured 0)");
        return total_;
    }

private:
    struct Instance {
        CTree tree;
        HdtDeco<F> deco;
        std::vector<std::pair<int, int>> leaf_map;  // (leaf id, CTree node)
        std::vector<int> hmark;                     // extraction stamps per HDT node
        int stamp = 0;
    };

    void prepare_t1() {
        const int m = t1_.num_nodes();
        size_.assign(m, 0);
        lo_.assign(m, 0);
        hi_.assign(m, 0);
        children_.assign(m, {});
        order_.clear();
        // iterative DFS computing Euler leaf ranges and subtree sizes
        std::vector<std::pair<int, int>> st{{t1_.root, 0}};
        while (!st.empty()) {
            auto& [v, k] = st.back();
            if (k == 0) lo_[v] = static_cast<int>(order_.size());
            if (k < static_cast<int>(t1_.nodes[v].children.size())) {
                int c = t1_.nodes[v].children[k++];
                st.emplace_back(c, 0);
            } else {
                if (t1_.is_leaf(v)) order_.push_back(t1_.nodes[v].leaf_id);
                hi_[v] = static_cast<int>(order_.size());
                size_[v] = hi_[v] - lo_[v];
                st.pop_back();
            }
        }
        for (int v = 0; v < m; ++v) {
            children_[v] = t1_.nodes[v].children;
            std::stable_sort(children_[v].begin(), children_[v].end(),
                             [&](int a, int b) { return size_[a] > size_[b]; });
        }
    }

    void finish_instance(Instance& inst) {
        for (auto [lid, tn] : inst.leaf_map) colour_[lid] = 1;
        inst.deco.build(inst.tree, colour_, asym_only_);
        inst.hmark.assign(inst.deco.hdt().nodes.size(), 0);
    }

    void activate(Instance& inst) {
        for (auto [lid, tn] : inst.leaf_map) {
            cur_tnode_[lid] = tn;
            colour_[lid] = 1;
        }
    }

    void recolour_range(Instance& inst, int lo, int hi, int colour) {
        for (int k = lo; k < hi; ++k) {
            const int lid = order_[k];
            colour_[lid] = colour;
            inst.deco.mark_leaf(cur_tnode_[lid]);
        }
    }

    // ---- Sec. 5, Fig. 7 ----
    // The recursion of Fig. 7 is as deep as T1 (a caterpillar gives depth n),
    // so it is run on an explicit stack of frames.
    struct Frame {
        int v;
        Instance* inst;                               // instance used by this call
        std::unique_ptr<Instance> owned;              // if this call owns it
        std::vector<std::unique_ptr<Instance>> subs;  // extracted trees for c_2..c_k
        int next = 0;                                 // next child index to visit
        bool started = false;
    };

    void count(int root, Instance& root_inst) {
        std::vector<Frame> stack;
        stack.push_back(Frame{root, &root_inst, nullptr, {}, 0, false});
        while (!stack.empty()) {
            Frame& f = stack.back();
            const int v = f.v;
            if (t1_.is_leaf(v)) {
                stack.pop_back();
                continue;
            }
            const std::vector<int>& ch = children_[v];
            const int k = static_cast<int>(ch.size());
            Instance& inst = *f.inst;
            if (!f.started) {
                f.started = true;
                for (int i = 1; i < k; ++i) recolour_range(inst, lo_[ch[i]], hi_[ch[i]], i + 1);
                inst.deco.refresh();
                total_ += inst.deco.root_sum();
                if (!count_no_compress()) {
                    // extract the trees for c_2..c_k "just after the query", while
                    // their leaves are still coloured i (Sec. 5)
                    f.subs.resize(k);
                    for (int i = 1; i < k; ++i)
                        if (!t1_.is_leaf(ch[i])) f.subs[i] = extract(inst, lo_[ch[i]], hi_[ch[i]]);
                }
                // Fig. 7: "for i = 2 to k: colour the leaves in the subtree of c_i by
                // the colour 0".  A leaf c_1 is coloured 0 here as well, which is what
                // COUNT(c_1) does for a leaf.  This re-establishes the Sec. 5 invariant
                // ("when exiting v, all leaves in T1 have the colour 0") on this
                // instance even when c_1 is contracted and the instance is never
                // queried again; the marks are consumed by the next refresh().
                for (int i = 1; i < k; ++i) recolour_range(inst, lo_[ch[i]], hi_[ch[i]], 0);
                if (t1_.is_leaf(ch[0])) recolour_range(inst, lo_[ch[0]], hi_[ch[0]], 0);
                f.next = 0;
            }
            if (f.next >= k) {
                stack.pop_back();
                continue;
            }
            const int i = f.next++;
            const int c = ch[i];
            if (t1_.is_leaf(c)) continue;
            if (count_no_compress()) {
                if (i >= 1) recolour_range(inst, lo_[c], hi_[c], 1);
                stack.push_back(Frame{c, &inst, nullptr, {}, 0, false});
                continue;
            }
            if (i == 0) {
                if (2 * size_[c] < inst.tree.num_designated) {
                    auto sub = contract(inst, lo_[c], hi_[c]);
                    activate(*sub);
                    Instance* p = sub.get();
                    stack.push_back(Frame{c, p, std::move(sub), {}, 0, false});
                } else {
                    // COUNT(c_1) on the shared instance (c_2..c_k are already 0)
                    stack.push_back(Frame{c, &inst, nullptr, {}, 0, false});
                }
            } else {
                std::unique_ptr<Instance> sub = std::move(f.subs[i]);
                activate(*sub);
                Instance* p = sub.get();
                stack.push_back(Frame{c, p, std::move(sub), {}, 0, false});
            }
        }
    }

    // ---- contract (Sec. 5) ----
    std::unique_ptr<Instance> contract(Instance& src, int lo, int hi) {
        ++desig_stamp_;
        if (desig_mark_.size() < colour_.size()) desig_mark_.assign(colour_.size(), 0);
        for (int k = lo; k < hi; ++k) desig_mark_[order_[k]] = desig_stamp_;
        auto out = std::make_unique<Instance>();
        const int stamp = desig_stamp_;
        out->tree = detail::contract_tree(
            src.tree, [&](int lid) { return desig_mark_[lid] == stamp; }, keep_zero_);
        collect_leaf_map(*out);
        finish_instance(*out);
        return out;
    }

    // ---- extract (Sec. 5) ----
    std::unique_ptr<Instance> extract(Instance& src, int lo, int hi) {
        const Hdt& h = src.deco.hdt();
        const int stamp = ++src.stamp;
        for (int k = lo; k < hi; ++k) {
            int x = src.deco.leaf_hnode(cur_tnode_[order_[k]]);
            while (x >= 0 && src.hmark[x] != stamp) {
                src.hmark[x] = stamp;
                x = h.nodes[x].parent;
            }
        }
        CTree ex;
        ExtractCtx ctx{src, ex, stamp};
        auto piece = build_C(ctx, h.root);
        ex.root = piece.first;
        auto out = std::make_unique<Instance>();
        out->tree = detail::contract_tree(
            ex, [](int) { return true; }, keep_zero_);
        collect_leaf_map(*out);
        finish_instance(*out);
        return out;
    }

    struct ExtractCtx {
        Instance& src;
        CTree& ex;
        int stamp;
    };

    // Rebuild the piece of T2 represented by the C component `x`: returns
    // (top, bottom) nodes of its external path in the new tree.
    std::pair<int, int> build_C(ExtractCtx& ctx, int x) {
        const Hdt& h = ctx.src.deco.hdt();
        const HNode& n = h.nodes[x];
        if (ctx.src.hmark[x] != ctx.stamp) {
            // side component: all its leaves have colour 0
            if (n.open) {
                const int p = ctx.ex.add_node();
                if (keep_zero_ && n.weight > 0) ctx.ex.add_zero_leaf(p, n.weight);
                return {p, p};
            }
            const int z = ctx.ex.add_node();
            ctx.ex.nodes[z].weight = keep_zero_ ? n.weight : 0;
            return {z, z};
        }
        switch (n.op) {
            case HOp::LEAF: {
                const CTree::Node& ln = ctx.src.tree.nodes[n.tnode];
                const int leaf = ctx.ex.add_node();
                ctx.ex.nodes[leaf].leaf_id = ln.leaf_id;
                ctx.ex.nodes[leaf].weight = 1;
                ++ctx.ex.num_designated;
                return {leaf, leaf};
            }
            case HOp::L2C:
                return build_C(ctx, n.c1);
            case HOp::CC2C: {
                auto up = build_C(ctx, n.c2);
                auto low = build_C(ctx, n.c1);
                ctx.ex.attach(up.second, low.first);
                return {up.first, low.second};
            }
            case HOp::IG2C: {
                const int p = ctx.ex.add_node();
                build_G(ctx, n.c1, p);
                return {p, p};
            }
            default:
                break;
        }
        return {-1, -1};  // unreachable
    }

    void build_G(ExtractCtx& ctx, int x, int parent) {
        const Hdt& h = ctx.src.deco.hdt();
        const HNode& n = h.nodes[x];
        if (ctx.src.hmark[x] != ctx.stamp) {
            if (keep_zero_ && n.weight > 0) ctx.ex.add_zero_leaf(parent, n.weight);
            return;
        }
        if (n.op == HOp::C2G) {
            auto piece = build_C(ctx, n.c1);
            ctx.ex.attach(parent, piece.first);
        } else {  // GG2G
            build_G(ctx, n.c1, parent);
            build_G(ctx, n.c2, parent);
        }
    }

    void collect_leaf_map(Instance& inst) {
        inst.leaf_map.clear();
        for (int v = 0; v < inst.tree.num_nodes(); ++v)
            if (inst.tree.is_leaf(v) && inst.tree.nodes[v].leaf_id >= 0)
                inst.leaf_map.emplace_back(inst.tree.nodes[v].leaf_id, v);
    }

    const Tree& t1_;
    const Tree& t2_;
    bool asym_only_;
    bool keep_zero_;
    std::vector<int> order_, lo_, hi_, size_;
    std::vector<std::vector<int>> children_;
    std::vector<int> colour_, cur_tnode_;
    std::vector<int> desig_mark_;
    int desig_stamp_ = 0;
    Contrib total_;
};

}  // namespace soda13
