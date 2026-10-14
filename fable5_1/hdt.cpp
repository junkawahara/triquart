// hdt.cpp — Sec. 4 HDT construction; see hdt.hpp.
#include "hdt.hpp"

#include <algorithm>
#include <stdexcept>

namespace soda13 {

namespace {

enum class CType : unsigned char { L, I, C, G };

struct Comp {
    CType type = CType::L;
    int hnode = -1;
    int parent = -1;             // component id
    std::vector<int> children;   // component ids
    bool alive = true;
    bool open = false;           // C: has a downward edge
    int redirect = -1;           // when dead: the component that absorbed it
    int tnode = -1;              // I: the T node
};

struct Builder {
    const CTree& t;
    Hdt h;
    std::vector<Comp> comps;
    std::vector<int> seen;  // dedupe stamps for child-list rebuilding
    int stamp = 0;

    explicit Builder(const CTree& tree) : t(tree) {}

    int new_hnode(HOp op, int c1, int c2, int tnode) {
        HNode n;
        n.op = op;
        n.c1 = c1;
        n.c2 = c2;
        n.tnode = tnode;
        const int id = static_cast<int>(h.nodes.size());
        if (c1 >= 0) h.nodes[c1].parent = id;
        if (c2 >= 0) h.nodes[c2].parent = id;
        switch (op) {
            case HOp::LEAF:
                n.weight = t.nodes[tnode].weight;
                n.top = n.bottom = tnode;
                break;
            case HOp::L2C:
            case HOp::C2G:
                n.weight = h.nodes[c1].weight;
                n.top = h.nodes[c1].top;
                n.bottom = h.nodes[c1].bottom;
                n.open = false;
                break;
            case HOp::CC2C:  // c1 = lower, c2 = upper
                n.weight = h.nodes[c1].weight + h.nodes[c2].weight;
                n.top = h.nodes[c2].top;
                n.bottom = h.nodes[c1].bottom;
                n.open = h.nodes[c1].open;
                break;
            case HOp::GG2G:
                n.weight = h.nodes[c1].weight + h.nodes[c2].weight;
                break;
            case HOp::IG2C:
                n.weight = h.nodes[c1].weight;
                n.top = n.bottom = tnode;
                break;
        }
        h.nodes.push_back(n);
        return id;
    }

    int resolve(int c) const {
        while (c >= 0 && !comps[c].alive) c = comps[c].redirect;
        return c;
    }

    void run() {
        const int m = t.num_nodes();
        comps.resize(m);
        h.leaf_hnode.assign(m, -1);
        seen.assign(m * 2 + 2, -1);
        // Before the first round: leaves -> L -> C, internal nodes -> I.
        for (int v = 0; v < m; ++v) {
            Comp& c = comps[v];
            c.parent = t.nodes[v].parent;
            c.children = t.nodes[v].children;
            c.tnode = v;
            if (t.is_leaf(v)) {
                const int leaf = new_hnode(HOp::LEAF, -1, -1, v);
                h.leaf_hnode[v] = leaf;
                c.hnode = new_hnode(HOp::L2C, leaf, -1, v);
                c.type = CType::C;
                c.open = false;
            } else {
                c.type = CType::I;
            }
        }
        std::vector<int> alive;
        for (int v = 0; v < m; ++v) alive.push_back(v);

        int rounds = 0;
        while (alive.size() > 1) {
            ++rounds;
            const size_t before = alive.size();
            // Components created in this round get ids >= first_new; they are the
            // only candidates for the next active set besides the survivors of
            // `alive` (see the end of the round).
            const size_t first_new = comps.size();
            // ---- step 1: C→G for downward-closed C children of I components ----
            for (int c : alive) {
                Comp& cc = comps[c];
                if (cc.type == CType::C && cc.children.empty() && cc.parent >= 0 && comps[cc.parent].type == CType::I) {
                    cc.hnode = new_hnode(HOp::C2G, cc.hnode, -1, -1);
                    cc.type = CType::G;
                }
            }
            // ---- step 2 (iii): CC→C along maximal paths of C components ----
            std::vector<int> chain;
            for (int c : alive) {
                const Comp& cc = comps[c];
                if (cc.type != CType::C) continue;
                if (cc.parent >= 0 && comps[cc.parent].type == CType::C) continue;  // not a chain top
                chain.clear();
                int cur = c;
                while (true) {
                    chain.push_back(cur);
                    const Comp& k = comps[cur];
                    if (k.children.size() == 1 && comps[k.children[0]].type == CType::C)
                        cur = k.children[0];
                    else
                        break;
                }
                for (size_t q = 0; q + 1 < chain.size(); q += 2) {
                    const int upper = chain[q], lower = chain[q + 1];
                    Comp n;
                    n.type = CType::C;
                    n.hnode = new_hnode(HOp::CC2C, comps[lower].hnode, comps[upper].hnode, -1);
                    n.parent = comps[upper].parent;
                    n.children = comps[lower].children;
                    n.open = comps[lower].open;
                    const int id = static_cast<int>(comps.size());
                    comps.push_back(n);
                    comps[upper].alive = false;
                    comps[upper].redirect = id;
                    comps[lower].alive = false;
                    comps[lower].redirect = id;
                }
            }
            // ---- step 2 (i) and (ii): compositions at I components ----
            std::vector<int> gs, opens;
            for (int c : alive) {
                if (comps[c].type != CType::I) continue;
                gs.clear();
                opens.clear();
                for (int k : comps[c].children) {
                    const int r = resolve(k);
                    if (comps[r].type == CType::G)
                        gs.push_back(r);
                    else
                        opens.push_back(r);
                }
                if (gs.size() >= 2) {
                    for (size_t q = 0; q + 1 < gs.size(); q += 2) {
                        Comp n;
                        n.type = CType::G;
                        n.hnode = new_hnode(HOp::GG2G, comps[gs[q]].hnode, comps[gs[q + 1]].hnode, -1);
                        n.parent = c;
                        const int id = static_cast<int>(comps.size());
                        comps.push_back(n);
                        comps[gs[q]].alive = false;
                        comps[gs[q]].redirect = id;
                        comps[gs[q + 1]].alive = false;
                        comps[gs[q + 1]].redirect = id;
                    }
                } else if (gs.size() == 1 && opens.size() <= 1) {
                    Comp n;
                    n.type = CType::C;
                    n.hnode = new_hnode(HOp::IG2C, comps[gs[0]].hnode, -1, comps[c].tnode);
                    n.parent = comps[c].parent;
                    n.children = opens;
                    n.open = !opens.empty();
                    h.nodes[n.hnode].open = n.open;  // the I node keeps the downward edge, if any
                    const int id = static_cast<int>(comps.size());
                    comps.push_back(n);
                    comps[c].alive = false;
                    comps[c].redirect = id;
                    comps[gs[0]].alive = false;
                    comps[gs[0]].redirect = id;
                }
            }
            // ---- rebuild parent / children pointers through redirects ----
            // Lemma 4.1 charges each round O(current component tree): the
            // proof sums a geometric series over the shrinking trees, so the
            // round may not touch the append-only history `comps` as a whole
            // (that costs Theta(m) per round, Theta(m log m) in total on a
            // complete binary tree; trial 8).  The next active set is the
            // survivors of `alive` followed by the components created in this
            // round; both lists are increasing in id and every new id exceeds
            // every old one, so `next` stays sorted and the HDT built here is
            // node-for-node the same as before the fix.
            std::vector<int> next;
            next.reserve(alive.size());
            for (int c : alive)
                if (comps[c].alive) next.push_back(c);
            for (size_t k = first_new; k < comps.size(); ++k)
                if (comps[k].alive) next.push_back(static_cast<int>(k));
            if (seen.size() < comps.size()) seen.resize(comps.size() * 2, -1);
            for (int c : next) {
                Comp& cc = comps[c];
                cc.parent = resolve(cc.parent);
                std::vector<int> nc;
                const int st = ++stamp;
                for (int k : cc.children) {
                    const int r = resolve(k);
                    if (r < 0 || r == c) continue;
                    if (seen[r] == st) continue;
                    seen[r] = st;
                    nc.push_back(r);
                }
                cc.children.swap(nc);
            }
            alive.swap(next);
            if (alive.size() >= before)
                throw std::logic_error("build_hdt: round made no progress (violates Lemma 4.1)");
        }
        h.root = comps[alive[0]].hnode;
        h.height = rounds;
    }
};

}  // namespace

Hdt build_hdt(const CTree& t) {
    Builder b(t);
    b.run();
    return std::move(b.h);
}

}  // namespace soda13
