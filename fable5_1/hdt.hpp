// hdt.hpp — Sec. 4 of soda13.pdf: the hierarchical decomposition tree (HDT)
// of a rooted tree of arbitrary degree.
//
// Components (Fig. 5):  L = a leaf,  I = an internal node,  C = a connected
// set of nodes with at most one edge leaving upwards and at most one leaving
// downwards,  G = a set of whole subtrees whose roots are siblings.
//
// Actions (Fig. 6): the transformations L→C and C→G, and the compositions
// CC→C, GG→G and IG→C.  The paper only materialises compositions as HDT
// nodes; we also keep the transformations as unary HDT nodes because the
// decoration rules of Fig. 12 / Fig. 15 are stated per action.  This at
// most doubles the height (each component is transformed at most once per
// round).
//
// Construction (Sec. 4, "Given a tree T with n nodes ..."): leaves become L
// then C components, internal nodes I components; then rounds of
//   1. C→G for every downward-closed C component that is a child of an I,
//   2. (i) GG→G on ⌊g/2⌋ arbitrary pairs of G children of every I with g ≥ 2,
//      (ii) IG→C for every non-forking I with exactly one G child,
//      (iii) CC→C on ⌊c/2⌋ consecutive pairs (topmost with second, …) of every
//            maximal path of c consecutive C components,
// until a single component remains.  Lemma 4.1: O(n) time, height ≤
// ⌈log_{10/9} |v|⌉ below any node v (local balance).
#pragma once

#include <cstdint>
#include <vector>

#include "ctree.hpp"

namespace soda13 {

enum class HOp : unsigned char { LEAF, L2C, C2G, CC2C, GG2G, IG2C };

struct HNode {
    HOp op = HOp::LEAF;
    int parent = -1;
    int c1 = -1;  // unary ops: the child.  CC2C: C1 (lower).  GG2G: G1.
    int c2 = -1;  // CC2C: C2 (upper).  GG2G: G2.
    int tnode = -1;          // LEAF/L2C: the CTree leaf;  IG2C: the I node
    bool open = false;       // C components: has a downward edge
    std::int64_t weight = 0; // total leaf weight (real leaves + zero-leaf weights)
    int top = -1, bottom = -1;  // C components: ends of the external path (CTree nodes)
    bool is_C() const { return op == HOp::L2C || op == HOp::CC2C || op == HOp::IG2C; }
    bool is_G() const { return op == HOp::C2G || op == HOp::GG2G; }
};

struct Hdt {
    std::vector<HNode> nodes;  // children always precede parents (post-order)
    int root = -1;
    std::vector<int> leaf_hnode;  // CTree node -> LEAF hnode, -1 for internal nodes
    int height = 0;
};

Hdt build_hdt(const CTree& t);

}  // namespace soda13
