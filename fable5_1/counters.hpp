// counters.hpp — Sec. 7 of soda13.pdf: the counters stored in every node
// (component) of the hierarchical decomposition of T2, and their
// composition rules for the five actions of Fig. 6
//
//     L→C,  C→G,  CC→C,  GG→G,  IG→C.
//
// The counters are those of Fig. 11 (shared triplets / quartets) and
// Fig. 14 (disagreeing quartets); the composition rules are Fig. 12 and
// Fig. 15.  A plain-text transcription of all figures, including the
// handful of printed typos we had to correct, is in `paper_figures.md`;
// every rule below quotes the figure row it implements.
//
// Naming of counters in code (see "Notation" in paper_figures.md):
//     0  = colour zero,  i = the colour of the entry,  b = • ,  s = □
//     n_xy      = n_{xy}        pair/triple branching at the path node / super root
//     p_xy      = n_{(xy)}      LCA not on the external path / not the super root
//     n_x_pyz   = n_{x(yz)}     x at the path node, (yz) below it off-path
//     p_x_pyz   = n_{(x(yz))}
//     br_x_pyz  = n_{[x(yz)]}   bracket: no path constraint
//     n_xUy     = n_{x↑y}       x branches off a lower path node than y
//     n_pxyUz   = n_{(xy)↑z},  n_xUpyz = n_{x↑(yz)},  n_xUyUz = n_{x↑y↑z}
//
// Per-colour counters (those whose subscript contains i) are stored
// sparsely: only colours that actually occur in the component have an
// entry, so composing a node costs O((1+x)^2) where x is the number of
// distinct non-zero colours in it (Sec. 7.1 / 7.4).  Two-colour counters
// n_{ij} and n_{(ij)} (Fig. 14, marked ★) are stored as a sorted sparse
// pair table.
//
// `TripletFields` keeps only the † counters of Fig. 11 (needed for the
// O(n log n) triplet distance); `QuartetFields` keeps everything.
#pragma once

#include <algorithm>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace soda13 {

using i64 = std::int64_t;

// ---------------------------------------------------------------------------
// Per-colour counter fields.
// ---------------------------------------------------------------------------
#define SODA13_TRIPLET_FIELDS(X) X(n_i) X(n_ib) X(p_ii) X(n_iUb)

#define SODA13_QUARTET_EXTRA_FIELDS(X)                                                                 \
    /* pairs (Fig. 11 col. 1-2) */                                                                     \
    X(n_0i) X(n_ii) X(p_0i) X(p_ib)                                                                    \
    /* triples at the path node / super root (Fig. 11 col. 2, Fig. 14 col. 1) */                     \
    X(n_0_pii) X(n_b_pii) X(n_i_p0b) X(n_i_pbb) X(n_i_pbs)                                             \
    X(n_i_p0i) X(n_b_p0i) X(n_0_pib) X(n_b_pis) X(n_i_pib) X(n_b_pib)                                  \
    /* fully parenthesised triples */                                                                  \
    X(p_0_pii) X(p_b_pii) X(p_i_p0b) X(p_i_pib) X(p_i_p0i) X(p_0_pib) X(p_b_p0i)                       \
    /* bracket triples */                                                                              \
    X(br_0_pii) X(br_b_pii) X(br_i_p0b) X(br_i_pib) X(br_0_pib) X(br_i_p0i) X(br_b_p0i)                \
    /* ↑ pairs (C only) */                                                                             \
    X(n_0Ui) X(n_iU0) X(n_iUi) X(n_bUi)                                                                \
    /* leaf ↑ pair-at-node (C only) */                                                                 \
    X(n_iU0b) X(n_bUii) X(n_iUbb) X(n_iUbs) X(n_0Uib) X(n_bU0i) X(n_iUib) X(n_0Uii) X(n_iU0i) X(n_bUis) X(n_bUib) X(n_bUpib) \
    /* leaf ↑ (pair) (C only) */                                                                       \
    X(n_0Upii) X(n_iUp0b) X(n_iUpbs) X(n_bUpii) X(n_iUpbb) X(n_0Upib) X(n_iUpib) X(n_bUp0i) X(n_iUp0i) \
    X(n_bUpis)                                                                                         \
    /* (pair) ↑ leaf (C only) */                                                                       \
    X(n_piiU0) X(n_piiUb) X(n_pbbUi) X(n_p0bUi) X(n_pibU0) X(n_pibUi) X(n_p0iUi) X(n_p0iUb)            \
    /* leaf ↑ leaf ↑ leaf (C only) */                                                                  \
    X(n_0UiUi) X(n_bUiUi) X(n_iUbU0) X(n_iU0Ub) X(n_bUiU0) X(n_iUiU0) X(n_0UbUi) X(n_bU0Ui) X(n_iU0Ui) \
    X(n_0UiUb) X(n_iUiUb) X(n_iUbUi)

struct TripletFields {
    static constexpr bool kQuartet = false;
#define X(f) i64 f = 0;
    SODA13_TRIPLET_FIELDS(X)
#undef X
};

struct QuartetFields {
    static constexpr bool kQuartet = true;
#define X(f) i64 f = 0;
    SODA13_TRIPLET_FIELDS(X)
    SODA13_QUARTET_EXTRA_FIELDS(X)
#undef X
};

// One row entry of the sparse pair table: colours (i, j) with i < j.
struct PairEntry {
    int j = 0;
    i64 n_ij = 0;  // n_{ij}   (Fig. 14 ★): i and j at the path node / super root
    i64 p_ij = 0;  // n_{(ij)} (Fig. 14 ★): LCA off the path / below the super root
};

// Colour-independent counters.  The triplet variant only needs n_0, n_• and
// n_{•□}; everything else is quartet-only (Fig. 11 / Fig. 14).
template <bool Quartet>
struct TotalsT;

template <>
struct TotalsT<false> {
    i64 n0 = 0;    // n_0
    i64 nb = 0;    // n_•   = Σ n_i
    i64 n_bs = 0;  // n_{•□} = ½ Σ n_{i•}
};

template <>
struct TotalsT<true> {
    i64 n0 = 0;       // n_0
    i64 nb = 0;       // n_•   = Σ n_i
    i64 n_bs = 0;     // n_{•□} = ½ Σ n_{i•}
    i64 n_0b = 0;     // n_{0•}
    i64 n_bb = 0;     // n_{••}
    i64 p_0b = 0;     // n_{(0•)}
    i64 p_bb = 0;     // n_{(••)}
    i64 p_bs = 0;     // n_{(•□)}
    i64 n_s_pbb = 0;  // n_{□(••)}
    i64 br_bs = 0;    // n_{[•□]}
    i64 n_0Ub = 0;    // n_{0↑•}
    i64 n_bU0 = 0;    // n_{•↑0}
    i64 n_0Ubb = 0;   // n_{0↑••}
    i64 n_bUpss = 0;  // n_{•↑(□□)} = n_{□↑(••)}
    i64 n_b_pbs = 0;  // n_{•(•□)}
    i64 n_bUbs = 0;   // n_{•↑•□}
    i64 n_pbsUb = 0;  // n_{(•□)↑•}
    i64 n_bUbUs = 0;  // n_{•↑•↑□}
    i64 n_bUsUb = 0;  // n_{•↑□↑•}
    i64 p_b_pbs = 0;  // n_{(•(•□))}
    i64 br_b_pbs = 0; // n_{[•(•□)]}
    i64 n_sUbb = 0;   // n_{□↑••}
    i64 n_bUpbs = 0;  // n_{•↑(•□)}
};

// Sparse pair table (quartet mode only): row k (colour colours[k]) is
// pairs[row_begin[k] .. row_begin[k+1]), sorted by j, holding only j > colours[k].
struct PairTable {
    std::vector<PairEntry> pairs;
    std::vector<int> row_begin;
};
struct NoPairTable {};

// ---------------------------------------------------------------------------
// The counter set of one component.
// ---------------------------------------------------------------------------
template <class F>
struct CountersT : std::conditional_t<F::kQuartet, PairTable, NoPairTable> {
    using Fields = F;

    using Totals = TotalsT<F::kQuartet>;
    Totals t;

    // Sparse per-colour entries, sorted by colour (non-zero colours only).
    std::vector<int> colours;
    std::vector<F> ent;

    void reset_keep_capacity() {
        t = Totals{};
        colours.clear();
        ent.clear();
        if constexpr (F::kQuartet) {
            this->pairs.clear();
            this->row_begin.clear();
        }
    }
    int find(int colour) const {
        auto it = std::lower_bound(colours.begin(), colours.end(), colour);
        if (it == colours.end() || *it != colour) return -1;
        return static_cast<int>(it - colours.begin());
    }
    i64 n_of(int colour) const {
        int k = find(colour);
        return k < 0 ? 0 : ent[k].n_i;
    }
};

using TripletCounters = CountersT<TripletFields>;
using QuartetCounters = CountersT<QuartetFields>;

// Scratch space shared by all compositions (indexed by colour).  Colours are
// bounded by the maximum degree d of T1 (Sec. 5), so this stays O(d).
struct Scratch {
    std::vector<i64> n1;    // n_i of the first operand, by colour
    std::vector<int> pos;   // position of a colour in the output entry list
    void ensure(int colour) {
        if (colour >= static_cast<int>(n1.size())) {
            n1.resize(colour + 1, 0);
            pos.resize(colour + 1, -1);
        }
    }
};

inline i64 c2(i64 n) { return n * (n - 1) / 2; }

// ---------------------------------------------------------------------------
// Internal helpers.
// ---------------------------------------------------------------------------
namespace detail {

// Walk the union of the two sorted colour lists; `fn(colour, a, b, out)` is
// called with zero entries for absent colours.
template <class F, class Fn>
inline void merge_colours(const CountersT<F>& A, const CountersT<F>& B, CountersT<F>& out, Fn fn) {
    static const F zero{};
    const size_t na = A.colours.size(), nbn = B.colours.size();
    out.colours.reserve(na + nbn);
    out.ent.reserve(na + nbn);
    size_t ia = 0, ib = 0;
    while (ia < na || ib < nbn) {
        int c;
        const F* a = &zero;
        const F* b = &zero;
        if (ib >= nbn || (ia < na && A.colours[ia] < B.colours[ib])) {
            c = A.colours[ia];
            a = &A.ent[ia++];
        } else if (ia >= na || B.colours[ib] < A.colours[ia]) {
            c = B.colours[ib];
            b = &B.ent[ib++];
        } else {
            c = A.colours[ia];
            a = &A.ent[ia++];
            b = &B.ent[ib++];
        }
        out.colours.push_back(c);
        out.ent.emplace_back();
        fn(c, *a, *b, out.ent.back());
    }
}

// Box sums of Fig. 12 / Fig. 15: the colour-independent counters that are
// defined as sums over the per-colour ones.
template <class F>
inline void finalize_box_sums(CountersT<F>& o) {
    i64 nb = 0, n_bs2 = 0;
    for (const F& e : o.ent) {
        nb += e.n_i;
        n_bs2 += e.n_ib;
    }
    o.t.nb = nb;
    o.t.n_bs = n_bs2 / 2;
    if constexpr (F::kQuartet) {
        i64 n_0b = 0, n_bb = 0, p_0b = 0, p_bb = 0, p_bs2 = 0, n_s_pbb = 0, br_bs2 = 0, n_0Ub = 0, n_bU0 = 0,
            n_bUpss = 0, n_b_pbs = 0, n_bUbs = 0, n_pbsUb = 0, n_bUbUs = 0, n_bUsUb = 0, p_b_pbs = 0, br_b_pbs = 0,
            n_sUbb = 0, n_bUpbs = 0;
        for (const F& e : o.ent) {
            n_0b += e.n_0i;
            n_bb += e.n_ii;
            p_0b += e.p_0i;
            p_bb += e.p_ii;
            p_bs2 += e.p_ib;
            n_s_pbb += e.n_i_pbb;
            br_bs2 += e.n_i * (nb - e.n_i);  // n_{[i•]} = n_i (n_• − n_i)
            n_0Ub += e.n_0Ui;
            n_bU0 += e.n_iU0;
            n_bUpss += e.n_bUpii;
            n_b_pbs += e.n_i_pib;
            n_bUbs += e.n_iUib;
            n_pbsUb += e.n_pibUi;
            n_bUbUs += e.n_iUiUb;
            n_bUsUb += e.n_iUbUi;
            p_b_pbs += e.p_i_pib;
            br_b_pbs += e.br_i_pib;
            n_sUbb += e.n_bUii;
            n_bUpbs += e.n_iUpib;
        }
        o.t.n_0b = n_0b;
        o.t.n_bb = n_bb;
        o.t.p_0b = p_0b;
        o.t.p_bb = p_bb;
        o.t.p_bs = p_bs2 / 2;
        o.t.n_s_pbb = n_s_pbb;
        o.t.br_bs = br_bs2 / 2;
        o.t.n_0Ub = n_0Ub;
        o.t.n_bU0 = n_bU0;
        o.t.n_bUpss = n_bUpss;
        o.t.n_b_pbs = n_b_pbs;
        o.t.n_bUbs = n_bUbs;
        o.t.n_pbsUb = n_pbsUb;
        o.t.n_bUbUs = n_bUbUs;
        o.t.n_bUsUb = n_bUsUb;
        o.t.p_b_pbs = p_b_pbs;
        o.t.br_b_pbs = br_b_pbs;
        o.t.n_sUbb = n_sUbb;
        o.t.n_bUpbs = n_bUpbs;
    }
}

// Merge the pair tables of A and B (sums), optionally adding the GG→G cross
// term  n_{ij} += n_i^{G1} n_j^{G2} + n_j^{G1} n_i^{G2}  (Fig. 15 group A).
template <class F>
inline void merge_pairs(const CountersT<F>& A, const CountersT<F>& B, CountersT<F>& out, bool cross) {
    if constexpr (!F::kQuartet) {
        (void)A;
        (void)B;
        (void)out;
        (void)cross;
        return;
    } else {
        out.row_begin.reserve(out.colours.size() + 1);
        size_t ia = 0, ib = 0;  // row cursors in A and B
        for (size_t k = 0; k < out.colours.size(); ++k) {
            const int c = out.colours[k];
            out.row_begin.push_back(static_cast<int>(out.pairs.size()));
            // rows of A / B for colour c (if present)
            const PairEntry *abeg = nullptr, *aend = nullptr, *bbeg = nullptr, *bend = nullptr;
            i64 nA_c = 0, nB_c = 0;
            if (ia < A.colours.size() && A.colours[ia] == c) {
                abeg = A.pairs.data() + A.row_begin[ia];
                aend = A.pairs.data() + A.row_begin[ia + 1];
                nA_c = A.ent[ia].n_i;
                ++ia;
            }
            if (ib < B.colours.size() && B.colours[ib] == c) {
                bbeg = B.pairs.data() + B.row_begin[ib];
                bend = B.pairs.data() + B.row_begin[ib + 1];
                nB_c = B.ent[ib].n_i;
                ++ib;
            }
            // cross lists: (j in B, j > c) weighted nA_c * n_j^B, and (j in A, j > c) weighted nB_c * n_j^A
            size_t xa = cross && nB_c ? static_cast<size_t>(std::upper_bound(A.colours.begin(), A.colours.end(), c) -
                                                            A.colours.begin())
                                      : A.colours.size();
            size_t xb = cross && nA_c ? static_cast<size_t>(std::upper_bound(B.colours.begin(), B.colours.end(), c) -
                                                            B.colours.begin())
                                      : B.colours.size();
            // 4-way merge by j
            while (true) {
                int j = 1 << 30;
                if (abeg != aend) j = std::min(j, abeg->j);
                if (bbeg != bend) j = std::min(j, bbeg->j);
                if (xa < A.colours.size()) j = std::min(j, A.colours[xa]);
                if (xb < B.colours.size()) j = std::min(j, B.colours[xb]);
                if (j == (1 << 30)) break;
                PairEntry e;
                e.j = j;
                if (abeg != aend && abeg->j == j) {
                    e.n_ij += abeg->n_ij;
                    e.p_ij += abeg->p_ij;
                    ++abeg;
                }
                if (bbeg != bend && bbeg->j == j) {
                    e.n_ij += bbeg->n_ij;
                    e.p_ij += bbeg->p_ij;
                    ++bbeg;
                }
                if (xa < A.colours.size() && A.colours[xa] == j) {
                    e.n_ij += nB_c * A.ent[xa].n_i;
                    ++xa;
                }
                if (xb < B.colours.size() && B.colours[xb] == j) {
                    e.n_ij += nA_c * B.ent[xb].n_i;
                    ++xb;
                }
                out.pairs.push_back(e);
            }
        }
        out.row_begin.push_back(static_cast<int>(out.pairs.size()));
    }
}

}  // namespace detail

// ---------------------------------------------------------------------------
// L (leaf) component: Fig. 12 row 1.  A zero-leaf of weight w has n_0 = w.
// ---------------------------------------------------------------------------
template <class F>
inline void compose_leaf(CountersT<F>& out, int colour, i64 weight) {
    out.reset_keep_capacity();
    if (colour == 0) {
        out.t.n0 = weight;
    } else {
        out.colours.push_back(colour);
        out.ent.emplace_back();
        out.ent.back().n_i = weight;
    }
    if constexpr (F::kQuartet) {
        out.row_begin.assign(out.colours.size() + 1, 0);
    }
    detail::finalize_box_sums(out);
}

// L→C is the identity on the counters (Fig. 12 column L→C).
template <class F>
inline void compose_L2C(CountersT<F>& out, const CountersT<F>& L) {
    out = L;
}

// ---------------------------------------------------------------------------
// C→G  (Fig. 12 / Fig. 15 column C→G): the C component becomes the single
// subtree below a new super root.  Counters "at the super root" vanish,
// parenthesised counters are filled from the bracket / plain counters of C,
// bracket counters are copied, ↑ counters are not defined (0).
// ---------------------------------------------------------------------------
template <class F>
inline void compose_C2G(CountersT<F>& out, const CountersT<F>& C) {
    out.reset_keep_capacity();
    out.t.n0 = C.t.n0;
    out.colours = C.colours;
    out.ent.resize(C.ent.size());
    const i64 nb = C.t.nb, n0 = C.t.n0;
    for (size_t k = 0; k < C.ent.size(); ++k) {
        const F& c = C.ent[k];
        F& o = out.ent[k];
        o.n_i = c.n_i;
        o.p_ii = c2(c.n_i);  // n_{(ii)} = C(n_i, 2)
        if constexpr (F::kQuartet) {
            o.p_0i = n0 * c.n_i;            // n_{(0i)} = n_0 n_i
            o.p_ib = c.n_i * (nb - c.n_i);  // n_{(i•)} = n_i (n_• − n_i)
            // Fig. 12 group 3 / Fig. 15 group B: n_{(x(yz))} = n_{[x(yz)]}^C
            o.p_0_pii = c.br_0_pii;
            o.p_b_pii = c.br_b_pii;
            o.p_i_p0b = c.br_i_p0b;
            o.p_i_pib = c.br_i_pib;
            o.p_i_p0i = c.br_i_p0i;
            o.p_0_pib = c.br_0_pib;
            o.p_b_p0i = c.br_b_p0i;
            // bracket counters: copied (C→G column "n^C")
            o.br_0_pii = c.br_0_pii;
            o.br_b_pii = c.br_b_pii;
            o.br_i_p0b = c.br_i_p0b;
            o.br_i_pib = c.br_i_pib;
            o.br_0_pib = c.br_0_pib;
            o.br_i_p0i = c.br_i_p0i;
            o.br_b_p0i = c.br_b_p0i;
        }
    }
    if constexpr (F::kQuartet) {
        // pair table: n_{(ij)} = n_i n_j for every pair (Fig. 15 group B,
        // PAPER-NOTE F15-a), n_{ij} = 0.
        out.row_begin.reserve(C.colours.size() + 1);
        for (size_t k = 0; k < C.colours.size(); ++k) {
            out.row_begin.push_back(static_cast<int>(out.pairs.size()));
            for (size_t l = k + 1; l < C.colours.size(); ++l) {
                PairEntry e;
                e.j = C.colours[l];
                e.p_ij = C.ent[k].n_i * C.ent[l].n_i;
                out.pairs.push_back(e);
            }
        }
        out.row_begin.push_back(static_cast<int>(out.pairs.size()));
    }
    detail::finalize_box_sums(out);
}

// ---------------------------------------------------------------------------
// IG→C  (Fig. 12 / Fig. 15 column IG→C): the I node becomes the (only) node
// of the external path; every counter is copied from G except the ↑
// counters, which are 0.  G components never have non-zero ↑ counters, so a
// plain copy implements the column.
// ---------------------------------------------------------------------------
template <class F>
inline void compose_IG2C(CountersT<F>& out, const CountersT<F>& G) {
    out = G;
}

// ---------------------------------------------------------------------------
// GG→G  (Fig. 12 / Fig. 15 column GG→G):  n = n^{G1} + n^{G2} + g(G1,G2) + g(G2,G1).
// ---------------------------------------------------------------------------
template <class F>
inline void compose_GG2G(CountersT<F>& out, const CountersT<F>& A, const CountersT<F>& B, Scratch& sc) {
    out.reset_keep_capacity();
    out.t.n0 = A.t.n0 + B.t.n0;
    const i64 nbA = A.t.nb, nbB = B.t.nb, n0A = A.t.n0, n0B = B.t.n0;
    detail::merge_colours(A, B, out, [&](int, const F& a, const F& b, F& o) {
        o.n_i = a.n_i + b.n_i;
        // n_{i•}: Sec. 7.1 worked example / Fig. 12 group 2
        o.n_ib = a.n_ib + b.n_ib + a.n_i * (nbB - b.n_i) + b.n_i * (nbA - a.n_i);
        o.p_ii = a.p_ii + b.p_ii;  // parenthesised: plain sum
        if constexpr (F::kQuartet) {
            // ---- Fig. 12 group 2 (sum + g + g) ----
            o.n_0i = a.n_0i + b.n_0i + n0A * b.n_i + n0B * a.n_i;
            o.n_ii = a.n_ii + b.n_ii + a.n_i * b.n_i;  // ½ n_i' n_i'' twice
            o.n_0_pii = a.n_0_pii + b.n_0_pii + n0A * b.p_ii + n0B * a.p_ii;
            o.n_b_pii = a.n_b_pii + b.n_b_pii + (nbA - a.n_i) * b.p_ii + (nbB - b.n_i) * a.p_ii;
            o.n_i_p0b = a.n_i_p0b + b.n_i_p0b + a.n_i * (B.t.p_0b - b.p_0i) + b.n_i * (A.t.p_0b - a.p_0i);
            o.n_i_pbb = a.n_i_pbb + b.n_i_pbb + a.n_i * (B.t.p_bb - b.p_ii) + b.n_i * (A.t.p_bb - a.p_ii);
            o.n_i_pbs = a.n_i_pbs + b.n_i_pbs + a.n_i * (B.t.p_bs - b.p_ib) + b.n_i * (A.t.p_bs - a.p_ib);
            // ---- Fig. 12 group 3: parenthesised, plain sums ----
            o.p_0i = a.p_0i + b.p_0i;
            o.p_ib = a.p_ib + b.p_ib;
            o.p_0_pii = a.p_0_pii + b.p_0_pii;
            o.p_b_pii = a.p_b_pii + b.p_b_pii;
            o.p_i_p0b = a.p_i_p0b + b.p_i_p0b;
            // ---- Fig. 12 group 4: brackets (sum + g + g) ----
            o.br_0_pii = a.br_0_pii + b.br_0_pii + n0A * b.p_ii + n0B * a.p_ii;
            o.br_b_pii = a.br_b_pii + b.br_b_pii + (nbA - a.n_i) * b.p_ii + (nbB - b.n_i) * a.p_ii;
            o.br_i_p0b = a.br_i_p0b + b.br_i_p0b + a.n_i * (B.t.p_0b - b.p_0i) + b.n_i * (A.t.p_0b - a.p_0i);
            // ---- Fig. 15 group A (sum + g + g) ----
            o.n_i_p0i = a.n_i_p0i + b.n_i_p0i + a.n_i * b.p_0i + b.n_i * a.p_0i;
            o.n_b_p0i = a.n_b_p0i + b.n_b_p0i + (nbA - a.n_i) * b.p_0i + (nbB - b.n_i) * a.p_0i;
            o.n_0_pib = a.n_0_pib + b.n_0_pib + n0A * b.p_ib + n0B * a.p_ib;
            o.n_i_pib = a.n_i_pib + b.n_i_pib + a.n_i * b.p_ib + b.n_i * a.p_ib;
            // n_{•(i•)} (same colour • at the node and in the pair, Fig. 14 picture): the printed g
            // (n_•' − n_i') n_{(i•)}'' is the sum n_{•(i•)} + n_{•(i□)}; the correct g is
            // Σ_{j≠i} n_j' n_{(ij)}'', accumulated from the pair table below (PAPER-NOTE F15-c).
            o.n_b_pib = a.n_b_pib + b.n_b_pib;
            // n_{•(i□)} ★: (n_•' − n_i') n_{(i•)}'' − Σ_{j≠i} n_j' n_{(ij)}'' (pair part added below)
            o.n_b_pis = a.n_b_pis + b.n_b_pis + (nbA - a.n_i) * b.p_ib + (nbB - b.n_i) * a.p_ib;
            // ---- Fig. 15 group B: parenthesised sums ----
            o.p_i_pib = a.p_i_pib + b.p_i_pib;
            o.p_i_p0i = a.p_i_p0i + b.p_i_p0i;
            o.p_0_pib = a.p_0_pib + b.p_0_pib;
            o.p_b_p0i = a.p_b_p0i + b.p_b_p0i;
            // ---- Fig. 15 group C: brackets (sum + g + g) ----
            o.br_i_pib = a.br_i_pib + b.br_i_pib + a.p_ib * b.n_i + b.p_ib * a.n_i;
            o.br_0_pib = a.br_0_pib + b.br_0_pib + n0A * b.p_ib + n0B * a.p_ib;
            o.br_i_p0i = a.br_i_p0i + b.br_i_p0i + a.n_i * b.p_0i + b.n_i * a.p_0i;
            o.br_b_p0i = a.br_b_p0i + b.br_b_p0i + a.p_0i * (nbB - b.n_i) + b.p_0i * (nbA - a.n_i);
            // ↑ counters are not defined for G (stay 0).
        }
    });
    if constexpr (F::kQuartet) {
        // pair-table part of n_{•(i□)} ★ (Fig. 15 group A):
        //   g(G',G'') = Σ_{j≠i} (n_•' − n_i' − n_j') n_{(ij)}''  — the
        //   "− n_j' n_{(ij)}''" part is accumulated here by iterating the
        //   pair table of G'' and looking up n_j' in scratch.
        auto sub_pairs = [&](const CountersT<F>& P /*G'*/, const CountersT<F>& Q /*G''*/) {
            for (size_t k = 0; k < P.colours.size(); ++k) {
                sc.ensure(P.colours[k]);
                sc.n1[P.colours[k]] = P.ent[k].n_i;
            }
            for (size_t k = 0; k < out.colours.size(); ++k) {
                sc.ensure(out.colours[k]);
                sc.pos[out.colours[k]] = static_cast<int>(k);
            }
            for (size_t k = 0; k < Q.colours.size(); ++k) {
                const int i = Q.colours[k];
                for (int r = Q.row_begin[k]; r < Q.row_begin[k + 1]; ++r) {
                    const int j = Q.pairs[r].j;
                    const i64 pij = Q.pairs[r].p_ij;
                    sc.ensure(j);
                    out.ent[sc.pos[i]].n_b_pis -= sc.n1[j] * pij;
                    out.ent[sc.pos[j]].n_b_pis -= sc.n1[i] * pij;
                    out.ent[sc.pos[i]].n_b_pib += sc.n1[j] * pij;
                    out.ent[sc.pos[j]].n_b_pib += sc.n1[i] * pij;
                }
            }
            for (size_t k = 0; k < P.colours.size(); ++k) sc.n1[P.colours[k]] = 0;
        };
        sub_pairs(A, B);
        sub_pairs(B, A);
        detail::merge_pairs(A, B, out, /*cross=*/true);
    }
    detail::finalize_box_sums(out);
}

// ---------------------------------------------------------------------------
// CC→C  (Fig. 12 / Fig. 15 column CC→C):  n = n^{C1} + n^{C2} + f(C1,C2),
// where C1 is the LOWER and C2 the UPPER component (paper_figures.md,
// "Notation").
// ---------------------------------------------------------------------------
template <class F>
inline void compose_CC2C(CountersT<F>& out, const CountersT<F>& A /*C1 lower*/, const CountersT<F>& B /*C2 upper*/,
                         Scratch& sc) {
    out.reset_keep_capacity();
    out.t.n0 = A.t.n0 + B.t.n0;
    const i64 nbA = A.t.nb, nbB = B.t.nb, n0A = A.t.n0, n0B = B.t.n0;
    (void)n0B;
    detail::merge_colours(A, B, out, [&](int, const F& a, const F& b, F& o) {
        o.n_i = a.n_i + b.n_i;
        o.n_ib = a.n_ib + b.n_ib;  // group 2: plain sum for CC→C
        o.p_ii = a.p_ii + b.p_ii;
        o.n_iUb = a.n_iUb + b.n_iUb + a.n_i * (nbB - b.n_i);  // Fig. 12 group 5
        if constexpr (F::kQuartet) {
            // ---- plain sums: Fig. 12 group 2-3, Fig. 15 group A-B ----
            o.n_0i = a.n_0i + b.n_0i;
            o.n_ii = a.n_ii + b.n_ii;
            o.p_0i = a.p_0i + b.p_0i;
            o.p_ib = a.p_ib + b.p_ib;
            o.n_0_pii = a.n_0_pii + b.n_0_pii;
            o.n_b_pii = a.n_b_pii + b.n_b_pii;
            o.n_i_p0b = a.n_i_p0b + b.n_i_p0b;
            o.n_i_pbb = a.n_i_pbb + b.n_i_pbb;
            o.n_i_pbs = a.n_i_pbs + b.n_i_pbs;
            o.n_i_p0i = a.n_i_p0i + b.n_i_p0i;
            o.n_b_p0i = a.n_b_p0i + b.n_b_p0i;
            o.n_0_pib = a.n_0_pib + b.n_0_pib;
            o.n_b_pis = a.n_b_pis + b.n_b_pis;
            o.n_i_pib = a.n_i_pib + b.n_i_pib;
            o.n_b_pib = a.n_b_pib + b.n_b_pib;
            o.p_0_pii = a.p_0_pii + b.p_0_pii;
            o.p_b_pii = a.p_b_pii + b.p_b_pii;
            o.p_i_p0b = a.p_i_p0b + b.p_i_p0b;
            o.p_i_pib = a.p_i_pib + b.p_i_pib;
            o.p_i_p0i = a.p_i_p0i + b.p_i_p0i;
            o.p_0_pib = a.p_0_pib + b.p_0_pib;
            o.p_b_p0i = a.p_b_p0i + b.p_b_p0i;
            // ---- Fig. 12 group 4: brackets, sum + f(C1,C2) ----
            o.br_0_pii = a.br_0_pii + b.br_0_pii + c2(a.n_i) * n0B + n0A * b.p_ii + a.n_i * b.n_iU0;
            o.br_b_pii = a.br_b_pii + b.br_b_pii + c2(a.n_i) * (nbB - b.n_i) + (nbA - a.n_i) * b.p_ii +
                         a.n_i * b.n_iUb;
            o.br_i_p0b = a.br_i_p0b + b.br_i_p0b + n0A * (nbA - a.n_i) * b.n_i + a.n_i * (B.t.p_0b - b.p_0i) +
                         n0A * b.n_bUi + (nbA - a.n_i) * b.n_0Ui;
            // ---- Fig. 15 group C: brackets, sum + f(C1,C2) (complete formulas, PAPER-NOTE F15-b) ----
            o.br_i_pib = a.br_i_pib + b.br_i_pib + a.n_i * (nbA - a.n_i) * b.n_i + a.n_i * b.n_bUi +
                         (nbA - a.n_i) * b.n_iUi + a.n_i * b.p_ib;
            o.br_0_pib = a.br_0_pib + b.br_0_pib + a.n_i * (nbA - a.n_i) * n0B + a.n_i * (B.t.n_bU0 - b.n_iU0) +
                         (nbA - a.n_i) * b.n_iU0 + n0A * b.p_ib;
            o.br_i_p0i = a.br_i_p0i + b.br_i_p0i + a.n_i * b.p_0i + a.n_i * b.n_0Ui + n0A * a.n_i * b.n_i +
                         n0A * b.n_iUi;
            o.br_b_p0i = a.br_b_p0i + b.br_b_p0i + a.n_i * n0A * (nbB - b.n_i) + (nbA - a.n_i) * b.p_0i +
                         n0A * b.n_iUb + a.n_i * (B.t.n_0Ub - b.n_0Ui);
            // ---- Fig. 12 group 5: ↑ counters, sum + f(C1,C2) ----
            o.n_0Ui = a.n_0Ui + b.n_0Ui + n0A * b.n_i;
            o.n_iU0 = a.n_iU0 + b.n_iU0 + a.n_i * n0B;
            o.n_iUi = a.n_iUi + b.n_iUi + a.n_i * b.n_i;
            o.n_bUi = a.n_bUi + b.n_bUi + (nbA - a.n_i) * b.n_i;
            o.n_0Upii = a.n_0Upii + b.n_0Upii + n0A * b.p_ii;
            o.n_iUp0b = a.n_iUp0b + b.n_iUp0b + a.n_i * (B.t.p_0b - b.p_0i);  // PAPER-NOTE F12-a
            o.n_iUpbs = a.n_iUpbs + b.n_iUpbs + a.n_i * (B.t.p_bs - b.p_ib);
            o.n_bUpii = a.n_bUpii + b.n_bUpii + (nbA - a.n_i) * b.p_ii;
            o.n_iUpbb = a.n_iUpbb + b.n_iUpbb + a.n_i * (B.t.p_bb - b.p_ii);  // PAPER-NOTE F11-a
            o.n_iU0b = a.n_iU0b + b.n_iU0b + a.n_i * (B.t.n_0b - b.n_0i);
            o.n_iUbb = a.n_iUbb + b.n_iUbb + a.n_i * (B.t.n_bb - b.n_ii);
            o.n_iUbs = a.n_iUbs + b.n_iUbs + a.n_i * (B.t.n_bs - b.n_ib);
            o.n_bUii = a.n_bUii + b.n_bUii + (nbA - a.n_i) * b.n_ii;
            o.n_piiU0 = a.n_piiU0 + b.n_piiU0 + a.p_ii * n0B;
            o.n_piiUb = a.n_piiUb + b.n_piiUb + a.p_ii * (nbB - b.n_i);
            o.n_pbbUi = a.n_pbbUi + b.n_pbbUi + (A.t.p_bb - a.p_ii) * b.n_i;
            o.n_0UiUi = a.n_0UiUi + b.n_0UiUi + n0A * b.n_iUi + a.n_0Ui * b.n_i;
            o.n_iUbU0 = a.n_iUbU0 + b.n_iUbU0 + a.n_i * (B.t.n_bU0 - b.n_iU0) + a.n_iUb * n0B;
            o.n_iU0Ub = a.n_iU0Ub + b.n_iU0Ub + a.n_i * (B.t.n_0Ub - b.n_0Ui) + a.n_iU0 * (nbB - b.n_i);
            o.n_bUiUi = a.n_bUiUi + b.n_bUiUi + (nbA - a.n_i) * b.n_iUi + a.n_bUi * b.n_i;
            // ---- Fig. 15 group D: ↑ counters, sum + f(C1,C2) ----
            o.n_0Uib = a.n_0Uib + b.n_0Uib + n0A * b.n_ib;
            o.n_bU0i = a.n_bU0i + b.n_bU0i + (nbA - a.n_i) * b.n_0i;
            o.n_iUib = a.n_iUib + b.n_iUib + a.n_i * b.n_ib;
            // Counters used by Fig. 13 / Fig. 16 but missing from Fig. 11 / 14 (PAPER-NOTE F11-a, F14-a):
            o.n_0Uii = a.n_0Uii + b.n_0Uii + n0A * b.n_ii;              // n_{0↑ii}  (Fig. 13 γβ)
            o.n_iU0i = a.n_iU0i + b.n_iU0i + a.n_i * b.n_0i;            // n_{i↑0i}  (Fig. 16 γβ)
            // n_{•↑i□} ★ (Fig. 16 ββ row 1): (n_•^{C1} − n_i^{C1}) n_{i•}^{C2} − Σ_{c≠i} n_c^{C1} n_{ic}^{C2}; pair part below
            o.n_bUis = a.n_bUis + b.n_bUis + (nbA - a.n_i) * b.n_ib;
            // n_{•↑i•} ★ and n_{•↑(i•)} ★ (Fig. 16 αα/βα/ββ, not listed in Fig. 14):
            //   f = Σ_{j≠i} n_j^{C1} n_{ij}^{C2}  and  Σ_{j≠i} n_j^{C1} n_{(ij)}^{C2}; pair parts below
            o.n_bUib = a.n_bUib + b.n_bUib;
            o.n_bUpib = a.n_bUpib + b.n_bUpib;
            o.n_pibU0 = a.n_pibU0 + b.n_pibU0 + a.p_ib * n0B;
            o.n_pibUi = a.n_pibUi + b.n_pibUi + a.p_ib * b.n_i;
            o.n_p0iUi = a.n_p0iUi + b.n_p0iUi + a.p_0i * b.n_i;
            o.n_p0iUb = a.n_p0iUb + b.n_p0iUb + a.p_0i * (nbB - b.n_i);
            o.n_p0bUi = a.n_p0bUi + b.n_p0bUi + (A.t.p_0b - a.p_0i) * b.n_i;  // PAPER-NOTE F11-a
            o.n_0Upib = a.n_0Upib + b.n_0Upib + n0A * b.p_ib;
            o.n_iUpib = a.n_iUpib + b.n_iUpib + a.n_i * b.p_ib;
            o.n_bUp0i = a.n_bUp0i + b.n_bUp0i + (nbA - a.n_i) * b.p_0i;
            o.n_iUp0i = a.n_iUp0i + b.n_iUp0i + a.n_i * b.p_0i;
            // n_{•↑(i□)} ★: Σ_{k≠i} (n_•^{C1} − n_k^{C1} − n_i^{C1}) n_{(ik)}^{C2}
            //   = (n_•^{C1} − n_i^{C1}) n_{(i•)}^{C2} − Σ_{k≠i} n_k^{C1} n_{(ik)}^{C2}; pair part below
            o.n_bUpis = a.n_bUpis + b.n_bUpis + (nbA - a.n_i) * b.p_ib;
            o.n_bUiU0 = a.n_bUiU0 + b.n_bUiU0 + a.n_bUi * n0B + (nbA - a.n_i) * b.n_iU0;
            o.n_iUiU0 = a.n_iUiU0 + b.n_iUiU0 + a.n_iUi * n0B + a.n_i * b.n_iU0;
            o.n_0UbUi = a.n_0UbUi + b.n_0UbUi + (A.t.n_0Ub - a.n_0Ui) * b.n_i + n0A * b.n_bUi;
            o.n_bU0Ui = a.n_bU0Ui + b.n_bU0Ui + (A.t.n_bU0 - a.n_iU0) * b.n_i + (nbA - a.n_i) * b.n_0Ui;
            o.n_iU0Ui = a.n_iU0Ui + b.n_iU0Ui + a.n_iU0 * b.n_i + a.n_i * b.n_0Ui;
            o.n_0UiUb = a.n_0UiUb + b.n_0UiUb + n0A * b.n_iUb + a.n_0Ui * (nbB - b.n_i);
            o.n_iUiUb = a.n_iUiUb + b.n_iUiUb + a.n_iUi * (nbB - b.n_i) + a.n_i * b.n_iUb;
            o.n_iUbUi = a.n_iUbUi + b.n_iUbUi + a.n_i * b.n_bUi + a.n_iUb * b.n_i;  // PAPER-NOTE F14-a (n_{•↑i} already excludes • = i)
        }
    });
    if constexpr (F::kQuartet) {
        // pair part of n_{•↑(i□)} ★ (Fig. 15 group D): − Σ_{k≠i} n_k^{C1} n_{(ik)}^{C2}
        for (size_t k = 0; k < A.colours.size(); ++k) {
            sc.ensure(A.colours[k]);
            sc.n1[A.colours[k]] = A.ent[k].n_i;
        }
        for (size_t k = 0; k < out.colours.size(); ++k) {
            sc.ensure(out.colours[k]);
            sc.pos[out.colours[k]] = static_cast<int>(k);
        }
        for (size_t k = 0; k < B.colours.size(); ++k) {
            const int i = B.colours[k];
            for (int r = B.row_begin[k]; r < B.row_begin[k + 1]; ++r) {
                const int j = B.pairs[r].j;
                const i64 pij = B.pairs[r].p_ij, nij = B.pairs[r].n_ij;
                sc.ensure(j);
                out.ent[sc.pos[i]].n_bUpis -= sc.n1[j] * pij;
                out.ent[sc.pos[j]].n_bUpis -= sc.n1[i] * pij;
                out.ent[sc.pos[i]].n_bUis -= sc.n1[j] * nij;
                out.ent[sc.pos[j]].n_bUis -= sc.n1[i] * nij;
                out.ent[sc.pos[i]].n_bUib += sc.n1[j] * nij;
                out.ent[sc.pos[j]].n_bUib += sc.n1[i] * nij;
                out.ent[sc.pos[i]].n_bUpib += sc.n1[j] * pij;
                out.ent[sc.pos[j]].n_bUpib += sc.n1[i] * pij;
            }
        }
        for (size_t k = 0; k < A.colours.size(); ++k) sc.n1[A.colours[k]] = 0;
        detail::merge_pairs(A, B, out, /*cross=*/false);
        // colour-independent ↑ counter composed directly (Fig. 12 group 5):
        out.t.n_0Ubb = A.t.n_0Ubb + B.t.n_0Ubb + n0A * B.t.n_bb;
    }
    detail::finalize_box_sums(out);
}

}  // namespace soda13
