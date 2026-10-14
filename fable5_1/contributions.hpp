// contributions.hpp — the contribution of one composition node of the HDT to
// the global counts (Sec. 6, Sec. 7.2–7.4 of soda13.pdf).
//
// Triplets and quartets "can only be anchored at nodes corresponding to
// CC→C and GG→G component compositions" (Sec. 6), so these two node kinds
// are the only ones that contribute:
//
//   * Fig. 10: shared resolved triplets (A) and shared unresolved triplets (E);
//   * Fig. 13: shared resolved quartets (A), nine (T1,T2)-type cases
//     αα, βα, ββ, γα, γβ, γγ directly and αβ, αγ, βγ by swapping T1 and T2
//     (Fig. 8); the swapped run only counts the asymmetric cases βα, γα, γβ;
//   * Fig. 16: resolved quartets with different splits (B), same nine cases
//     (Fig. 9).
//
// Every term below is labelled with its figure row ("αα 1,2,5,6" etc.) so it
// can be compared with `paper_figures.md`.  In CC→C, A is the lower
// component C1 and B the upper component C2.
#pragma once

#include "count_type.hpp"
#include "counters.hpp"

namespace soda13 {

struct Contrib {
    count_t a = 0;  // triplets: A (shared resolved);     quartets: A (shared resolved)
    count_t b = 0;  // triplets: E (shared unresolved);   quartets: B (resolved, different split)
    Contrib& operator+=(const Contrib& o) {
        a += o.a;
        b += o.b;
        return *this;
    }
};

// Debug aid: restrict the quartet contributions to a subset of the nine
// (T1-type, T2-type) cases.  Bits: AA=1, BA=2, BB=4, GA=8, GB=16, GG=32.
// Default: everything.  Only consulted in quartet mode.
enum CaseBit : unsigned { CASE_AA = 1, CASE_BA = 2, CASE_BB = 4, CASE_GA = 8, CASE_GB = 16, CASE_GG = 32, CASE_ALL = 63 };
inline unsigned& contrib_case_mask() {
    static unsigned mask = CASE_ALL;
    return mask;
}

namespace detail {

// Σ_{i<j} n_i^{A} n_j^{A} · (field of the pair (i,j) in B's pair table).
template <class F, class Get>
inline count_t sum_pairs_A_weights(const CountersT<F>& A, const CountersT<F>& B, Scratch& sc, Get get) {
    count_t s = 0;
    for (size_t k = 0; k < A.colours.size(); ++k) {
        sc.ensure(A.colours[k]);
        sc.n1[A.colours[k]] = A.ent[k].n_i;
    }
    for (size_t k = 0; k < B.colours.size(); ++k) {
        const int i = B.colours[k];
        sc.ensure(i);
        for (int r = B.row_begin[k]; r < B.row_begin[k + 1]; ++r) {
            const int j = B.pairs[r].j;
            sc.ensure(j);
            s += static_cast<count_t>(sc.n1[i]) * sc.n1[j] * get(B.pairs[r]);
        }
    }
    for (size_t k = 0; k < A.colours.size(); ++k) sc.n1[A.colours[k]] = 0;
    return s;
}

// Σ_{i<j} getA(pair_A) · getB(pair_B) over pairs present in both tables.
template <class F, class GetA, class GetB>
inline count_t sum_pairs_both(const CountersT<F>& A, const CountersT<F>& B, GetA ga, GetB gb) {
    count_t s = 0;
    size_t ka = 0, kb = 0;
    while (ka < A.colours.size() && kb < B.colours.size()) {
        if (A.colours[ka] < B.colours[kb]) {
            ++ka;
            continue;
        }
        if (B.colours[kb] < A.colours[ka]) {
            ++kb;
            continue;
        }
        int ra = A.row_begin[ka], rae = A.row_begin[ka + 1];
        int rb = B.row_begin[kb], rbe = B.row_begin[kb + 1];
        while (ra < rae && rb < rbe) {
            if (A.pairs[ra].j < B.pairs[rb].j) {
                ++ra;
            } else if (B.pairs[rb].j < A.pairs[ra].j) {
                ++rb;
            } else {
                s += static_cast<count_t>(ga(A.pairs[ra])) * gb(B.pairs[rb]);
                ++ra;
                ++rb;
            }
        }
        ++ka;
        ++kb;
    }
    return s;
}

// Iterate the union of colours of A and B with zero entries for absences.
template <class F, class Fn>
inline void for_colours(const CountersT<F>& A, const CountersT<F>& B, Fn fn) {
    static const F zero{};
    size_t ia = 0, ib = 0;
    const size_t na = A.colours.size(), nbn = B.colours.size();
    while (ia < na || ib < nbn) {
        const F* a = &zero;
        const F* b = &zero;
        if (ib >= nbn || (ia < na && A.colours[ia] < B.colours[ib])) {
            a = &A.ent[ia++];
        } else if (ia >= na || B.colours[ib] < A.colours[ia]) {
            b = &B.ent[ib++];
        } else {
            a = &A.ent[ia++];
            b = &B.ent[ib++];
        }
        fn(*a, *b);
    }
}

struct CaseFlags {
    bool AA, BA, BB, GA, GB, GG;
    explicit CaseFlags(bool asym_only) {
        const unsigned m = contrib_case_mask();
        AA = !asym_only && (m & CASE_AA);
        BB = !asym_only && (m & CASE_BB);
        GG = !asym_only && (m & CASE_GG);
        BA = (m & CASE_BA) != 0;
        GA = (m & CASE_GA) != 0;
        GB = (m & CASE_GB) != 0;
    }
};

}  // namespace detail

// ---------------------------------------------------------------------------
// CC→C node.  A = C1 (lower), B = C2 (upper).
// ---------------------------------------------------------------------------
template <class F>
inline Contrib contrib_CC2C(const CountersT<F>& A, const CountersT<F>& B, bool asym_only, Scratch& sc) {
    Contrib r;
    const count_t nbA = A.t.nb, nbB = B.t.nb, n0A = A.t.n0, n0B = B.t.n0;
    if constexpr (!F::kQuartet) {
        (void)asym_only;
        (void)sc;
        (void)n0A;
        (void)n0B;
        // Fig. 10, CC→C column.
        detail::for_colours(A, B, [&](const F& a, const F& b) {
            const count_t ni = a.n_i;
            // resolved–resolved, edges 1,2 / 3 / 4
            r.a += ni * b.n_iUb;
            r.a += static_cast<count_t>(c2(a.n_i)) * (nbB - b.n_i);
            r.a += (nbA - a.n_i) * b.p_ii;
            // unresolved–unresolved, edges 1,2,3
            r.b += ni * (B.t.n_bs - b.n_ib);
        });
        return r;
    } else {
        const detail::CaseFlags cf(asym_only);
        count_t sa = 0, sb = 0;
        detail::for_colours(A, B, [&](const F& a, const F& b) {
            const count_t ni = a.n_i, nbi = nbA - a.n_i, c2i = c2(a.n_i);
            // ================= Fig. 13 (shared resolved quartets, A) =================
            if (cf.AA) {
                // αα 1,2,5,6 ; 3,4
                sa += ni * b.n_iUpbb;
                sa += c2i * (B.t.p_bb - b.p_ii);
            }
            if (cf.BA) {
                // βα: 1,2 ; 3 ; 4 ; 5,6
                sa += ni * b.n_iUpbs;
                sa += c2i * (B.t.p_bs - b.p_ib);
                sa += (A.t.br_bs - ni * nbi) * b.p_ii;
                sa += ni * (B.t.n_bUpss - b.n_iUpbb - b.n_bUpii);
            }
            if (cf.BB) {
                // ββ, T2 = ((i,i),j,k): 1,2 ; 3 ; 4,5
                sa += ni * b.n_iUbs;
                sa += c2i * (B.t.n_bs - b.n_ib);
                sa += ni * (B.t.n_s_pbb - b.n_i_pbb - b.n_b_pii);
                // ββ, T2 = ((j,k),i,i): 1,2 ; 3 ; 4,5
                sa += ni * (B.t.n_sUbb - b.n_iUbb - b.n_bUii);
                sa += (A.t.br_bs - ni * nbi) * b.n_ii;
                sa += ni * b.n_i_pbs;
            }
            if (cf.GA) {
                // γα: 1,2 ; 3 ; 4 ; 5 ; 6
                sa += ni * b.n_iUp0b;
                sa += c2i * (B.t.p_0b - b.p_0i);
                sa += n0A * nbi * b.p_ii;
                sa += nbi * b.n_0Upii;
                sa += n0A * b.n_bUpii;
            }
            if (cf.GB) {
                // γβ, T2 = ((i,i),j,0): 3 ; 1,2 ; 5 ; 4
                sa += c2i * (B.t.n_0b - b.n_0i);
                sa += ni * b.n_iU0b;
                sa += n0A * b.n_b_pii;
                sa += nbi * b.n_0_pii;
                // γβ, T2 = ((j,0),i,i): 1 ; 2 ; 3 ; 4,5
                sa += nbi * b.n_0Uii;
                sa += n0A * b.n_bUii;
                sa += n0A * nbi * b.n_ii;
                sa += ni * b.n_i_p0b;
            }
            if (cf.GG) {
                // γγ, T2 = (((i,i),j),0): 1,2 ; 3 ; 4 ; 5 ; 6
                sa += ni * b.n_iUbU0;
                sa += c2i * (B.t.n_bU0 - b.n_iU0);
                sa += nbi * b.n_piiU0;
                sa += a.br_b_pii * n0B;
                sa += n0A * b.p_b_pii;
                // γγ, T2 = (((i,i),0),j): 1,2 ; 3 ; 4 ; 5 ; 6
                sa += ni * b.n_iU0Ub;
                sa += c2i * (B.t.n_0Ub - b.n_0Ui);
                sa += n0A * b.n_piiUb;
                sa += a.br_0_pii * (nbB - b.n_i);
                sa += nbi * b.p_0_pii;
                // γγ, T2 = (((0,j),i),i): 1 ; 2 ; 3 ; 4 ; 5 ; 6
                sa += n0A * b.n_bUiUi;
                sa += nbi * b.n_0UiUi;
                sa += n0A * nbi * b.n_iUi;
                sa += ni * b.n_p0bUi;
                sa += static_cast<count_t>(a.br_i_p0b) * b.n_i;
                sa += ni * b.p_i_p0b;
            }
            // ================= Fig. 16 (disagreeing resolved quartets, B) =================
            if (cf.AA) {
                // αα 1,2,5,6   (3,4 ★ handled with the pair table below)
                sb += ni * b.n_bUpib;
            }
            if (cf.BA) {
                // βα: 1,5 ; 2,6 ; 3,4 ★(non-pair part)
                sb += ni * b.n_bUpis;
                sb += ni * (B.t.n_bUpbs - b.n_iUpib - b.n_bUpib);
                sb += ni * nbi * b.p_ib;
            }
            if (cf.BB) {
                // ββ: 1 ; 2 ; 3 ★(non-pair part) ; 4 ; 5
                sb += ni * b.n_bUis;
                sb += ni * (B.t.n_bUbs - b.n_bUib - b.n_iUib);
                sb += ni * nbi * b.n_ib;
                sb += ni * b.n_b_pis;
                sb += ni * (B.t.n_b_pbs - b.n_b_pib - b.n_i_pib);
            }
            if (cf.GA) {
                // γα: 1 ; 2 ; 3 ; 4 ; 5 ; 6
                sb += ni * b.n_0Upib;
                sb += n0A * b.n_iUpib;
                sb += n0A * ni * b.p_ib;
                sb += ni * nbi * b.p_0i;
                sb += ni * b.n_bUp0i;
                sb += nbi * b.n_iUp0i;
            }
            if (cf.GB) {
                // γβ, T2 = ((i,0),i,j): 1 ; 3 ; 4 ; 5   (2 colour independent)
                sb += ni * b.n_0Uib;
                sb += n0A * ni * b.n_ib;
                sb += ni * b.n_b_p0i;
                sb += nbi * b.n_i_p0i;
                // γβ, T2 = ((i,j),i,0): 1 ; 2 ; 3 ; 4   (5 colour independent)
                sb += ni * b.n_bU0i;
                sb += nbi * b.n_iU0i;
                sb += ni * nbi * b.n_0i;
                sb += ni * b.n_0_pib;
            }
            if (cf.GG) {
                // γγ (a) T2 = (((i,0),i),j): 1 ; 3 ; 4 ; 5 ; 6   (2 is colour independent)
                sb += ni * b.n_0UiUb;
                sb += n0A * ni * b.n_iUb;
                sb += ni * b.n_p0iUb;
                sb += a.br_i_p0i * (nbB - b.n_i);
                sb += nbi * b.p_i_p0i;
                // γγ (b) T2 = (((i,0),j),i): 1 ; 3 ; 4 ; 5 ; 6
                sb += ni * b.n_0UbUi;
                sb += n0A * ni * b.n_bUi;
                sb += nbi * b.n_p0iUi;
                sb += static_cast<count_t>(a.br_b_p0i) * b.n_i;
                sb += ni * b.p_b_p0i;
                // γγ (c) T2 = (((i,j),i),0): 1 ; 2 ; 3 ; 4
                sb += ni * b.n_bUiU0;
                sb += nbi * b.n_iUiU0;
                sb += ni * nbi * b.n_iU0;
                sb += ni * b.n_pibU0;
                // γγ (d) T2 = (((i,j),0),i): 1 ; 2 ; 3 ; 5 ; 6
                sb += ni * b.n_bU0Ui;
                sb += nbi * b.n_iU0Ui;
                sb += ni * nbi * b.n_0Ui;
                sb += static_cast<count_t>(a.br_0_pib) * b.n_i;
                sb += ni * b.p_0_pib;
            }
        });
        // colour-independent rows of Fig. 16
        if (cf.GG) {
            sb += n0A * B.t.n_bUbUs;   // γγ (a) 2
            sb += n0A * B.t.n_bUsUb;   // γγ (b) 2
            sb += n0B * A.t.br_b_pbs;  // γγ (c) 5
            sb += n0A * B.t.p_b_pbs;   // γγ (c) 6
            sb += n0A * B.t.n_pbsUb;   // γγ (d) 4
        }
        if (cf.GB) {
            sb += n0A * B.t.n_bUbs;   // γβ ((i,0),i,j) row 2
            sb += n0A * B.t.n_b_pbs;  // γβ ((i,j),i,0) row 5
        }
        // ★ pair-table rows of Fig. 16
        if (cf.AA || cf.BA) {
            const count_t pair_p = detail::sum_pairs_A_weights(A, B, sc, [](const PairEntry& e) { return e.p_ij; });
            if (cf.AA) sb += pair_p;       // αα 3,4
            if (cf.BA) sb -= 2 * pair_p;   // βα 3,4: Σ_{i≠j} n_i n_j (n_{(i•)} − n_{(ij)})
        }
        if (cf.BB) {
            const count_t pair_n = detail::sum_pairs_A_weights(A, B, sc, [](const PairEntry& e) { return e.n_ij; });
            sb -= 2 * pair_n;  // ββ 3: Σ_{i≠j} n_i n_j (n_{i•} − n_{ij})
        }
        r.a = sa;
        r.b = sb;
        return r;
    }
}

// ---------------------------------------------------------------------------
// GG→G node.  σ[f](G1,G2) = f(G1,G2) + f(G2,G1).
// ---------------------------------------------------------------------------
template <class F>
inline Contrib contrib_GG2G(const CountersT<F>& G1, const CountersT<F>& G2, bool asym_only, Scratch& sc) {
    (void)sc;
    Contrib r;
    if constexpr (!F::kQuartet) {
        (void)asym_only;
        auto h = [&](const CountersT<F>& P, const CountersT<F>& Q) {
            detail::for_colours(P, Q, [&](const F& p, const F& q) {
                r.a += static_cast<count_t>(p.p_ii) * (Q.t.nb - q.n_i);   // Fig. 10 resolved
                r.b += static_cast<count_t>(p.n_i) * (Q.t.n_bs - q.n_ib);  // Fig. 10 unresolved
            });
        };
        h(G1, G2);
        h(G2, G1);
        return r;
    } else {
        const detail::CaseFlags cf(asym_only);
        count_t sa = 0, sb = 0;
        // σ rows: evaluated for (P,Q) = (G1,G2) and (G2,G1)
        auto h = [&](const CountersT<F>& P, const CountersT<F>& Q) {
            const count_t nbP = P.t.nb, n0P = P.t.n0;
            detail::for_colours(P, Q, [&](const F& pp, const F& q) {
                // 128-bit copies of the G' factors so that every product is wide
                struct {
                    count_t p_ii, n_ii, n_i, p_ib, p_0i, n_b_pis, n_i_p0i, n_b_p0i, n_0_pib;
                } p{pp.p_ii, pp.n_ii, pp.n_i, pp.p_ib, pp.p_0i, pp.n_b_pis, pp.n_i_p0i, pp.n_b_p0i, pp.n_0_pib};
                const count_t nbi = nbP - p.n_i;
                // ---------- Fig. 13 ----------
                if (cf.BA) sa += p.p_ii * (Q.t.p_bs - q.p_ib);  // βα 3:4
                if (cf.BB) {
                    // ββ ((i,i),j,k): 3:45 ; 34:5,35:4
                    sa += p.p_ii * (Q.t.n_bs - q.n_ib);
                    sa += p.n_i * (Q.t.n_s_pbb - q.n_i_pbb - q.n_b_pii);
                    // ββ ((j,k),i,i): 3:45 ; 34:5,35:4
                    sa += p.n_ii * (Q.t.p_bs - q.p_ib);
                    sa += p.n_i * q.n_i_pbs;
                }
                if (cf.GA) sa += p.p_ii * (Q.t.p_0b - q.p_0i);  // γα 3:4
                if (cf.GB) {
                    // γβ ((i,i),j,0): 3:45 ; 35:4   (34:5 colour independent)
                    sa += p.p_ii * (Q.t.n_0b - q.n_0i);
                    sa += nbi * q.n_0_pii;
                    // γβ ((j,0),i,i): 3:45 ; 34:5,35:4
                    sa += p.n_ii * (Q.t.p_0b - q.p_0i);
                    sa += p.n_i * q.n_i_p0b;
                }
                if (cf.GG) {
                    // γγ 5:6 for the three sub-cases
                    sa += n0P * q.p_b_pii;
                    sa += nbi * q.p_0_pii;
                    sa += p.n_i * q.p_i_p0b;
                }
                // ---------- Fig. 16 ----------
                if (cf.BB) {
                    // ββ: 3:45 ★ non-pair part ; 35:4 ; 34:5
                    sb += p.p_ib * q.n_ib;
                    sb += p.n_b_pis * q.n_i;
                    sb += p.n_i * (Q.t.n_b_pbs - q.n_i_pib - q.n_b_pib);
                }
                if (cf.GA) sb += p.p_0i * q.p_ib;  // γα 3:4
                if (cf.GB) {
                    // γβ ((i,0),i,j): 3:45 ; 34:5 ; 35:4
                    sb += p.p_0i * q.n_ib;
                    sb += p.n_i_p0i * (Q.t.nb - q.n_i);
                    sb += p.n_b_p0i * q.n_i;
                    // γβ ((i,j),i,0): 3:45 ; 35:4   (34:5 colour independent)
                    sb += p.p_ib * q.n_0i;
                    sb += p.n_0_pib * q.n_i;
                }
                if (cf.GG) {
                    // γγ (a),(b),(d) 5:6
                    sb += nbi * q.p_i_p0i;
                    sb += p.n_i * q.p_b_p0i;
                    sb += p.n_i * q.p_0_pib;
                }
            });
            // colour-independent σ rows
            if (cf.GB) {
                sa += n0P * Q.t.n_s_pbb;     // Fig. 13 γβ ((i,i),j,0) 34:5
                sb += static_cast<count_t>(P.t.n_b_pbs) * Q.t.n0;  // Fig. 16 γβ ((i,j),i,0) 34:5
            }
            if (cf.GG) sb += n0P * Q.t.p_b_pbs;  // Fig. 16 γγ (c) 5:6
            if (cf.BB) {
                // Fig. 16 ββ 3:45 ★ pair part: − 2 Σ_{i<j} n_{(ij)}' n_{ij}''
                sb -= 2 * detail::sum_pairs_both(P, Q, [](const PairEntry& e) { return e.p_ij; },
                                                 [](const PairEntry& e) { return e.n_ij; });
            }
        };
        h(G1, G2);
        h(G2, G1);
        // rows without σ
        if (cf.AA) {
            // Fig. 13 αα 3:4: Σ_i n_{(ii)}^{G1} (n_{(••)}^{G2} − n_{(ii)}^{G2})
            detail::for_colours(G1, G2, [&](const F& p, const F& q) { sa += static_cast<count_t>(p.p_ii) * (G2.t.p_bb - q.p_ii); });
            // Fig. 16 αα 3:4 ★: Σ_{i<j} n_{(ij)}^{G1} n_{(ij)}^{G2}
            sb += detail::sum_pairs_both(G1, G2, [](const PairEntry& e) { return e.p_ij; },
                                         [](const PairEntry& e) { return e.p_ij; });
        }
        if (cf.BA) {
            // Fig. 16 βα 3:4 ★: Σ_i Σ_{j≠i} n_{(ij)}^{G1} (n_{(i•)}^{G2} − n_{(ij)}^{G2})
            //   = Σ_i n_{(i•)}^{G1} n_{(i•)}^{G2} − 2 Σ_{i<j} n_{(ij)}^{G1} n_{(ij)}^{G2}
            detail::for_colours(G1, G2, [&](const F& p, const F& q) { sb += static_cast<count_t>(p.p_ib) * q.p_ib; });
            sb -= 2 * detail::sum_pairs_both(G1, G2, [](const PairEntry& e) { return e.p_ij; },
                                             [](const PairEntry& e) { return e.p_ij; });
        }
        r.a = sa;
        r.b = sb;
        return r;
    }
}

}  // namespace soda13
