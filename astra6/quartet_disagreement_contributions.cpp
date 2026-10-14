// soda13.pdf, Figure 16: each numbered CC cut and GG partition is retained in
// printed order. Sums run over present colors only (Section 7); the starred
// two-color sums use the symmetric pair matrix in O((1+x)^2) time.
#include "quartet_disagreement_contributions.hpp"

#include <initializer_list>

namespace soda2013 {
namespace {
using S = SharedScalar;
using I = SharedIndexed;
using D = DisagreementIndexed;
using R = DisagreementScalar;
using T = QuartetType;

void terms(Count& target, std::initializer_list<Count> values) {
    for (Count value : values) target = add(target, value);
}
Count exclude(Count total, Count first, Count second) { return subtract(subtract(total, first), second); }

struct Palette {
    const QuartetCounters& a;
    const QuartetCounters& b;
    std::vector<std::array<std::size_t, 2>> indices;
    static constexpr std::size_t absent = static_cast<std::size_t>(-1);
    const QuartetSharedColor empty_shared{};
    const DisagreementColor empty_disagreement{};

    Palette(const QuartetCounters& first, const QuartetCounters& second) : a(first), b(second) {
        std::size_t i = 0, j = 0;
        while (i < a.shared.colors.size() || j < b.shared.colors.size()) {
            const Color color = j == b.shared.colors.size() ? a.shared.colors[i].color :
                i == a.shared.colors.size() ? b.shared.colors[j].color :
                std::min(a.shared.colors[i].color, b.shared.colors[j].color);
            const auto x = i < a.shared.colors.size() && a.shared.colors[i].color == color ? i++ : absent;
            const auto y = j < b.shared.colors.size() && b.shared.colors[j].color == color ? j++ : absent;
            indices.push_back({x, y});
        }
    }
    const QuartetSharedColor& shared(std::size_t side, std::size_t color) const {
        const auto index = indices[color][side];
        return index == absent ? empty_shared : (side ? b : a).shared.colors[index];
    }
    const DisagreementColor& disagreement(std::size_t side, std::size_t color) const {
        const auto index = indices[color][side];
        return index == absent ? empty_disagreement : (side ? b : a).colors[index];
    }
    DisagreementPair pair(std::size_t side, std::size_t i, std::size_t j) const {
        const auto x = indices[i][side], y = indices[j][side];
        return x == absent || y == absent ? DisagreementPair{} : (side ? b : a).pair(x, y);
    }
};
}  // namespace

QuartetCases disagreement_quartet_cc(const QuartetCounters& a, const QuartetCounters& b) {
    const Palette palette(a, b);
    const auto& sa = a.shared;
    const auto& sb = b.shared;
    QuartetCases result;
    for (std::size_t i = 0; i < palette.indices.size(); ++i) {
        const auto& x = palette.shared(0, i);
        const auto& y = palette.shared(1, i);
        const auto& dx = palette.disagreement(0, i);
        const auto& dy = palette.disagreement(1, i);
        const Count xi = x[I::i], yi = y[I::i];
        const Count other_a = subtract(sa[S::b], xi), other_b = subtract(sb[S::b], yi);
        // alpha-alpha: cuts 1,2,5,6; starred 3,4 below.
        terms(result.at(T::alpha, T::alpha), {multiply(xi, dy.b_up_p_i_b)});
        // beta-alpha: 1,5; 2,6; starred 3,4 below.
        terms(result.at(T::beta, T::alpha), {
            multiply(xi, dy[D::b_up_p_i_c]),
            multiply(xi, exclude(b.b_up_p_b_c, dy[D::i_up_p_i_b], dy.b_up_p_i_b))
        });
        // beta-beta: 1; 2; starred 3 below; 4; 5.
        terms(result.at(T::beta, T::beta), {
            multiply(xi, dy.b_up_i_c),
            multiply(xi, exclude(b[R::b_up_b_c], dy.b_up_i_b, dy[D::i_up_i_b])),
            multiply(xi, dy[D::b_p_i_c]),
            multiply(xi, exclude(b[R::b_p_b_c], dy[D::b_p_i_b], dy[D::i_p_i_b]))
        });
        // gamma-alpha: 1; 2; 3; 4; 5; 6.
        terms(result.at(T::gamma, T::alpha), {
            multiply(xi, dy[D::zero_up_p_i_b]), multiply(sa[S::zero], dy[D::i_up_p_i_b]),
            multiply(multiply(sa[S::zero], xi), y[I::p_i_b]),
            multiply(multiply(xi, other_a), y[I::p_zero_i]),
            multiply(xi, dy[D::b_up_p_zero_i]), multiply(other_a, dy[D::i_up_p_zero_i])
        });
        // gamma-beta, ((i0),i,j): 1; scalar 2 below; 3; 4; 5.
        terms(result.at(T::gamma, T::beta), {
            multiply(xi, dy[D::zero_up_i_b]), multiply(multiply(sa[S::zero], xi), y[I::i_b]),
            multiply(xi, dy[D::b_p_zero_i]), multiply(other_a, dy[D::i_p_zero_i])
        });
        // gamma-beta, ((ij),i,0): 1; 2; 3; 4; scalar 5 below.
        terms(result.at(T::gamma, T::beta), {
            multiply(xi, dy[D::b_up_zero_i]), multiply(other_a, dy.i_up_zero_i),
            multiply(multiply(xi, other_a), y[I::zero_i]), multiply(xi, dy[D::zero_p_i_b])
        });
        // gamma-gamma, (((i0),i),j): 1; scalar 2 below; 3; 4; 5; 6.
        terms(result.at(T::gamma, T::gamma), {
            multiply(xi, dy[D::zero_up_i_up_b]), multiply(multiply(sa[S::zero], xi), y[I::i_up_b]),
            multiply(xi, dy[D::p_zero_i_up_b]), multiply(dx[D::any_i_p_zero_i], other_b),
            multiply(other_a, dy[D::p_i_p_zero_i])
        });
        // gamma-gamma, (((i0),j),i): 1; scalar 2 below; 3; 4; 5; 6.
        terms(result.at(T::gamma, T::gamma), {
            multiply(xi, dy[D::zero_up_b_up_i]), multiply(multiply(sa[S::zero], xi), y[I::b_up_i]),
            multiply(other_a, dy[D::p_zero_i_up_i]), multiply(dx[D::any_b_p_zero_i], yi),
            multiply(xi, dy[D::p_b_p_zero_i])
        });
        // gamma-gamma, (((ij),i),0): 1; 2; 3; 4; 5; scalar 6 below.
        terms(result.at(T::gamma, T::gamma), {
            multiply(xi, dy[D::b_up_i_up_zero]), multiply(other_a, dy[D::i_up_i_up_zero]),
            multiply(multiply(xi, other_a), y[I::i_up_zero]), multiply(xi, dy[D::p_i_b_up_zero]),
            multiply(sb[S::zero], dx[D::any_i_p_i_b])
        });
        // gamma-gamma, (((ij),0),i): 1; 2; 3; scalar 4 below; 5; 6.
        terms(result.at(T::gamma, T::gamma), {
            multiply(xi, dy[D::b_up_zero_up_i]), multiply(other_a, dy[D::i_up_zero_up_i]),
            multiply(multiply(xi, other_a), y[I::zero_up_i]), multiply(dx[D::any_zero_p_i_b], yi),
            multiply(xi, dy[D::p_zero_p_i_b])
        });

        for (std::size_t j = 0; j < palette.indices.size(); ++j) {
            if (i == j) continue;
            const Count leaves = multiply(xi, palette.shared(0, j)[I::i]);
            const auto pair = palette.pair(1, i, j);
            if (i < j) terms(result.at(T::alpha, T::alpha), {multiply(leaves, pair.inside)});
            terms(result.at(T::beta, T::alpha), {multiply(leaves, subtract(y[I::p_i_b], pair.inside))});
            terms(result.at(T::beta, T::beta), {multiply(leaves, subtract(y[I::i_b], pair.at))});
        }
    }
    // Scalar rows already sum over the repeated color. In gamma-beta ((ij),i,0)
    // cut 5, the printed extra sum_i would count every quartet d times; the
    // Figure 14 definition of n_b(bc) and the adjacent GG cell give one product.
    terms(result.at(T::gamma, T::beta), {
        multiply(sa[S::zero], b[R::b_up_b_c]), multiply(sa[S::zero], b[R::b_p_b_c])
    });
    terms(result.at(T::gamma, T::gamma), {
        multiply(sa[S::zero], b[R::b_up_b_up_c]), multiply(sa[S::zero], b.b_up_c_up_b),
        multiply(sa[S::zero], b[R::p_b_p_b_c]), multiply(sa[S::zero], b[R::p_b_c_up_b])
    });
    return result;
}

QuartetCases disagreement_quartet_gg(const QuartetCounters& a, const QuartetCounters& b) {
    const Palette palette(a, b);
    QuartetCases result;
    for (std::size_t i = 0; i < palette.indices.size(); ++i) {
        const auto& x = palette.shared(0, i);
        const auto& y = palette.shared(1, i);
        const auto& dx = palette.disagreement(0, i);
        const auto& dy = palette.disagreement(1, i);
        for (std::size_t j = 0; j < palette.indices.size(); ++j) {
            if (i == j) continue;
            const auto ap = palette.pair(0, i, j), bp = palette.pair(1, i, j);
            // alpha-alpha / beta-alpha 3:4 have NO sigma in Figure 16.
            if (i < j) terms(result.at(T::alpha, T::alpha), {multiply(ap.inside, bp.inside)});
            terms(result.at(T::beta, T::alpha), {multiply(ap.inside, subtract(y[I::p_i_b], bp.inside))});
            // beta-beta 3:45 uses both argument orders.
            terms(result.at(T::beta, T::beta), {
                multiply(ap.inside, subtract(y[I::i_b], bp.at)),
                multiply(bp.inside, subtract(x[I::i_b], ap.at))
            });
        }
        const auto ordered = [&](const QuartetCounters& first, const QuartetCounters& second,
                                 const QuartetSharedColor& u, const QuartetSharedColor& v,
                                 const DisagreementColor& du, const DisagreementColor& dv) {
            // beta-beta 35:4; 34:5.
            terms(result.at(T::beta, T::beta), {
                multiply(du[D::b_p_i_c], v[I::i]),
                multiply(u[I::i], exclude(second[R::b_p_b_c], dv[D::i_p_i_b], dv[D::b_p_i_b]))
            });
            terms(result.at(T::gamma, T::alpha), {multiply(u[I::p_zero_i], v[I::p_i_b])});
            // gamma-beta ((i0),i,j): 3:45; 34:5; 35:4.
            terms(result.at(T::gamma, T::beta), {
                multiply(u[I::p_zero_i], v[I::i_b]),
                multiply(du[D::i_p_zero_i], subtract(second.shared[S::b], v[I::i])),
                multiply(du[D::b_p_zero_i], v[I::i])
            });
            // gamma-beta ((ij),i,0): 3:45; scalar 34:5 below; 35:4.
            terms(result.at(T::gamma, T::beta), {
                multiply(u[I::p_i_b], v[I::zero_i]), multiply(du[D::zero_p_i_b], v[I::i])
            });
            // gamma-gamma: the four T2 colorings, partition 5:6 in each.
            terms(result.at(T::gamma, T::gamma), {
                multiply(subtract(first.shared[S::b], u[I::i]), dv[D::p_i_p_zero_i]),
                multiply(u[I::i], dv[D::p_b_p_zero_i]), multiply(u[I::i], dv[D::p_zero_p_i_b])
            });
        };
        ordered(a, b, x, y, dx, dy);
        ordered(b, a, y, x, dy, dx);
    }
    terms(result.at(T::gamma, T::beta), {
        multiply(a[R::b_p_b_c], b.shared[S::zero]), multiply(b[R::b_p_b_c], a.shared[S::zero])
    });
    terms(result.at(T::gamma, T::gamma), {
        multiply(a.shared[S::zero], b[R::p_b_p_b_c]), multiply(b.shared[S::zero], a[R::p_b_p_b_c])
    });
    return result;
}

}  // namespace soda2013
