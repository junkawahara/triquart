// soda13.pdf, Figure 13, all CC and GG rows, in printed case/row order.
// Section 6 / Figure 8 uses the opposite tree ordering for the upper triangle.
// The individual sums are kept separate to make every edge split auditable.
#include "quartet_contributions.hpp"

#include <initializer_list>

namespace soda2013 {
namespace {
using S = SharedScalar;
using I = SharedIndexed;
using T = QuartetType;

void terms(Count& target, std::initializer_list<Count> values) {
    for (Count value : values) target = add(target, value);
}
Count exclude(Count total, Count first, Count second) { return subtract(subtract(total, first), second); }

template<class Function>
void each_color(const QuartetSharedCounters& a, const QuartetSharedCounters& b, Function function) {
    const QuartetSharedColor zero;
    std::size_t i = 0, j = 0;
    while (i < a.colors.size() || j < b.colors.size()) {
        const Color color = j == b.colors.size() ? a.colors[i].color :
                            i == a.colors.size() ? b.colors[j].color :
                            std::min(a.colors[i].color, b.colors[j].color);
        const auto& x = i < a.colors.size() && a.colors[i].color == color ? a.colors[i++] : zero;
        const auto& y = j < b.colors.size() && b.colors[j].color == color ? b.colors[j++] : zero;
        function(x, y);
    }
}
}  // namespace

void QuartetCases::accumulate(const QuartetCases& other) {
    for (std::size_t i = 0; i < values.size(); ++i) values[i] = add(values[i], other.values[i]);
}
Count QuartetCases::total() const {
    Count result = 0;
    for (Count value : values) result = add(result, value);
    return result;
}
QuartetCases complete_quartet_cases(const QuartetCases& forward, const QuartetCases& reverse) {
    QuartetCases result;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            if (i == j && forward.values[3 * i + j] != reverse.values[3 * j + i])
                throw std::logic_error("forward and reverse quartet diagonal counts disagree");
            result.values[3 * i + j] = i >= j ? forward.values[3 * i + j] : reverse.values[3 * j + i];
        }
    return result;
}

QuartetCases shared_quartet_cc(const QuartetSharedCounters& a, const QuartetSharedCounters& b) {
    QuartetCases result;
    each_color(a, b, [&](const QuartetSharedColor& x, const QuartetSharedColor& y) {
        const Count xi = x[I::i], yi = y[I::i];
        const Count other_a = subtract(a[S::b], xi), other_b = subtract(b[S::b], yi);
        const Count pair_a = choose(xi, 2);
        // alpha-alpha: splits 1,2,5,6; 3,4.
        terms(result.at(T::alpha, T::alpha), {
            multiply(xi, y.i_up_p_bb), multiply(pair_a, subtract(b[S::p_bb], y[I::p_ii]))
        });
        // beta-alpha: splits 1,2; 3; 4; 5,6 (last row uses k as index).
        terms(result.at(T::beta, T::alpha), {
            multiply(xi, y[I::i_up_p_bc]), multiply(pair_a, subtract(b[S::p_bc], y[I::p_i_b])),
            multiply(subtract(a[S::any_bc], x[I::any_i_b]), y[I::p_ii]),
            multiply(xi, exclude(b[S::b_up_p_cc], y.i_up_p_bb, y[I::b_up_p_ii]))
        });
        // beta-beta, T2 coloring ((ii),j,k): splits 1,2; 3; 4,5.
        terms(result.at(T::beta, T::beta), {
            multiply(xi, y[I::i_up_bc]), multiply(pair_a, subtract(b[S::bc], y[I::i_b])),
            multiply(xi, exclude(b[S::c_p_bb], y[I::i_p_bb], y[I::b_p_ii]))
        });
        // beta-beta, T2 coloring ((jk),i,i): splits 1,2; 3; 4,5.
        terms(result.at(T::beta, T::beta), {
            multiply(xi, exclude(b.c_up_bb, y[I::i_up_bb], y[I::b_up_ii])),
            multiply(subtract(a[S::any_bc], x[I::any_i_b]), y[I::ii]),
            multiply(xi, y[I::i_p_bc])
        });
        // gamma-alpha: splits 1,2; 3; 4; 5; 6.
        terms(result.at(T::gamma, T::alpha), {
            multiply(xi, y[I::i_up_p_zero_b]), multiply(pair_a, subtract(b[S::p_zero_b], y[I::p_zero_i])),
            multiply(multiply(a[S::zero], other_a), y[I::p_ii]),
            multiply(other_a, y[I::zero_up_p_ii]), multiply(a[S::zero], y[I::b_up_p_ii])
        });
        // gamma-beta, T2 coloring ((ii),j,0): splits 3; 1,2; 5; 4.
        terms(result.at(T::gamma, T::beta), {
            multiply(pair_a, subtract(b[S::zero_b], y[I::zero_i])), multiply(xi, y[I::i_up_zero_b]),
            multiply(a[S::zero], y[I::b_p_ii]), multiply(other_a, y[I::zero_p_ii])
        });
        // gamma-beta, T2 coloring ((j0),i,i): splits 1; 2; 3; 4,5.
        terms(result.at(T::gamma, T::beta), {
            multiply(other_a, y.zero_up_ii),
            multiply(a[S::zero], y[I::b_up_ii]), multiply(multiply(a[S::zero], other_a), y[I::ii]),
            multiply(xi, y[I::i_p_zero_b])
        });
        // gamma-gamma, T2 coloring (((ii),j),0): splits 1,2; 3; 4; 5; 6.
        terms(result.at(T::gamma, T::gamma), {
            multiply(xi, y[I::i_up_b_up_zero]), multiply(pair_a, subtract(b[S::b_up_zero], y[I::i_up_zero])),
            multiply(other_a, y[I::p_ii_up_zero]), multiply(x[I::any_b_p_ii], b[S::zero]),
            multiply(a[S::zero], y[I::p_b_p_ii])
        });
        // gamma-gamma, T2 coloring (((ii),0),j): splits 1,2; 3; 4; 5; 6.
        terms(result.at(T::gamma, T::gamma), {
            multiply(xi, y[I::i_up_zero_up_b]), multiply(pair_a, subtract(b[S::zero_up_b], y[I::zero_up_i])),
            multiply(a[S::zero], y[I::p_ii_up_b]), multiply(x[I::any_zero_p_ii], other_b),
            multiply(other_a, y[I::p_zero_p_ii])
        });
        // gamma-gamma, T2 coloring (((0j),i),i): splits 1; 2; 3; 4; 5; 6.
        terms(result.at(T::gamma, T::gamma), {
            multiply(a[S::zero], y[I::b_up_i_up_i]), multiply(other_a, y[I::zero_up_i_up_i]),
            multiply(multiply(a[S::zero], other_a), y[I::i_up_i]), multiply(xi, y.p_zero_b_up_i),
            multiply(x[I::any_i_p_zero_b], yi), multiply(xi, y[I::p_i_p_zero_b])
        });
    });
    return result;
}

QuartetCases shared_quartet_gg(const QuartetSharedCounters& a, const QuartetSharedCounters& b) {
    QuartetCases result;
    each_color(a, b, [&](const QuartetSharedColor& x, const QuartetSharedColor& y) {
        // alpha-alpha has NO sigma: each unordered pair of colors is already
        // counted once with the color of the pair in G1 as the summation index.
        terms(result.at(T::alpha, T::alpha), {
            multiply(x[I::p_ii], subtract(b[S::p_bb], y[I::p_ii]))
        });
        const auto ordered = [&](const QuartetSharedCounters& first, const QuartetSharedCounters& second,
                                 const QuartetSharedColor& u, const QuartetSharedColor& v) {
            const Count ui = u[I::i], other_first = subtract(first[S::b], ui);
            terms(result.at(T::beta, T::alpha), {
                multiply(u[I::p_ii], subtract(second[S::p_bc], v[I::p_i_b]))
            });
            terms(result.at(T::beta, T::beta), {
                multiply(u[I::p_ii], subtract(second[S::bc], v[I::i_b])),
                multiply(ui, exclude(second[S::c_p_bb], v[I::i_p_bb], v[I::b_p_ii])),
                multiply(u[I::ii], subtract(second[S::p_bc], v[I::p_i_b])), multiply(ui, v[I::i_p_bc])
            });
            terms(result.at(T::gamma, T::alpha), {
                multiply(u[I::p_ii], subtract(second[S::p_zero_b], v[I::p_zero_i]))
            });
            terms(result.at(T::gamma, T::beta), {
                multiply(u[I::p_ii], subtract(second[S::zero_b], v[I::zero_i])),
                multiply(other_first, v[I::zero_p_ii]),
                multiply(u[I::ii], subtract(second[S::p_zero_b], v[I::p_zero_i])), multiply(ui, v[I::i_p_zero_b])
            });
            terms(result.at(T::gamma, T::gamma), {
                multiply(first[S::zero], v[I::p_b_p_ii]), multiply(other_first, v[I::p_zero_p_ii]),
                multiply(ui, v[I::p_i_p_zero_b])
            });
        };
        ordered(a, b, x, y);
        ordered(b, a, y, x);
    });
    // gamma-beta ((ii),j,0), partition 34:5; already a scalar sum over colors.
    terms(result.at(T::gamma, T::beta), {
        multiply(a[S::zero], b[S::c_p_bb]), multiply(b[S::zero], a[S::c_p_bb])
    });
    return result;
}

}  // namespace soda2013
