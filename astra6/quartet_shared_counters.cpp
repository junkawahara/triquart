// soda13.pdf, Figure 12: complete shared-counter composition table. Every row
// is present, including counters unnecessary for triplets. Three apparent
// printed errors are resolved from Section 7 / Figure 11, as annotated below
// and checked by the independent semantic oracle (try_and_error.md, Trial 6).
#include "quartet_shared_counters.hpp"

namespace soda2013 {
namespace {

using S = SharedScalar;
using I = SharedIndexed;

Count half(Count n) {
    if (n % 2) throw std::logic_error("nonintegral symmetric shared counter");
    return n / 2;
}

void derived(QuartetSharedCounters& result) {
    // Figure 12, boxed sums. Never iterate over absent colors.
    result[S::b] = result[S::bb] = result[S::bc] = 0;
    result[S::p_bb] = result[S::p_bc] = result[S::c_p_bb] = 0;
    result[S::b_up_p_cc] = result[S::any_bc] = 0;
    result.c_up_bb = 0;
    for (auto& row : result.colors) result[S::b] = add(result[S::b], row[I::i]);
    for (auto& row : result.colors) {
        result[S::bb] = add(result[S::bb], row[I::ii]);
        result[S::bc] = add(result[S::bc], row[I::i_b]);
        result[S::p_bb] = add(result[S::p_bb], row[I::p_ii]);
        result[S::p_bc] = add(result[S::p_bc], row[I::p_i_b]);
        result[S::c_p_bb] = add(result[S::c_p_bb], row[I::i_p_bb]);
        // The printed RHS omits the upward arrow. Figure 11 requires the
        // singleton to branch BELOW the off-path equal-color pair.
        result[S::b_up_p_cc] = add(result[S::b_up_p_cc], row[I::b_up_p_ii]);
        result.c_up_bb = add(result.c_up_bb, row[I::b_up_ii]);
        row[I::any_i_b] = multiply(row[I::i], subtract(result[S::b], row[I::i]));
        // The printed RHS omits brackets on n_i•. Section 7 defines brackets
        // as having NO external-path restriction, so sum n_[i•] here.
        result[S::any_bc] = add(result[S::any_bc], row[I::any_i_b]);
    }
    result[S::bc] = half(result[S::bc]);
    result[S::p_bc] = half(result[S::p_bc]);
    result[S::any_bc] = half(result[S::any_bc]);
}

void clear_path(QuartetSharedCounters& result) {
    result[S::zero_up_b] = result[S::b_up_zero] = result[S::zero_up_bb] = 0;
    for (auto& row : result.colors) {
        for (std::size_t i = static_cast<std::size_t>(I::zero_up_i);
             i <= static_cast<std::size_t>(I::p_bb_up_i); ++i) row.values[i] = 0;
        row.i_up_p_bb = row.p_zero_b_up_i = row.zero_up_ii = 0;
    }
}

// Common CC/GG initialization and the linear merge of sorted color tables.
template<class Function>
QuartetSharedCounters merge(const QuartetSharedCounters& a, const QuartetSharedCounters& b, Function function) {
    QuartetSharedCounters result;
    for (std::size_t i = 0; i < shared_scalar_count; ++i)
        result.values[i] = add(a.values[i], b.values[i]);
    const QuartetSharedColor zero;
    std::size_t i = 0, j = 0;
    result.colors.reserve(a.colors.size() + b.colors.size());
    while (i < a.colors.size() || j < b.colors.size()) {
        const Color color = j == b.colors.size() ? a.colors[i].color :
                            i == a.colors.size() ? b.colors[j].color :
                            std::min(a.colors[i].color, b.colors[j].color);
        const auto& x = i < a.colors.size() && a.colors[i].color == color ? a.colors[i++] : zero;
        const auto& y = j < b.colors.size() && b.colors[j].color == color ? b.colors[j++] : zero;
        QuartetSharedColor row;
        row.color = color;
        for (std::size_t c = 0; c < shared_indexed_count; ++c) row.values[c] = add(x.values[c], y.values[c]);
        row.i_up_p_bb = add(x.i_up_p_bb, y.i_up_p_bb);
        row.p_zero_b_up_i = add(x.p_zero_b_up_i, y.p_zero_b_up_i);
        row.zero_up_ii = add(x.zero_up_ii, y.zero_up_ii);
        function(row, x, y);
        result.colors.push_back(row);
    }
    return result;
}

void increment(QuartetSharedColor& row, I counter, Count extra) { row[counter] = add(row[counter], extra); }

}  // namespace

QuartetSharedCounters quartet_shared_leaf(Color color, Count zero_multiplicity) {
    QuartetSharedCounters result;
    if (!color) result[S::zero] = zero_multiplicity;
    else {
        if (zero_multiplicity != 1) throw std::invalid_argument("a positive leaf has multiplicity one");
        QuartetSharedColor row;
        row.color = color;
        row[I::i] = 1;
        result.colors.push_back(row);
        result[S::b] = 1;
    }
    return result;
}

QuartetSharedCounters quartet_shared_cc(const QuartetSharedCounters& a, const QuartetSharedCounters& b) {
    auto result = merge(a, b, [&](QuartetSharedColor& row, const QuartetSharedColor& x, const QuartetSharedColor& y) {
        const Count other_a = subtract(a[S::b], x[I::i]);
        const Count other_b = subtract(b[S::b], y[I::i]);
        // Figure 12, bracketed triplets: additive terms plus f(C1,C2).
        increment(row, I::any_zero_p_ii,
            add(add(multiply(choose(x[I::i], 2), b[S::zero]), multiply(a[S::zero], y[I::p_ii])),
                multiply(x[I::i], y[I::i_up_zero])));
        increment(row, I::any_b_p_ii,
            add(add(multiply(choose(x[I::i], 2), other_b), multiply(other_a, y[I::p_ii])),
                multiply(x[I::i], y[I::i_up_b])));
        increment(row, I::any_i_p_zero_b,
            add(add(multiply(multiply(a[S::zero], other_a), y[I::i]),
                    multiply(x[I::i], subtract(b[S::p_zero_b], y[I::p_zero_i]))),
                add(multiply(a[S::zero], y[I::b_up_i]), multiply(other_a, y[I::zero_up_i]))));

        // Figure 12, C-only path rows, in their printed order.
        increment(row, I::zero_up_i, multiply(a[S::zero], y[I::i]));
        increment(row, I::i_up_zero, multiply(x[I::i], b[S::zero]));
        increment(row, I::i_up_i, multiply(x[I::i], y[I::i]));
        increment(row, I::i_up_b, multiply(x[I::i], other_b));
        increment(row, I::b_up_i, multiply(other_a, y[I::i]));
        increment(row, I::zero_up_p_ii, multiply(a[S::zero], y[I::p_ii]));
        // Printed row i↑(0•) duplicates the 0↑(ii) RHS. From Figure 11,
        // choose a lower i and an upper off-path (0,j) pair with j != i.
        increment(row, I::i_up_p_zero_b, multiply(x[I::i], subtract(b[S::p_zero_b], y[I::p_zero_i])));
        increment(row, I::i_up_p_bc, multiply(x[I::i], subtract(b[S::p_bc], y[I::p_i_b])));
        increment(row, I::b_up_p_ii, multiply(other_a, y[I::p_ii]));
        increment(row, I::i_up_zero_b, multiply(x[I::i], subtract(b[S::zero_b], y[I::zero_i])));
        increment(row, I::i_up_bb, multiply(x[I::i], subtract(b[S::bb], y[I::ii])));
        increment(row, I::i_up_bc, multiply(x[I::i], subtract(b[S::bc], y[I::i_b])));
        increment(row, I::b_up_ii, multiply(other_a, y[I::ii]));
        increment(row, I::p_ii_up_zero, multiply(x[I::p_ii], b[S::zero]));
        increment(row, I::p_ii_up_b, multiply(x[I::p_ii], other_b));
        increment(row, I::p_bb_up_i, multiply(subtract(a[S::p_bb], x[I::p_ii]), y[I::i]));
        increment(row, I::zero_up_i_up_i,
            add(multiply(a[S::zero], y[I::i_up_i]), multiply(x[I::zero_up_i], y[I::i])));
        increment(row, I::i_up_b_up_zero,
            add(multiply(x[I::i], subtract(b[S::b_up_zero], y[I::i_up_zero])),
                multiply(x[I::i_up_b], b[S::zero])));
        increment(row, I::i_up_zero_up_b,
            add(multiply(x[I::i], subtract(b[S::zero_up_b], y[I::zero_up_i])),
                multiply(x[I::i_up_zero], other_b)));
        increment(row, I::b_up_i_up_i,
            add(multiply(other_a, y[I::i_up_i]), multiply(x[I::b_up_i], y[I::i])));
        // Figure 13's supplementary rows: choose the singleton below/above
        // an off-path pair. Each split has exactly this one cross term.
        row.i_up_p_bb = add(row.i_up_p_bb, multiply(x[I::i], subtract(b[S::p_bb], y[I::p_ii])));
        row.p_zero_b_up_i = add(row.p_zero_b_up_i, multiply(subtract(a[S::p_zero_b], x[I::p_zero_i]), y[I::i]));
        row.zero_up_ii = add(row.zero_up_ii, multiply(a[S::zero], y[I::ii]));
        // All other rows are additive for CC; derived rows are set below.
    });
    result[S::zero_up_b] = add(result[S::zero_up_b], multiply(a[S::zero], b[S::b]));
    result[S::b_up_zero] = add(result[S::b_up_zero], multiply(a[S::b], b[S::zero]));
    result[S::zero_up_bb] = add(result[S::zero_up_bb], multiply(a[S::zero], b[S::bb]));
    derived(result);
    return result;
}

QuartetSharedCounters quartet_shared_gg(const QuartetSharedCounters& a, const QuartetSharedCounters& b) {
    auto result = merge(a, b, [&](QuartetSharedColor& row, const QuartetSharedColor& x, const QuartetSharedColor& y) {
        // Figure 12: sigma[g] means BOTH complete operand orders.
        const auto ordered = [&](const QuartetSharedCounters& first, const QuartetSharedCounters& second,
                                 const QuartetSharedColor& u, const QuartetSharedColor& v) {
            increment(row, I::zero_i, multiply(first[S::zero], v[I::i]));
            increment(row, I::i_b, multiply(u[I::i], subtract(second[S::b], v[I::i])));
            const Count zero_p_ii = multiply(first[S::zero], v[I::p_ii]);
            const Count b_p_ii = multiply(subtract(first[S::b], u[I::i]), v[I::p_ii]);
            const Count i_p_zero_b = multiply(u[I::i], subtract(second[S::p_zero_b], v[I::p_zero_i]));
            increment(row, I::zero_p_ii, zero_p_ii);
            increment(row, I::b_p_ii, b_p_ii);
            increment(row, I::i_p_zero_b, i_p_zero_b);
            increment(row, I::i_p_bb, multiply(u[I::i], subtract(second[S::p_bb], v[I::p_ii])));
            increment(row, I::i_p_bc, multiply(u[I::i], subtract(second[S::p_bc], v[I::p_i_b])));
            increment(row, I::any_zero_p_ii, zero_p_ii);
            increment(row, I::any_b_p_ii, b_p_ii);
            increment(row, I::any_i_p_zero_b, i_p_zero_b);
        };
        ordered(a, b, x, y);
        ordered(b, a, y, x);
        // sigma[1/2 n_i' n_i''] = n_i' n_i''. Combine before integer division.
        increment(row, I::ii, multiply(x[I::i], y[I::i]));
    });
    result[S::zero_b] = add(result[S::zero_b],
        add(multiply(a[S::zero], b[S::b]), multiply(b[S::zero], a[S::b])));
    clear_path(result);  // These counters are undefined for G.
    derived(result);
    return result;
}

QuartetSharedCounters quartet_shared_ig(const QuartetSharedCounters& group) {
    auto result = group;
    clear_path(result);  // Figure 12: a one-vertex path admits no upward pattern.
    derived(result);
    return result;
}

void quartet_shared_c_to_g(QuartetSharedCounters& result) {
    // Figure 12: all leaves of C lie in one child subtree of the new super root.
    result[S::zero_b] = 0;
    result[S::p_zero_b] = multiply(result[S::zero], result[S::b]);
    for (auto& row : result.colors) {
        row[I::zero_i] = row[I::ii] = row[I::i_b] = 0;
        row[I::zero_p_ii] = row[I::b_p_ii] = row[I::i_p_zero_b] = 0;
        row[I::i_p_bb] = row[I::i_p_bc] = 0;
        row[I::p_zero_i] = multiply(result[S::zero], row[I::i]);
        row[I::p_ii] = choose(row[I::i], 2);
        row[I::p_i_b] = multiply(row[I::i], subtract(result[S::b], row[I::i]));
        row[I::p_zero_p_ii] = row[I::any_zero_p_ii];
        row[I::p_b_p_ii] = row[I::any_b_p_ii];
        row[I::p_i_p_zero_b] = row[I::any_i_p_zero_b];
    }
    clear_path(result);
    derived(result);
}

}  // namespace soda2013
