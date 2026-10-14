// soda13.pdf, Figures 14–15. All five operations and every distinct counter.
// Ambiguous/incorrect printed cells are resolved from their Figure 14 topology
// constraints, with an independent enumeration oracle; see Trial 8.
#include "quartet_disagreement_counters.hpp"

#include <utility>

namespace soda2013 {
namespace {
using S = SharedScalar;
using I = SharedIndexed;
using DS = DisagreementScalar;
using DI = DisagreementIndexed;

std::size_t pair_offset(std::size_t size, std::size_t i, std::size_t j) {
    if (i >= size || j >= size || i == j) throw std::out_of_range("distinct present-color indices required");
    if (i > j) std::swap(i, j);
    return i * (2 * size - i - 1) / 2 + (j - i - 1);
}

void increment(DisagreementColor& row, DI counter, Count value) { row[counter] = add(row[counter], value); }

void derived(QuartetCounters& result) {
    result[DS::b_p_b_c] = result[DS::any_b_p_b_c] = result[DS::b_up_b_c] = 0;
    result[DS::p_b_c_up_b] = result[DS::b_up_b_up_c] = result.b_up_p_b_c = 0;
    for (const auto& row : result.colors) {
        result[DS::b_p_b_c] = add(result[DS::b_p_b_c], row[DI::i_p_i_b]);
        result[DS::any_b_p_b_c] = add(result[DS::any_b_p_b_c], row[DI::any_i_p_i_b]);
        result[DS::b_up_b_c] = add(result[DS::b_up_b_c], row[DI::i_up_i_b]);
        result[DS::p_b_c_up_b] = add(result[DS::p_b_c_up_b], row[DI::p_i_b_up_i]);
        result[DS::b_up_b_up_c] = add(result[DS::b_up_b_up_c], row[DI::i_up_i_up_b]);
        result.b_up_p_b_c = add(result.b_up_p_b_c, row[DI::i_up_p_i_b]);
    }
}

void clear_path(QuartetCounters& result) {
    result.b_up_c_up_b = 0;
    for (auto& row : result.colors) {
        for (std::size_t i = static_cast<std::size_t>(DI::zero_up_i_b); i < disagreement_indexed_count; ++i)
            row.values[i] = 0;
        row.b_up_p_i_b = row.b_up_i_b = row.b_up_i_c = row.i_up_zero_i = 0;
    }
}

struct Aligned {
    std::vector<NodeId> a, b;
};

Aligned initialize(QuartetCounters& result, const QuartetCounters& a, const QuartetCounters& b) {
    const auto size = result.shared.colors.size();
    result.colors.resize(size);
    result.pairs.resize(static_cast<std::size_t>(choose(size, 2)));
    for (std::size_t c = 0; c < disagreement_scalar_count; ++c) result.values[c] = add(a.values[c], b.values[c]);
    result.b_up_c_up_b = add(a.b_up_c_up_b, b.b_up_c_up_b);
    Aligned aligned{std::vector<NodeId>(size, no_node), std::vector<NodeId>(size, no_node)};
    std::size_t i = 0, j = 0;
    for (std::size_t c = 0; c < size; ++c) {
        const Color color = result.shared.colors[c].color;
        if (i < a.shared.colors.size() && a.shared.colors[i].color == color) aligned.a[c] = i++;
        if (j < b.shared.colors.size() && b.shared.colors[j].color == color) aligned.b[c] = j++;
        auto& row = result.colors[c];
        for (const auto& operand : {std::pair<const QuartetCounters*, NodeId>{&a, aligned.a[c]}, {&b, aligned.b[c]}}) {
            if (operand.second == no_node) continue;
            const auto& source = operand.first->colors[operand.second];
            for (std::size_t k = 0; k < disagreement_indexed_count; ++k) row.values[k] = add(row.values[k], source.values[k]);
            row.b_up_p_i_b = add(row.b_up_p_i_b, source.b_up_p_i_b);
            row.b_up_i_b = add(row.b_up_i_b, source.b_up_i_b);
            row.b_up_i_c = add(row.b_up_i_c, source.b_up_i_c);
            row.i_up_zero_i = add(row.i_up_zero_i, source.i_up_zero_i);
        }
    }
    return aligned;
}

const QuartetSharedColor zero_color;
const DisagreementPair zero_pair;
const QuartetSharedColor& color_at(const QuartetCounters& counters, NodeId index) {
    return index == no_node ? zero_color : counters.shared.colors[index];
}
const DisagreementPair& pair_at(const QuartetCounters& counters, NodeId first, NodeId second) {
    return first == no_node || second == no_node ? zero_pair : counters.pair(first, second);
}

void compose_pairs(QuartetCounters& result, const QuartetCounters& a, const QuartetCounters& b,
                   const Aligned& aligned, bool gg) {
    for (std::size_t i = 0; i < result.colors.size(); ++i)
        for (std::size_t j = i + 1; j < result.colors.size(); ++j) {
            const auto& x = pair_at(a, aligned.a[i], aligned.a[j]);
            const auto& y = pair_at(b, aligned.b[i], aligned.b[j]);
            auto& row = result.pair(i, j);
            row = {add(x.at, y.at), add(x.inside, y.inside)};
            if (gg) {
                const auto& ai = color_at(a, aligned.a[i]); const auto& aj = color_at(a, aligned.a[j]);
                const auto& bi = color_at(b, aligned.b[i]); const auto& bj = color_at(b, aligned.b[j]);
                row.at = add(row.at, add(multiply(ai[I::i], bj[I::i]), multiply(bi[I::i], aj[I::i])));
            }
        }
}
}  // namespace

const DisagreementPair& QuartetCounters::pair(std::size_t i, std::size_t j) const {
    return pairs.at(pair_offset(colors.size(), i, j));
}
DisagreementPair& QuartetCounters::pair(std::size_t i, std::size_t j) {
    return pairs.at(pair_offset(colors.size(), i, j));
}

QuartetCounters quartet_leaf(Color color, Count zero_multiplicity) {
    QuartetCounters result;
    result.shared = quartet_shared_leaf(color, zero_multiplicity);
    result.colors.resize(result.shared.colors.size());
    return result;
}

QuartetCounters quartet_cc(const QuartetCounters& a, const QuartetCounters& b) {
    QuartetCounters result;
    result.shared = quartet_shared_cc(a.shared, b.shared);
    const auto aligned = initialize(result, a, b);
    compose_pairs(result, a, b, aligned, false);
    const Count a0 = a.shared[S::zero], b0 = b.shared[S::zero];
    for (std::size_t i = 0; i < result.colors.size(); ++i) {
        const auto& x = color_at(a, aligned.a[i]); const auto& y = color_at(b, aligned.b[i]);
        const Count xi = x[I::i], yi = y[I::i];
        const Count other_a = subtract(a.shared[S::b], xi), other_b = subtract(b.shared[S::b], yi);
        auto& row = result.colors[i];
        // Figure 15 bracketed rows. Each resolved triple has four crossing
        // cuts: inner pair below, first inner leaf below, second inner leaf
        // below, or outer leaf below with an off-path inner pair above.
        increment(row, DI::any_i_p_i_b, add(add(multiply(multiply(xi, other_a), yi), multiply(xi, y[I::b_up_i])),
            add(multiply(other_a, y[I::i_up_i]), multiply(xi, y[I::p_i_b]))));
        increment(row, DI::any_zero_p_i_b, add(add(multiply(multiply(xi, other_a), b0),
            multiply(xi, subtract(b.shared[S::b_up_zero], y[I::i_up_zero]))),
            add(multiply(other_a, y[I::i_up_zero]), multiply(a0, y[I::p_i_b]))));
        increment(row, DI::any_i_p_zero_i, add(add(multiply(multiply(a0, xi), yi), multiply(a0, y[I::i_up_i])),
            add(multiply(xi, y[I::zero_up_i]), multiply(xi, y[I::p_zero_i]))));
        increment(row, DI::any_b_p_zero_i, add(add(multiply(multiply(a0, xi), other_b), multiply(a0, y[I::i_up_b])),
            add(multiply(xi, subtract(b.shared[S::zero_up_b], y[I::zero_up_i])), multiply(other_a, y[I::p_zero_i]))));

        // Figure 15, all C-only path rows.
        increment(row, DI::zero_up_i_b, multiply(a0, y[I::i_b]));
        increment(row, DI::b_up_zero_i, multiply(other_a, y[I::zero_i]));
        increment(row, DI::i_up_i_b, multiply(xi, y[I::i_b]));
        increment(row, DI::p_i_b_up_zero, multiply(x[I::p_i_b], b0));
        increment(row, DI::p_i_b_up_i, multiply(x[I::p_i_b], yi));
        increment(row, DI::p_zero_i_up_i, multiply(x[I::p_zero_i], yi));
        increment(row, DI::p_zero_i_up_b, multiply(x[I::p_zero_i], other_b));
        increment(row, DI::zero_up_p_i_b, multiply(a0, y[I::p_i_b]));
        increment(row, DI::i_up_p_i_b, multiply(xi, y[I::p_i_b]));
        increment(row, DI::b_up_p_zero_i, multiply(other_a, y[I::p_zero_i]));
        increment(row, DI::i_up_p_zero_i, multiply(xi, y[I::p_zero_i]));
        increment(row, DI::b_up_i_up_zero, add(multiply(x[I::b_up_i], b0), multiply(other_a, y[I::i_up_zero])));
        increment(row, DI::i_up_i_up_zero, add(multiply(x[I::i_up_i], b0), multiply(xi, y[I::i_up_zero])));
        increment(row, DI::zero_up_b_up_i, add(multiply(subtract(a.shared[S::zero_up_b], x[I::zero_up_i]), yi),
            multiply(a0, y[I::b_up_i])));
        increment(row, DI::b_up_zero_up_i, add(multiply(subtract(a.shared[S::b_up_zero], x[I::i_up_zero]), yi),
            multiply(other_a, y[I::zero_up_i])));
        increment(row, DI::i_up_zero_up_i, add(multiply(x[I::i_up_zero], yi), multiply(xi, y[I::zero_up_i])));
        increment(row, DI::zero_up_i_up_b, add(multiply(a0, y[I::i_up_b]), multiply(x[I::zero_up_i], other_b)));
        increment(row, DI::i_up_i_up_b, add(multiply(x[I::i_up_i], other_b), multiply(xi, y[I::i_up_b])));
        row.i_up_zero_i = add(row.i_up_zero_i, multiply(xi, y[I::zero_i]));

        // Figure 15's starred b↑(ic), plus Figure 16's analogous missing rows.
        for (std::size_t j = 0; j < result.colors.size(); ++j) {
            if (i == j) continue;
            const Count aj = color_at(a, aligned.a[j])[I::i];
            const auto& pair = pair_at(b, aligned.b[i], aligned.b[j]);
            increment(row, DI::b_up_p_i_c, multiply(subtract(other_a, aj), pair.inside));
            row.b_up_p_i_b = add(row.b_up_p_i_b, multiply(aj, pair.inside));
            row.b_up_i_c = add(row.b_up_i_c, multiply(subtract(other_a, aj), pair.at));
            row.b_up_i_b = add(row.b_up_i_b, multiply(aj, pair.at));
        }
        // Figure 16's b↑c↑b: the equal-color endpoints lie on opposite sides
        // of the cut, with the unlike middle leaf on either side.
        result.b_up_c_up_b = add(result.b_up_c_up_b,
            add(multiply(x[I::i_up_b], yi), multiply(xi, y[I::b_up_i])));
    }
    derived(result);
    return result;
}

QuartetCounters quartet_gg(const QuartetCounters& a, const QuartetCounters& b) {
    QuartetCounters result;
    result.shared = quartet_shared_gg(a.shared, b.shared);
    const auto aligned = initialize(result, a, b);
    compose_pairs(result, a, b, aligned, true);
    for (std::size_t i = 0; i < result.colors.size(); ++i) {
        auto& row = result.colors[i];
        const auto& x = color_at(a, aligned.a[i]); const auto& y = color_at(b, aligned.b[i]);
        const auto ordered = [&](const QuartetCounters& first, const QuartetCounters& second,
                                 const QuartetSharedColor& u, const QuartetSharedColor& v) {
            const Count other = subtract(first.shared[S::b], u[I::i]);
            const Count i0i = multiply(u[I::i], v[I::p_zero_i]);
            const Count b0i = multiply(other, v[I::p_zero_i]);
            const Count zero_ib = multiply(first.shared[S::zero], v[I::p_i_b]);
            const Count i_ib = multiply(u[I::i], v[I::p_i_b]);
            increment(row, DI::i_p_zero_i, i0i);
            increment(row, DI::b_p_zero_i, b0i);
            increment(row, DI::zero_p_i_b, zero_ib);
            increment(row, DI::i_p_i_b, i_ib);
            increment(row, DI::any_i_p_i_b, i_ib);
            increment(row, DI::any_zero_p_i_b, zero_ib);
            increment(row, DI::any_i_p_zero_i, i0i);
            increment(row, DI::any_b_p_zero_i, b0i);
            (void)second;
        };
        ordered(a, b, x, y);
        ordered(b, a, y, x);
        for (std::size_t j = 0; j < result.colors.size(); ++j) {
            if (j == i) continue;
            const Count aj = color_at(a, aligned.a[j])[I::i], bj = color_at(b, aligned.b[j])[I::i];
            const Count ap = pair_at(a, aligned.a[i], aligned.a[j]).inside;
            const Count bp = pair_at(b, aligned.b[i], aligned.b[j]).inside;
            increment(row, DI::b_p_i_c,
                add(multiply(subtract(subtract(a.shared[S::b], x[I::i]), aj), bp),
                    multiply(subtract(subtract(b.shared[S::b], y[I::i]), bj), ap)));
            // In b(ib) BOTH b leaves have the same color. The printed
            // (n_b-n_i)n_(ib) also includes b(ic); separate that contribution.
            increment(row, DI::b_p_i_b, add(multiply(aj, bp), multiply(bj, ap)));
        }
    }
    clear_path(result);
    derived(result);
    return result;
}

QuartetCounters quartet_ig(const QuartetCounters& group) {
    auto result = group;
    result.shared = quartet_shared_ig(group.shared);
    clear_path(result);
    derived(result);
    return result;
}

void quartet_c_to_g(QuartetCounters& result) {
    for (std::size_t i = 0; i < result.colors.size(); ++i) {
        auto& row = result.colors[i];
        for (std::size_t c = 0; c <= static_cast<std::size_t>(DI::b_p_i_b); ++c) row.values[c] = 0;
        row[DI::p_i_p_i_b] = row[DI::any_i_p_i_b];
        row[DI::p_i_p_zero_i] = row[DI::any_i_p_zero_i];
        row[DI::p_zero_p_i_b] = row[DI::any_zero_p_i_b];
        row[DI::p_b_p_zero_i] = row[DI::any_b_p_zero_i];
        for (std::size_t j = i + 1; j < result.colors.size(); ++j) {
            // Figure 14 fixes BOTH i and j; the printed n_i(n_b-n_i) would
            // instead be the aggregate (ib) counter of Figure 11.
            result.pair(i, j) = {0, multiply(result.shared.colors[i][I::i], result.shared.colors[j][I::i])};
        }
    }
    result[DS::p_b_p_b_c] = result[DS::any_b_p_b_c];
    quartet_shared_c_to_g(result.shared);
    clear_path(result);
    derived(result);
}

}  // namespace soda2013
