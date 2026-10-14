// soda13.pdf, Section 7.1–7.2, Figures 10–12. In CC, a=C1 is BELOW b=C2.
// In GG, both argument orders of sigma[f] are evaluated using the original
// operands. C -> G changes the pair counters after computing contributions.
#include "triplet_counters.hpp"

namespace soda2013 {
namespace {

TripletCounters compose(const TripletCounters& a, const TripletCounters& b, bool cc) {
    TripletCounters result;
    result.zero = add(a.zero, b.zero);
    result.positive = add(a.positive, b.positive);
    result.shared = {add(a.shared.resolved, b.shared.resolved),
                     add(a.shared.unresolved, b.shared.unresolved)};
    result.colors.reserve(a.colors.size() + b.colors.size());
    const TripletColorCounters zero;
    std::size_t i = 0, j = 0;
    Count twice_different_at = 0;
    while (i < a.colors.size() || j < b.colors.size()) {
        const Color color = j == b.colors.size() ? a.colors[i].color :
                            i == a.colors.size() ? b.colors[j].color :
                            std::min(a.colors[i].color, b.colors[j].color);
        const auto& x = i < a.colors.size() && a.colors[i].color == color ? a.colors[i++] : zero;
        const auto& y = j < b.colors.size() && b.colors[j].color == color ? b.colors[j++] : zero;
        const Count other_a = subtract(a.positive, x.leaves);
        const Count other_b = subtract(b.positive, y.leaves);
        TripletColorCounters row;
        row.color = color;
        row.leaves = add(x.leaves, y.leaves);
        row.same_inside = add(x.same_inside, y.same_inside);
        row.different_at = add(x.different_at, y.different_at);

        if (cc) {
            // Figure 12: n_i↑• and additive n_(ii), n_i•.
            row.i_up_other = add(add(x.i_up_other, y.i_up_other), multiply(x.leaves, other_b));
            // Figure 10, resolved CC splits 1,2; 3; 4 respectively.
            const Count splits12 = multiply(x.leaves, y.i_up_other);
            const Count split3 = multiply(choose(x.leaves, 2), other_b);
            const Count split4 = multiply(other_a, y.same_inside);
            result.shared.resolved = add(result.shared.resolved, add(add(splits12, split3), split4));
            // Figure 10, unresolved CC: the two upper leaves branch at the
            // SAME external-path node into two off-path subtrees.
            result.shared.unresolved = add(result.shared.unresolved,
                multiply(x.leaves, subtract(b.different_at, y.different_at)));
        } else {
            // Figure 12, n_i•: g(G1,G2) + g(G2,G1).
            row.different_at = add(row.different_at,
                add(multiply(x.leaves, other_b), multiply(y.leaves, other_a)));
            // Figure 10, resolved GG, both orders of sigma.
            result.shared.resolved = add(result.shared.resolved,
                add(multiply(x.same_inside, other_b), multiply(y.same_inside, other_a)));
            // Figure 10, unresolved GG, both orders of sigma.
            result.shared.unresolved = add(result.shared.unresolved,
                add(multiply(x.leaves, subtract(b.different_at, y.different_at)),
                    multiply(y.leaves, subtract(a.different_at, x.different_at))));
        }
        twice_different_at = add(twice_different_at, row.different_at);
        result.colors.push_back(row);
    }
    if (twice_different_at % 2) throw std::logic_error("nonintegral unordered pair counter");
    result.different_at = twice_different_at / 2;
    return result;
}

void transform_c_to_g(TripletCounters& counters) {
    // Figure 12: all leaves of C are now in ONE child of the super root of G.
    // No pair is rooted at that super root; every equal-color pair is inside.
    counters.different_at = 0;
    for (auto& row : counters.colors) {
        row.same_inside = choose(row.leaves, 2);
        row.different_at = 0;
        row.i_up_other = 0;  // Undefined for G; never consumed on a G input.
    }
}

}  // namespace

ColoredTripletHdt::ColoredTripletHdt(const Tree& tree, Color initial_color)
    : ColoredTripletHdt(ComponentTree::from_tree(tree), initial_color) {}

ColoredTripletHdt::ColoredTripletHdt(ComponentTree tree, Color initial_color)
    : ColoredHdt(std::move(tree), initial_color), counters_(hdt_.nodes().size()) {
    for (NodeId uid = 0; uid < hdt_.nodes().size(); ++uid) recompute(uid);
}

void ColoredTripletHdt::recompute(NodeId uid) {
    const auto& node = hdt_.node(uid);
    TripletCounters result;
    switch (node.composition) {
    case Composition::Atom:
        if (tree_.atom(node.source).active_leaf == no_node) {
            result.zero = tree_.atom(node.source).leaves;
            break;
        }
        [[fallthrough]];  // An incontractible G atom can also be a colored leaf.
    case Composition::LeafToC:
        if (colors_[node.source] == 0) result.zero = 1;
        else {
            result.positive = 1;
            result.colors.push_back({colors_[node.source], 1, 0, 0, 0});
        }
        break;
    case Composition::IdentityI:
        break;  // An internal vertex contains no leaves.
    case Composition::CC:
        result = compose(counters_[node.left], counters_[node.right], true);
        break;
    case Composition::GG:
        result = compose(counters_[node.left], counters_[node.right], false);
        break;
    case Composition::IG:
        // Figure 12: all relevant counters are copied from G; the newly
        // introduced external path is just the I vertex, so no upward pair.
        result = counters_[node.right];
        break;
    }
    if (node.c_to_g) transform_c_to_g(result);
    ++stats_.recomputed_components;
    stats_.color_rows += result.colors.size();
    counters_[uid] = std::move(result);
}

}  // namespace soda2013
