// soda13.pdf, Section 7: contributions are evaluated before C -> G changes
// the boundary counters. All-zero contracted atoms contain only n_0; every
// other counter and every contribution requires a positive color.
#include "quartet_hdt.hpp"

namespace soda2013 {

ColoredQuartetHdt::ColoredQuartetHdt(const Tree& tree, Color initial_color)
    : ColoredQuartetHdt(ComponentTree::from_tree(tree), initial_color) {}

ColoredQuartetHdt::ColoredQuartetHdt(ComponentTree tree, Color initial_color)
    : ColoredHdt(std::move(tree), initial_color), counters_(hdt_.nodes().size()),
      totals_(hdt_.nodes().size()) {
    for (NodeId uid = 0; uid < hdt_.nodes().size(); ++uid) recompute(uid);
}

void ColoredQuartetHdt::recompute(NodeId uid) {
    const auto& node = hdt_.node(uid);
    QuartetCounters result;
    ResolvedQuartets total;
    if (node.left != no_node) {
        total = totals_[node.left];
        total.accumulate(totals_[node.right]);
    }
    switch (node.composition) {
    case Composition::Atom:
        if (tree_.atom(node.source).active_leaf == no_node) {
            result = quartet_leaf(0, tree_.atom(node.source).leaves);
            break;
        }
        [[fallthrough]];
    case Composition::LeafToC:
        result = quartet_leaf(colors_[node.source]);
        break;
    case Composition::IdentityI:
        break;
    case Composition::CC:
        total.agreeing.accumulate(shared_quartet_cc(counters_[node.left].shared, counters_[node.right].shared));
        total.disagreeing.accumulate(disagreement_quartet_cc(counters_[node.left], counters_[node.right]));
        result = quartet_cc(counters_[node.left], counters_[node.right]);
        break;
    case Composition::GG:
        total.agreeing.accumulate(shared_quartet_gg(counters_[node.left].shared, counters_[node.right].shared));
        total.disagreeing.accumulate(disagreement_quartet_gg(counters_[node.left], counters_[node.right]));
        result = quartet_gg(counters_[node.left], counters_[node.right]);
        break;
    case Composition::IG:
        result = quartet_ig(counters_[node.right]);
        break;
    }
    if (node.c_to_g) quartet_c_to_g(result);
    ++stats_.recomputed_components;
    stats_.color_rows += result.colors.size();
    stats_.color_pairs += result.pairs.size();
    counters_[uid] = std::move(result);
    totals_[uid] = std::move(total);
}

}  // namespace soda2013
