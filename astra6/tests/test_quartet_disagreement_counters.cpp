// Figures 14–15: semantic checks of all scalar, per-color and fixed-color-pair
// counters. Repeated bullets mean the SAME color, including across parentheses.
// Figure 16's missing path dependencies are tested by the same definitions.
#include "quartet_counter_oracle.hpp"
#include "quartet_disagreement_counters.hpp"

using namespace soda2013;
using namespace test_support;
using namespace counter_oracle;

int main() try {
    {
        const auto tree = Tree::from_newick("((a,b),c);");
        const Hdt hdt(tree);
        const ComponentGeometry geometry(tree, hdt);
        std::vector<Color> colors(tree.size());
        colors[tree.leaf("a")] = colors[tree.leaf("c")] = 1;
        colors[tree.leaf("b")] = 2;
        CHECK(Pattern("[b(bc)]").count(geometry[hdt.root()], colors) == 1);
        CHECK(Pattern("[bc]").count(geometry[hdt.root()], colors) == 2);
    }
    std::vector<Pattern> scalar_patterns, indexed_patterns;
    for (auto subscript : disagreement_scalar_subscripts) scalar_patterns.emplace_back(subscript);
    for (auto subscript : disagreement_indexed_subscripts) indexed_patterns.emplace_back(subscript);
    const Pattern ij("ij"), p_ij("(ij)"), b_up_p_bc("b^(bc)"), b_up_c_up_b("b^c^b");
    const std::array<Pattern, 4> indexed_extra{Pattern("b^(ib)"), Pattern("b^ib"), Pattern("b^ic"), Pattern("i^0i")};
    std::array<bool, disagreement_scalar_count> scalar_seen{};
    std::array<bool, disagreement_indexed_count> indexed_seen{};
    std::array<bool, 8> extra_seen{};  // Two pair types, two scalar/four indexed dependencies.
    std::array<bool, 5> actions{};
    std::size_t checked = 0;
    const auto check = [&](const Tree& tree, const std::vector<Color>& colors, NodeId uid,
                           const Geometry& geometry, const QuartetCounters& actual) {
        for (std::size_t c = 0; c < disagreement_scalar_count; ++c) {
            const Count expected = scalar_patterns[c].count(geometry, colors);
            if (actual.values[c] != expected)
                throw std::runtime_error(diagnostic(tree, colors, uid, disagreement_scalar_subscripts[c], 0, actual.values[c], expected));
            scalar_seen[c] = scalar_seen[c] || expected;
            ++checked;
        }
        const Count first_extra = b_up_p_bc.count(geometry, colors), second_extra = b_up_c_up_b.count(geometry, colors);
        CHECK(first_extra == actual.b_up_p_b_c && second_extra == actual.b_up_c_up_b);
        extra_seen[2] = extra_seen[2] || first_extra;
        extra_seen[3] = extra_seen[3] || second_extra;
        CHECK(actual.colors.size() == actual.shared.colors.size());
        CHECK(actual.pairs.size() == choose(actual.colors.size(), 2));
        for (std::size_t i = 0; i < actual.colors.size(); ++i) {
            const Color color = actual.shared.colors[i].color;
            const auto& row = actual.colors[i];
            for (std::size_t c = 0; c < disagreement_indexed_count; ++c) {
                const Count expected = indexed_patterns[c].count(geometry, colors, color);
                if (row.values[c] != expected)
                    throw std::runtime_error(diagnostic(tree, colors, uid, disagreement_indexed_subscripts[c], color, row.values[c], expected));
                indexed_seen[c] = indexed_seen[c] || expected;
                ++checked;
            }
            const std::array<Count, 4> extras{row.b_up_p_i_b, row.b_up_i_b, row.b_up_i_c, row.i_up_zero_i};
            for (std::size_t c = 0; c < extras.size(); ++c) {
                const Count expected = indexed_extra[c].count(geometry, colors, color);
                if (expected != extras[c])
                    throw std::runtime_error(diagnostic(tree, colors, uid, indexed_extra[c].expression(), color, extras[c], expected));
                extra_seen[4 + c] = extra_seen[4 + c] || expected;
            }
            for (std::size_t j = i + 1; j < actual.colors.size(); ++j) {
                const Color second_color = actual.shared.colors[j].color;
                const Count at = ij.count(geometry, colors, color, second_color);
                const Count inside = p_ij.count(geometry, colors, color, second_color);
                if (actual.pair(i, j).at != at || actual.pair(i, j).inside != inside)
                    throw std::runtime_error(diagnostic(tree, colors, uid, "(ij) j=" + std::to_string(second_color),
                        color, actual.pair(i, j).inside, inside));
                CHECK(actual.pair(i, j).at == actual.pair(j, i).at);
                CHECK(actual.pair(i, j).inside == actual.pair(j, i).inside);
                extra_seen[0] = extra_seen[0] || at;
                extra_seen[1] = extra_seen[1] || inside;
                checked += 2;
            }
        }
    };
    std::mt19937_64 rng(14012013);
    for (std::size_t n = 1; n <= 12; ++n)
        for (int trial = 0; trial < 50; ++trial) {
            const auto tree = trial == 0 ? star(n) : trial == 1 ? caterpillar(n) : random_tree(n, rng, trial % 2 ? 2 : 12);
            const Hdt hdt(tree);
            const ComponentGeometry geometry(tree, hdt);
            std::vector<Color> colors(tree.size());
            for (NodeId leaf : tree.leaves()) colors[leaf] = (rng() % 4) * (trial % 3 ? 1 : 1000000);
            std::vector<QuartetCounters> counters(hdt.nodes().size());
            for (NodeId uid = 0; uid < hdt.nodes().size(); ++uid) {
                const auto& node = hdt.node(uid);
                switch (node.composition) {
                case Composition::LeafToC: counters[uid] = quartet_leaf(colors[node.source]); actions[0] = true; break;
                case Composition::IdentityI: break;
                case Composition::CC: counters[uid] = quartet_cc(counters[node.left], counters[node.right]); actions[1] = true; break;
                case Composition::GG: counters[uid] = quartet_gg(counters[node.left], counters[node.right]); actions[2] = true; break;
                case Composition::IG: counters[uid] = quartet_ig(counters[node.right]); actions[3] = true; break;
                case Composition::Atom: throw std::logic_error("unexpected original-tree atom");
                }
                if (node.c_to_g) {
                    auto before = geometry[uid];
                    before.type = ComponentType::C;
                    for (NodeId leaf : before.leaves) {
                        NodeId v = leaf;
                        while (!before.path[v]) { v = tree.parent(v); CHECK(v != no_node); }
                        before.branch[leaf] = v;
                    }
                    check(tree, colors, uid, before, counters[uid]);
                    quartet_c_to_g(counters[uid]);
                    actions[4] = true;
                }
                check(tree, colors, uid, geometry[uid], counters[uid]);
            }
        }
    for (std::size_t c = 0; c < scalar_seen.size(); ++c)
        if (!scalar_seen[c]) throw std::runtime_error("unexercised scalar " + std::string(disagreement_scalar_subscripts[c]));
    for (std::size_t c = 0; c < indexed_seen.size(); ++c)
        if (!indexed_seen[c]) throw std::runtime_error("unexercised indexed " + std::string(disagreement_indexed_subscripts[c]));
    for (bool seen : extra_seen) CHECK(seen);
    for (bool action : actions) CHECK(action);
    std::cout << "disagreement counters: OK (41 types + 6 Figure 16 dependencies, " << checked << " base comparisons)\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
