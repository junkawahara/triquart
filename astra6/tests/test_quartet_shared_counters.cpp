// Each of Figure 11's 53 types is compared to the subscript's direct meaning
// on EVERY HDT component, covering all five Figure 12 actions and arbitrary
// colors. No production enumeration or assumed agreement with Figure 12.
#include "quartet_counter_oracle.hpp"

using namespace soda2013;
using namespace test_support;
using namespace counter_oracle;

int main() try {
    std::vector<Pattern> scalar_patterns, indexed_patterns;
    for (auto subscript : shared_scalar_subscripts) scalar_patterns.emplace_back(subscript);
    for (auto subscript : shared_indexed_subscripts) indexed_patterns.emplace_back(subscript);
    std::array<bool, shared_scalar_count> scalar_seen{};
    std::array<bool, shared_indexed_count> indexed_seen{};
    const Pattern c_up_bb("c^bb"), i_up_p_bb("i^(bb)"), p_zero_b_up_i("(0b)^i"), zero_up_ii("0^ii");
    std::array<bool, 4> supplementary_seen{};
    std::mt19937_64 rng(11012013);
    std::size_t checked = 0;
    std::array<bool, 5> actions{};
    const auto check = [&](const Tree& tree, const std::vector<Color>& colors, NodeId uid,
                           const Geometry& geometry, const QuartetSharedCounters& actual) {
        for (std::size_t i = 0; i < shared_scalar_count; ++i) {
            const Count expected = scalar_patterns[i].count(geometry, colors);
            if (actual.values[i] != expected)
                throw std::runtime_error(diagnostic(tree, colors, uid, shared_scalar_subscripts[i], 0, actual.values[i], expected));
            scalar_seen[i] = scalar_seen[i] || expected;
            ++checked;
        }
        const Count supplementary_scalar = c_up_bb.count(geometry, colors);
        CHECK(actual.c_up_bb == supplementary_scalar);
        supplementary_seen[0] = supplementary_seen[0] || supplementary_scalar;
        std::set<Color> present;
        for (NodeId leaf : geometry.leaves) if (colors[leaf]) present.insert(colors[leaf]);
        CHECK(actual.colors.size() == present.size());
        auto color = present.begin();
        for (const auto& row : actual.colors) {
            CHECK(row.color == *color++);
            for (std::size_t i = 0; i < shared_indexed_count; ++i) {
                const Count expected = indexed_patterns[i].count(geometry, colors, row.color);
                if (row.values[i] != expected)
                    throw std::runtime_error(diagnostic(tree, colors, uid, shared_indexed_subscripts[i], row.color, row.values[i], expected));
                indexed_seen[i] = indexed_seen[i] || expected;
                ++checked;
            }
            const std::array<Count, 3> expected_extra{
                i_up_p_bb.count(geometry, colors, row.color), p_zero_b_up_i.count(geometry, colors, row.color),
                zero_up_ii.count(geometry, colors, row.color)};
            CHECK(row.i_up_p_bb == expected_extra[0]);
            CHECK(row.p_zero_b_up_i == expected_extra[1]);
            CHECK(row.zero_up_ii == expected_extra[2]);
            for (std::size_t i = 0; i < expected_extra.size(); ++i)
                supplementary_seen[i + 1] = supplementary_seen[i + 1] || expected_extra[i];
        }
    };
    for (std::size_t n = 1; n <= 12; ++n) {
        for (int trial = 0; trial < 50; ++trial) {
            const auto tree = trial == 0 ? star(n) : trial == 1 ? caterpillar(n) : random_tree(n, rng, trial % 2 ? 2 : 12);
            const Hdt hdt(tree);
            const ComponentGeometry geometry(tree, hdt);
            std::vector<Color> colors(tree.size());
            for (NodeId leaf : tree.leaves()) colors[leaf] = (rng() % 4) * (trial % 3 ? 1 : 1000000);
            std::vector<QuartetSharedCounters> counters(hdt.nodes().size());
            for (NodeId uid = 0; uid < hdt.nodes().size(); ++uid) {
                const auto& node = hdt.node(uid);
                switch (node.composition) {
                case Composition::LeafToC: counters[uid] = quartet_shared_leaf(colors[node.source]); actions[0] = true; break;
                case Composition::IdentityI: break;
                case Composition::CC: counters[uid] = quartet_shared_cc(counters[node.left], counters[node.right]); actions[1] = true; break;
                case Composition::GG: counters[uid] = quartet_shared_gg(counters[node.left], counters[node.right]); actions[2] = true; break;
                case Composition::IG: counters[uid] = quartet_shared_ig(counters[node.right]); actions[3] = true; break;
                case Composition::Atom: throw std::logic_error("unexpected original-tree atom");
                }
                if (node.c_to_g) {
                    // Check the native C table before C -> G can reset rows.
                    auto before = geometry[uid];
                    before.type = ComponentType::C;
                    for (NodeId leaf : before.leaves) {
                        NodeId v = leaf;
                        while (!before.path[v]) { v = tree.parent(v); CHECK(v != no_node); }
                        before.branch[leaf] = v;
                    }
                    check(tree, colors, uid, before, counters[uid]);
                    quartet_shared_c_to_g(counters[uid]);
                    actions[4] = true;
                }
                check(tree, colors, uid, geometry[uid], counters[uid]);
            }
        }
    }
    for (bool action : actions) CHECK(action);
    for (bool seen : supplementary_seen) CHECK(seen);
    for (std::size_t i = 0; i < scalar_seen.size(); ++i)
        if (!scalar_seen[i]) throw std::runtime_error("unexercised scalar counter " + std::string(shared_scalar_subscripts[i]));
    for (std::size_t i = 0; i < indexed_seen.size(); ++i)
        if (!indexed_seen[i]) throw std::runtime_error("unexercised indexed counter " + std::string(shared_indexed_subscripts[i]));
    const auto zeros = quartet_shared_leaf(0, 1000000);
    CHECK(zeros[SharedScalar::zero] == 1000000 && zeros.colors.empty());
    expect_exception<std::invalid_argument>([] { (void)quartet_shared_leaf(1, 2); });
    std::cout << "quartet shared counters: OK (53 types + 4 Figure 13 dependencies, " << checked << " base comparisons)\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
