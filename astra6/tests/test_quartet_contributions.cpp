// Figures 13/16: compare every local contribution against four-leaf enumeration,
// and compare all nine cases after two tree orientations. The slow per-node
// coloring driver below is TEST ONLY, not a public quartet backend.
#include "quartet_topology_oracle.hpp"
#include "quartet_disagreement_contributions.hpp"

using namespace soda2013;
using namespace test_support;

struct Decorated {
    std::vector<QuartetCounters> counters;
    std::vector<QuartetCases> local, subtree;
};

Decorated decorate(const Tree& tree, const Hdt& hdt, const std::vector<Color>& colors, bool disagreement) {
    Decorated result{std::vector<QuartetCounters>(hdt.nodes().size()),
                     std::vector<QuartetCases>(hdt.nodes().size()), std::vector<QuartetCases>(hdt.nodes().size())};
    for (NodeId uid = 0; uid < hdt.nodes().size(); ++uid) {
        const auto& node = hdt.node(uid);
        auto& counters = result.counters[uid];
        if (node.composition == Composition::LeafToC) counters = quartet_leaf(colors[node.source]);
        else if (node.composition == Composition::CC) {
            result.local[uid] = disagreement ?
                disagreement_quartet_cc(result.counters[node.left], result.counters[node.right]) :
                shared_quartet_cc(result.counters[node.left].shared, result.counters[node.right].shared);
            counters = quartet_cc(result.counters[node.left], result.counters[node.right]);
        } else if (node.composition == Composition::GG) {
            result.local[uid] = disagreement ?
                disagreement_quartet_gg(result.counters[node.left], result.counters[node.right]) :
                shared_quartet_gg(result.counters[node.left].shared, result.counters[node.right].shared);
            counters = quartet_gg(result.counters[node.left], result.counters[node.right]);
        } else if (node.composition == Composition::IG) counters = quartet_ig(result.counters[node.right]);
        if (node.c_to_g) quartet_c_to_g(counters);
        result.subtree[uid] = result.local[uid];
        if (node.left != no_node) {
            result.subtree[uid].accumulate(result.subtree[node.left]);
            result.subtree[uid].accumulate(result.subtree[node.right]);
        }
    }
    (void)tree;
    return result;
}

QuartetCases local_oracle(const Tree& tree, const Hdt& hdt, const counter_oracle::ComponentGeometry& geometry,
                         const std::vector<Color>& colors, NodeId uid, bool disagreement) {
    const auto& node = hdt.node(uid);
    if (node.left == no_node) return {};
    std::vector<bool> left(tree.size());
    for (NodeId leaf : geometry[node.left].leaves) left[leaf] = true;
    const auto& leaves = geometry[uid].leaves;
    QuartetCases result;
    for (std::size_t i = 0; i < leaves.size(); ++i)
        for (std::size_t j = i + 1; j < leaves.size(); ++j)
            for (std::size_t k = j + 1; k < leaves.size(); ++k)
                for (std::size_t l = k + 1; l < leaves.size(); ++l) {
                    const std::array<NodeId, 4> selected{leaves[i], leaves[j], leaves[k], leaves[l]};
                    const int in_left = static_cast<int>(left[selected[0]]) + static_cast<int>(left[selected[1]]) +
                                        static_cast<int>(left[selected[2]]) + static_cast<int>(left[selected[3]]);
                    if (in_left == 0 || in_left == 4) continue;
                    const auto first = quartet_oracle::coloring(colors, selected);
                    if (!first.split) continue;
                    const auto second = quartet_oracle::topology(tree, selected);
                    if (!second.split || ((first.split != second.split) != disagreement) || first.type < second.type)
                        continue;
                    ++result.at(first.type, second.type);
                }
    return result;
}

QuartetCases test_coloring_pass(const Tree& first, const Tree& second, bool disagreement) {
    const Hdt hdt(second);
    QuartetCases result;
    for (NodeId v : first.preorder()) {
        if (first.is_leaf(v)) continue;
        std::vector<Color> colors(second.size());
        Color color = 0;
        for (NodeId child : first.node(v).children) {
            ++color;
            for (std::size_t i = first.leaf_begin(child); i < first.leaf_end(child); ++i)
                colors[second.leaf(first.node(first.leaves()[i]).label)] = color;
        }
        result.accumulate(decorate(second, hdt, colors, disagreement).subtree[hdt.root()]);
    }
    return result;
}

void check_contributions(bool disagreement) {
    std::mt19937_64 rng(13012013);
    std::array<bool, 9> cc_seen{}, gg_seen{}, all_cases_seen{};
    for (std::size_t n = 4; n <= 12; ++n) {
        for (int trial = 0; trial < 100; ++trial) {
            const auto tree = trial == 0 ? star(n) : trial == 1 ? caterpillar(n) : random_tree(n, rng, trial % 2 ? 2 : 12);
            const Hdt hdt(tree);
            const counter_oracle::ComponentGeometry geometry(tree, hdt);
            std::vector<Color> colors(tree.size());
            for (NodeId leaf : tree.leaves()) colors[leaf] = rng() % 5;
            const auto actual = decorate(tree, hdt, colors, disagreement);
            for (NodeId uid = 0; uid < hdt.nodes().size(); ++uid) {
                const auto expected = local_oracle(tree, hdt, geometry, colors, uid, disagreement);
                for (std::size_t c = 0; c < 9; ++c) {
                    if (actual.local[uid].values[c] != expected.values[c])
                        throw std::runtime_error(counter_oracle::diagnostic(tree, colors, uid, "case" + std::to_string(c),
                            0, actual.local[uid].values[c], expected.values[c]));
                    if (hdt.node(uid).composition == Composition::CC) cc_seen[c] = cc_seen[c] || expected.values[c];
                    if (hdt.node(uid).composition == Composition::GG) gg_seen[c] = gg_seen[c] || expected.values[c];
                }
            }
        }
    }
    for (std::size_t row = 0; row < 3; ++row)
        for (std::size_t column = 0; column <= row; ++column) {
            CHECK(cc_seen[3 * row + column]);
            CHECK(gg_seen[3 * row + column]);
        }

    // The five rooted four-leaf shapes, with every label permutation, give
    // all 26 rooted multifurcating trees (including both unresolved shapes).
    std::set<std::string> newicks;
    std::array<std::string, 4> labels{"a", "b", "c", "d"};
    do {
        const auto& a = labels[0]; const auto& b = labels[1];
        const auto& c = labels[2]; const auto& d = labels[3];
        // Child orders do not affect topology; duplicates here only broaden
        // checks of sibling ordering in the Section 4 decomposition.
        newicks.insert("((" + a + "," + b + "),(" + c + "," + d + "));" );
        newicks.insert("((" + a + "," + b + ")," + c + "," + d + ");");
        newicks.insert("(((" + a + "," + b + ")," + c + ")," + d + ");");
        newicks.insert("(" + a + "," + b + "," + c + "," + d + ");");
        newicks.insert("((" + a + "," + b + "," + c + ")," + d + ");");
    } while (std::next_permutation(labels.begin(), labels.end()));
    std::vector<Tree> trees;
    for (const auto& newick : newicks) trees.push_back(Tree::from_newick(newick));
    std::set<std::string> canonical_shapes;
    for (const auto& tree : trees) {
        std::vector<std::string> signature(tree.size());
        for (auto it = tree.preorder().rbegin(); it != tree.preorder().rend(); ++it) {
            if (tree.is_leaf(*it)) signature[*it] = tree.node(*it).label;
            else {
                std::vector<std::string> children;
                for (NodeId child : tree.node(*it).children) children.push_back(signature[child]);
                std::sort(children.begin(), children.end());
                signature[*it] = "(";
                for (const auto& child : children) signature[*it] += child + ",";
                signature[*it].back() = ')';
            }
        }
        canonical_shapes.insert(signature[tree.root()]);
    }
    CHECK(canonical_shapes.size() == 26);
    const auto check_pair = [&](const Tree& first, const Tree& second) {
        const auto forward = test_coloring_pass(first, second, disagreement);
        const auto reverse = test_coloring_pass(second, first, disagreement);
        const auto actual = complete_quartet_cases(forward, reverse);
        const auto expected = quartet_oracle::resolved_cases(first, second, disagreement);
        CHECK(actual.values == expected.values);
        const auto oracle_counts = oracle::quartet_counts(first, second);
        CHECK(actual.total() == (disagreement ? oracle_counts.B : oracle_counts.A));
        for (std::size_t c = 0; c < 9; ++c) all_cases_seen[c] = all_cases_seen[c] || actual.values[c];
    };
    for (const auto& first : trees) for (const auto& second : trees) check_pair(first, second);
    for (std::size_t n = 4; n <= 12; ++n)
        for (int trial = 0; trial < 20; ++trial) check_pair(random_tree(n, rng), random_tree(n, rng));
    for (bool seen : all_cases_seen) CHECK(seen);
    std::cout << "Figure " << (disagreement ? 16 : 13) << ": OK (900 local-coloring tests; " << trees.size() * trees.size()
              << " four-leaf ordered pairs; 180 random pairs; all nine cases)\n";
}

int main() try {
    // Figure 16 gamma-beta, ((ij),i,0), CC cut 5: the scalar already
    // sums over colors. One zero below and ((1,2),1) above give ONE quartet.
    auto one = quartet_leaf(1), two = quartet_leaf(2);
    quartet_c_to_g(one);
    quartet_c_to_g(two);
    auto paired = quartet_ig(quartet_gg(one, two));
    quartet_c_to_g(paired);
    const auto upper = quartet_ig(quartet_gg(paired, one));
    const auto regression = disagreement_quartet_cc(quartet_leaf(0), upper);
    CHECK(regression.at(QuartetType::gamma, QuartetType::beta) == 1);
    CHECK(regression.total() == 1);
    check_contributions(false);
    check_contributions(true);
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
