// tests/test_tree.cpp — Newick parsing / printing, unary contraction, deep
// trees (the parser, printer and contraction must not recurse).
#include <cstdio>
#include <string>

#include "../distance.hpp"
#include "../tree.hpp"

using namespace soda13;

static int fails = 0;
#define CHECK(c)                                                              \
    do {                                                                      \
        if (!(c)) {                                                           \
            ++fails;                                                          \
            std::printf("CHECK failed at line %d: %s\n", __LINE__, #c);       \
        }                                                                     \
    } while (0)

int main() {
    // round trip
    Tree t = Tree::from_newick("((a,b),(c,d,e));");
    CHECK(t.to_newick() == "((a,b),(c,d,e));");
    CHECK(t.num_leaves() == 5);
    CHECK(t.degree_max() == 3);
    // branch lengths, internal labels, quotes, whitespace
    Tree u = Tree::from_newick(" ( ( 'a a':0.1 , b:2 )x:3 , (c,d,e)y ) root ; ");
    CHECK(u.to_newick() == "((a a,b),(c,d,e));");
    // unary contraction
    Tree w = Tree::from_newick("(((a)),((b,c)));");
    w.contract_unary();
    CHECK(w.to_newick() == "(a,(b,c));");
    // errors
    bool threw = false;
    try {
        Tree::from_newick("((a,b),(c,d);");
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    CHECK(threw);
    threw = false;
    try {
        Tree a = Tree::from_newick("(a,b,c);"), b = Tree::from_newick("(a,b,d);");
        assign_leaf_ids(a, b);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    CHECK(threw);
    // deep caterpillar: 200k nested parentheses
    {
        const int n = 200000;
        std::string s;
        for (int i = 0; i < n - 1; ++i) s += '(';
        s += "L0";
        for (int i = 1; i < n; ++i) s += ",L" + std::to_string(i) + ")";
        s += ';';
        Tree deep = Tree::from_newick(s);
        CHECK(deep.num_leaves() == n);
        std::string back = deep.to_newick();
        CHECK(back == s);
        deep.contract_unary();
        CHECK(deep.num_leaves() == n);
        // distance of a tree with itself: zero
        Tree same = Tree::from_newick(s);
        CHECK(triplet_distance(deep, same) == 0);
    }
    // small hand-checked example: two trees differing in one quartet
    {
        Tree a = Tree::from_newick("((a,b),(c,d));"), b = Tree::from_newick("((a,c),(b,d));");
        CHECK(quartet_distance(a, b) == 1);
        CHECK(triplet_distance(a, b) == 4);
        Tree c = Tree::from_newick("(a,b,c,d);");
        CaseCounts q = quartet_case_counts(a, c);
        CHECK(q.A == 0 && q.B == 0 && q.C == 1 && q.D == 0 && q.E == 0);
        CHECK(parameterized_quartet_distance(a, c, 0.5) == 0.5);
    }
    std::printf("test_tree: %s (%d failures)\n", fails ? "FAILED" : "ok", fails);
    return fails ? 1 : 0;
}
