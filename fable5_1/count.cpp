// count.cpp — instantiates the Sec. 5 engine for triplets and quartets and
// assembles the five cells of Fig. 3 (see distance.hpp).
#include "count.hpp"

#include "distance.hpp"
#include "single_tree.hpp"

namespace soda13 {

namespace {

void normalise(Tree& a, Tree& b) {
    assign_leaf_ids(a, b);
    a.contract_unary();
    b.contract_unary();
}

}  // namespace

CaseCounts triplet_case_counts(const Tree& t1in, const Tree& t2in) {
    Tree t1 = t1in, t2 = t2in;
    normalise(t1, t2);
    // Sec. 2: A and E from the HDT algorithm, the rest from the Sec. 3 sums.
    Engine<TripletFields> eng(t1, t2, /*asym_only=*/false);
    const Contrib c = eng.run();
    const SingleTreeCounts s1 = count_single_tree(t1), s2 = count_single_tree(t2);
    CaseCounts r;
    r.A = c.a;
    r.E = c.b;
    r.C = s2.unresolved_triplets - r.E;  // C + E = unresolved in T2
    r.D = s1.unresolved_triplets - r.E;  // D + E = unresolved in T1
    r.B = s1.resolved_triplets - r.A - r.C;
    return r;
}

CaseCounts quartet_case_counts(const Tree& t1in, const Tree& t2in) {
    Tree t1 = t1in, t2 = t2in;
    normalise(t1, t2);
    // Sec. 6 / Fig. 8–9: the nine (T1-type, T2-type) cases are covered by a
    // direct run (αα, βα, ββ, γα, γβ, γγ) plus a run with the trees swapped
    // that only counts the asymmetric cases (βα, γα, γβ = αβ, αγ, βγ).
    Engine<QuartetFields> direct(t1, t2, /*asym_only=*/false);
    Contrib c = direct.run();
    Engine<QuartetFields> swapped(t2, t1, /*asym_only=*/true);
    c += swapped.run();
    const SingleTreeCounts s1 = count_single_tree(t1), s2 = count_single_tree(t2);
    CaseCounts r;
    r.A = c.a;
    r.B = c.b;
    r.C = s1.resolved_quartets - r.A - r.B;  // A + B + C = resolved in T1
    r.D = s2.resolved_quartets - r.A - r.B;  // A + B + D = resolved in T2
    r.E = s1.unresolved_quartets - r.D;      // D + E = unresolved in T1
    return r;
}

count_t triplet_distance(const Tree& t1, const Tree& t2) { return triplet_case_counts(t1, t2).distance(); }
count_t quartet_distance(const Tree& t1, const Tree& t2) { return quartet_case_counts(t1, t2).distance(); }
double parameterized_triplet_distance(const Tree& t1, const Tree& t2, double p) {
    return triplet_case_counts(t1, t2).parameterized(p);
}
double parameterized_quartet_distance(const Tree& t1, const Tree& t2, double p) {
    return quartet_case_counts(t1, t2).parameterized(p);
}

}  // namespace soda13
