#pragma once

// soda13.pdf, Figures 10–12: only the counters needed by shared triplets.
// Each component stores entries only for positive colors present inside it.
// Sorted sparse vectors allow composition in O(1+x), without a dense d scan
// or a logarithmic dictionary lookup for each color (Section 7.1).
#include "colored_hdt.hpp"

namespace soda2013 {

struct TripletColorCounters {
    Color color = 0;
    Count leaves = 0;        // n_i
    Count same_inside = 0;   // n_(ii): LCA off the external path / super root
    Count different_at = 0;  // n_i•: different colors, branching at the path / super root
    Count i_up_other = 0;    // n_i↑•: i branches below the other positive color
};

struct SharedTriplets {
    Count resolved = 0;
    Count unresolved = 0;
    void accumulate(const SharedTriplets& other) {
        resolved = add(resolved, other.resolved);
        unresolved = add(unresolved, other.unresolved);
    }
};

struct TripletCounters {
    Count zero = 0;          // n_0 (kept explicitly for future quartet counters)
    Count positive = 0;      // n_•
    Count different_at = 0;  // n_•□ = 1/2 sum_i n_i•; unordered pairs
    std::vector<TripletColorCounters> colors;
    SharedTriplets shared;  // Sum of Figure 10 contributions in this HDT subtree
};

using TripletUpdateStats = ColorUpdateStats;

class ColoredTripletHdt : public ColoredHdt {
public:
    explicit ColoredTripletHdt(const Tree& tree, Color initial_color = 1);
    explicit ColoredTripletHdt(ComponentTree tree, Color initial_color = 1);
    SharedTriplets query() const { return counters_[hdt_.root()].shared; }
    const TripletCounters& counters(NodeId uid) const { return counters_.at(uid); }

private:
    std::vector<TripletCounters> counters_;
    void recompute(NodeId uid) override;
};

}  // namespace soda2013
