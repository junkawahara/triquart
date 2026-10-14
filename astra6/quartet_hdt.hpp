#pragma once

// soda13.pdf, Section 7 and Figures 11–16. A full quartet decoration carries
// all counters and both six-case contribution sums, updated through Section 5.
#include "quartet_disagreement_contributions.hpp"

namespace soda2013 {

struct ResolvedQuartets {
    QuartetCases agreeing;
    QuartetCases disagreeing;
    void accumulate(const ResolvedQuartets& other) {
        agreeing.accumulate(other.agreeing);
        disagreeing.accumulate(other.disagreeing);
    }
};

class ColoredQuartetHdt : public ColoredHdt {
public:
    explicit ColoredQuartetHdt(const Tree& tree, Color initial_color = 1);
    explicit ColoredQuartetHdt(ComponentTree tree, Color initial_color = 1);
    ResolvedQuartets query() const { return totals_[hdt_.root()]; }
    const QuartetCounters& counters(NodeId uid) const { return counters_.at(uid); }

private:
    std::vector<QuartetCounters> counters_;
    std::vector<ResolvedQuartets> totals_;
    void recompute(NodeId uid) override;
};

}  // namespace soda2013
