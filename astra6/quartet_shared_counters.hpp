#pragma once

// soda13.pdf, Section 7.1 and Figures 11–12: ALL 53 shared-counter types.
// ASCII subscripts: b=bullet, c=square, p=parentheses, any=square brackets,
// up=bottom-to-top order along the C external path. If i appears, b and c
// exclude i; b and c also differ from one another. Repeated symbols agree.
// The actual figure has 14 scalar + 39 indexed types (the prose says 15+38).
#include "colored_hdt.hpp"

#include <array>
#include <string_view>

namespace soda2013 {

enum class SharedScalar {
    zero, b, zero_b, bb, bc,
    p_zero_b, p_bb, p_bc, c_p_bb,
    zero_up_b, b_up_zero, zero_up_bb, b_up_p_cc, any_bc,
    size
};

enum class SharedIndexed {
    i, zero_i, ii, i_b,
    p_zero_i, p_ii, p_i_b,
    zero_p_ii, b_p_ii, i_p_zero_b, i_p_bb, i_p_bc,
    p_zero_p_ii, p_b_p_ii, p_i_p_zero_b,
    zero_up_i, i_up_zero, i_up_i, i_up_b, b_up_i,
    i_up_zero_b, b_up_ii, i_up_bb, i_up_bc,
    zero_up_i_up_i, b_up_i_up_i, i_up_b_up_zero, i_up_zero_up_b,
    zero_up_p_ii, i_up_p_zero_b, i_up_p_bc, b_up_p_ii,
    p_ii_up_zero, p_ii_up_b, p_bb_up_i,
    any_i_b, any_zero_p_ii, any_b_p_ii, any_i_p_zero_b,
    size
};

inline constexpr std::size_t shared_scalar_count = static_cast<std::size_t>(SharedScalar::size);
inline constexpr std::size_t shared_indexed_count = static_cast<std::size_t>(SharedIndexed::size);
static_assert(shared_scalar_count == 14 && shared_indexed_count == 39);

// Literal Figure 11 subscripts; ^ is an upward arrow. Also used by the test
// oracle to interpret leaf color/topology constraints, independently of Fig. 12.
inline constexpr std::array<std::string_view, shared_scalar_count> shared_scalar_subscripts{
    "0", "b", "0b", "bb", "bc", "(0b)", "(bb)", "(bc)", "c(bb)",
    "0^b", "b^0", "0^bb", "b^(cc)", "[bc]"
};
inline constexpr std::array<std::string_view, shared_indexed_count> shared_indexed_subscripts{
    "i", "0i", "ii", "ib", "(0i)", "(ii)", "(ib)",
    "0(ii)", "b(ii)", "i(0b)", "i(bb)", "i(bc)",
    "(0(ii))", "(b(ii))", "(i(0b))",
    "0^i", "i^0", "i^i", "i^b", "b^i", "i^0b", "b^ii", "i^bb", "i^bc",
    "0^i^i", "b^i^i", "i^b^0", "i^0^b",
    "0^(ii)", "i^(0b)", "i^(bc)", "b^(ii)", "(ii)^0", "(ii)^b", "(bb)^i",
    "[ib]", "[0(ii)]", "[b(ii)]", "[i(0b)]"
};

struct QuartetSharedColor {
    Color color = 0;
    std::array<Count, shared_indexed_count> values{};
    // Used explicitly by Figure 13, but absent from the Figures 11–12 list.
    // Their meanings follow the same Section 7 notation; see Trial 7.
    Count i_up_p_bb = 0;       // n_i↑(••)
    Count p_zero_b_up_i = 0;   // n_(0•)↑i
    Count zero_up_ii = 0;      // n_0↑ii
    Count& operator[](SharedIndexed counter) { return values[static_cast<std::size_t>(counter)]; }
    Count operator[](SharedIndexed counter) const { return values[static_cast<std::size_t>(counter)]; }
};

struct QuartetSharedCounters {
    std::array<Count, shared_scalar_count> values{};
    std::vector<QuartetSharedColor> colors;  // Sorted; positive colors present only.
    Count c_up_bb = 0;  // Figure 13: n_□↑•• = sum_i n_•↑ii.
    Count& operator[](SharedScalar counter) { return values[static_cast<std::size_t>(counter)]; }
    Count operator[](SharedScalar counter) const { return values[static_cast<std::size_t>(counter)]; }
};

// Figure 12's five actions. For CC, first is the LOWER C1 and second the UPPER C2.
// The zero leaf multiplicity also initializes frozen all-zero C/G atoms from
// Section 5. A positive leaf always has multiplicity one.
QuartetSharedCounters quartet_shared_leaf(Color color, Count zero_multiplicity = 1);
QuartetSharedCounters quartet_shared_cc(const QuartetSharedCounters& lower, const QuartetSharedCounters& upper);
QuartetSharedCounters quartet_shared_gg(const QuartetSharedCounters& first, const QuartetSharedCounters& second);
QuartetSharedCounters quartet_shared_ig(const QuartetSharedCounters& group);
void quartet_shared_c_to_g(QuartetSharedCounters& counters);

}  // namespace soda2013
