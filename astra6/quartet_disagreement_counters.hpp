#pragma once

// soda13.pdf, Figures 14–15. Figure 14 repeats 0(i•) and [0(i•)]; identical
// definitions share storage. This leaves 6 scalar, 33 indexed, and 2 pair types.
// i,j are distinct fixed positive colors; b,c are the bullet/square of Sec. 7.
#include "quartet_shared_counters.hpp"

namespace soda2013 {

enum class DisagreementScalar {
    b_p_b_c, p_b_p_b_c, any_b_p_b_c, b_up_b_c, p_b_c_up_b, b_up_b_up_c, size
};
enum class DisagreementIndexed {
    i_p_zero_i, b_p_zero_i, zero_p_i_b, b_p_i_c, i_p_i_b, b_p_i_b,
    p_i_p_i_b, p_i_p_zero_i, p_zero_p_i_b, p_b_p_zero_i,
    any_i_p_i_b, any_zero_p_i_b, any_i_p_zero_i, any_b_p_zero_i,
    zero_up_i_b, b_up_zero_i, i_up_i_b,
    p_i_b_up_zero, p_i_b_up_i, p_zero_i_up_i, p_zero_i_up_b,
    zero_up_p_i_b, i_up_p_i_b, b_up_p_zero_i, i_up_p_zero_i, b_up_p_i_c,
    b_up_i_up_zero, i_up_i_up_zero, zero_up_b_up_i, b_up_zero_up_i,
    i_up_zero_up_i, zero_up_i_up_b, i_up_i_up_b, size
};
inline constexpr std::size_t disagreement_scalar_count = static_cast<std::size_t>(DisagreementScalar::size);
inline constexpr std::size_t disagreement_indexed_count = static_cast<std::size_t>(DisagreementIndexed::size);
static_assert(disagreement_scalar_count == 6 && disagreement_indexed_count == 33);
inline constexpr std::array<std::string_view, disagreement_scalar_count> disagreement_scalar_subscripts{
    "b(bc)", "(b(bc))", "[b(bc)]", "b^bc", "(bc)^b", "b^b^c"
};
inline constexpr std::array<std::string_view, disagreement_indexed_count> disagreement_indexed_subscripts{
    "i(0i)", "b(0i)", "0(ib)", "b(ic)", "i(ib)", "b(ib)",
    "(i(ib))", "(i(0i))", "(0(ib))", "(b(0i))",
    "[i(ib)]", "[0(ib)]", "[i(0i)]", "[b(0i)]",
    "0^ib", "b^0i", "i^ib", "(ib)^0", "(ib)^i", "(0i)^i", "(0i)^b",
    "0^(ib)", "i^(ib)", "b^(0i)", "i^(0i)", "b^(ic)",
    "b^i^0", "i^i^0", "0^b^i", "b^0^i", "i^0^i", "0^i^b", "i^i^b"
};

struct DisagreementColor {
    std::array<Count, disagreement_indexed_count> values{};
    // Figure 16 refers to these additional C-only rows, omitted in Figure 14.
    Count b_up_p_i_b = 0;  // n_•↑(i•)
    Count b_up_i_b = 0;    // n_•↑i•
    Count b_up_i_c = 0;    // n_•↑i□
    Count i_up_zero_i = 0; // n_i↑0i
    Count& operator[](DisagreementIndexed counter) { return values[static_cast<std::size_t>(counter)]; }
    Count operator[](DisagreementIndexed counter) const { return values[static_cast<std::size_t>(counter)]; }
};

struct DisagreementPair {
    Count at = 0;      // n_ij, i != j, unordered fixed-color pair at boundary.
    Count inside = 0;  // n_(ij), LCA off the external path / super root.
};

struct QuartetCounters {
    QuartetSharedCounters shared;
    std::array<Count, disagreement_scalar_count> values{};
    std::vector<DisagreementColor> colors;  // Same sorted palette as shared.colors.
    // Symmetric pair counters occupy the upper triangle of the PRESENT-color
    // matrix. There is no allocation over globally possible/absent colors.
    std::vector<DisagreementPair> pairs;
    Count b_up_p_b_c = 0;  // Figure 16: n_•↑(•□) = sum_i n_i↑(i•).
    Count b_up_c_up_b = 0; // Figure 16: n_•↑□↑•.
    Count& operator[](DisagreementScalar counter) { return values[static_cast<std::size_t>(counter)]; }
    Count operator[](DisagreementScalar counter) const { return values[static_cast<std::size_t>(counter)]; }
    const DisagreementPair& pair(std::size_t i, std::size_t j) const;
    DisagreementPair& pair(std::size_t i, std::size_t j);
};

QuartetCounters quartet_leaf(Color color, Count zero_multiplicity = 1);
QuartetCounters quartet_cc(const QuartetCounters& lower, const QuartetCounters& upper);
QuartetCounters quartet_gg(const QuartetCounters& first, const QuartetCounters& second);
QuartetCounters quartet_ig(const QuartetCounters& group);
void quartet_c_to_g(QuartetCounters& counters);

}  // namespace soda2013
