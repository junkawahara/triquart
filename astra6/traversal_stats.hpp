#pragma once

// soda13.pdf, Section 5: work charged to the same Count/extract/contract
// traversal for triplets and quartets (one tree orientation).
#include "colored_hdt.hpp"

namespace soda2013 {

struct CountTraversalStats {
    std::size_t queries = 0;
    std::size_t hdt_builds = 0;
    std::size_t hdt_input_components = 0;
    std::size_t hdt_round_visits = 0;
    std::size_t heavy_contractions = 0;
    std::size_t recolored_leaves = 0;
    std::size_t paired_leaf_updates = 0;
    std::size_t maximum_components_per_active_leaf = 0;
    ColorUpdateStats updates;
    CompressionStats compression;
};

}  // namespace soda2013
