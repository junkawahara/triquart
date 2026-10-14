# soda2013

C++17 implementation of Brodal et al., *Efficient Algorithms for
Computing the Triplet and Quartet Distance Between Trees of Arbitrary Degree*
(SODA 2013). The specification is the paper itself, kept locally as
`../soda13.pdf` (not committed); its figures are transcribed in
`../paper_figures.md`.

There are no external dependencies or build system. Working notes
(`status.md`, `try_and_error.md`) are local files and are not committed.

Compile and run the checks from the repository root:

```sh
mkdir -p build
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. \
    tree.cpp distance_statistics.cpp naive_distance.cpp tests/test_foundations.cpp \
    -o build/test_foundations
build/test_foundations
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. \
    tree.cpp component_tree.cpp hdt.cpp tests/test_hdt.cpp -o build/test_hdt
build/test_hdt
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. \
    tree.cpp component_tree.cpp hdt.cpp distance_statistics.cpp naive_distance.cpp \
    colored_hdt.cpp triplet_counters.cpp compression.cpp \
    tests/test_triplet_counters.cpp -o build/test_triplet_counters
build/test_triplet_counters
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. \
    tree.cpp component_tree.cpp hdt.cpp colored_hdt.cpp triplet_counters.cpp compression.cpp \
    tests/test_compression.cpp -o build/test_compression
build/test_compression
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. \
    tree.cpp component_tree.cpp hdt.cpp distance_statistics.cpp naive_distance.cpp \
    colored_hdt.cpp triplet_counters.cpp compression.cpp triplet_distance.cpp \
    tests/test_triplet_distance.cpp -o build/test_triplet_distance
build/test_triplet_distance
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. \
    tree.cpp component_tree.cpp hdt.cpp quartet_shared_counters.cpp \
    tests/test_quartet_shared_counters.cpp -o build/test_quartet_shared_counters
build/test_quartet_shared_counters
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. \
    tree.cpp component_tree.cpp hdt.cpp distance_statistics.cpp naive_distance.cpp \
    quartet_shared_counters.cpp quartet_contributions.cpp \
    quartet_disagreement_counters.cpp quartet_disagreement_contributions.cpp \
    tests/test_quartet_contributions.cpp -o build/test_quartet_contributions
build/test_quartet_contributions
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. \
    tree.cpp component_tree.cpp hdt.cpp quartet_shared_counters.cpp \
    quartet_disagreement_counters.cpp tests/test_quartet_disagreement_counters.cpp \
    -o build/test_quartet_disagreement_counters
build/test_quartet_disagreement_counters
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. \
    tree.cpp component_tree.cpp hdt.cpp colored_hdt.cpp compression.cpp \
    distance_statistics.cpp naive_distance.cpp quartet_shared_counters.cpp \
    quartet_contributions.cpp quartet_disagreement_counters.cpp \
    quartet_disagreement_contributions.cpp quartet_hdt.cpp \
    tests/test_quartet_hdt.cpp -o build/test_quartet_hdt
build/test_quartet_hdt
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -I. \
    tree.cpp component_tree.cpp hdt.cpp colored_hdt.cpp compression.cpp \
    distance_statistics.cpp naive_distance.cpp quartet_shared_counters.cpp \
    quartet_contributions.cpp quartet_disagreement_counters.cpp \
    quartet_disagreement_contributions.cpp quartet_hdt.cpp quartet_distance.cpp \
    tests/test_quartet_distance.cpp -o build/test_quartet_distance
build/test_quartet_distance
```

`Tree::from_newick` accepts rooted trees with unique, nonempty leaf labels and at
least two children at every internal node. A single labelled leaf is also valid.
Labels can be single-quoted; comments, internal labels and finite edge lengths
are ignored for topology. A terminal semicolon is required. Tree pairs must have
identical leaf-label sets. Invalid inputs throw `std::invalid_argument`.

Counts use checked `std::uint64_t` arithmetic; an unrepresentable counter or
intermediate value throws `std::overflow_error`. `DistanceCounts` describes the
five categories A–E of Figure 3 and computes `B + p(C + D)` for `0 <= p <= 1`.

The production triplet API uses the Section 4 HDT, Figure 10 counters and the
Section 5 `Count` traversal with `extract` / `contract`:

```cpp
#include "triplet_distance.hpp"

const auto a = soda2013::Tree::from_newick("((a,b),c);");
const auto b = soda2013::Tree::from_newick("(a,b,c);");
const auto counts = soda2013::triplet_counts(a, b);  // C = 1
const auto distance = counts.distance();           // 1
const auto weighted = counts.parameterized(0.5L);  // 0.5
```

Compile an application with `tree.cpp component_tree.cpp hdt.cpp
distance_statistics.cpp colored_hdt.cpp triplet_counters.cpp compression.cpp triplet_distance.cpp`.
The worst-case time is O(n log n), including O(n log n) label comparisons during
preprocessing. `Count` uses direct paired-leaf indices. `TripletDistanceStats`
exposes operation counts, including pointer updates at reconstruction.

The quartet API in `quartet_distance.hpp` supplies `quartet_counts`,
`quartet_distance`, and `parameterized_quartet_distance`. It uses the same
Section 5 traversal with the complete Figures 11–16 decorations and both tree
orientations. `QuartetDistanceStats` contains all nine A/B cases and work counts.
Compile with `tree.cpp component_tree.cpp hdt.cpp colored_hdt.cpp compression.cpp
distance_statistics.cpp quartet_shared_counters.cpp quartet_contributions.cpp
quartet_disagreement_counters.cpp quartet_disagreement_contributions.cpp
quartet_hdt.cpp quartet_distance.cpp`. Its worst-case time is O(d n log n).
Both bounds use the paper's unit-cost counting model; counts and component IDs
are fixed-width integers. Comparing string labels additionally costs the number
of characters examined.

The optional `ColoredTripletHdt` / `ColoredQuartetHdt` interfaces accept original
leaf IDs through a lazily built sorted index (O(m + k log k) construction and
O(log k) lookup for m component atoms and k active leaves). The production distance APIs use the
`recolor_atoms` / `extract_atoms` / `contract_atoms` variants and constant-time
paired indices, so this convenience index is outside their running-time bound.

`naive_distance.cpp` is a small-tree test oracle only. It must not be linked into
production as a replacement for the paper's HDT algorithm.
