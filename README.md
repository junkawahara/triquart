# triquart

Exact triplet and quartet distances between two phylogenetic trees of
arbitrary degree, in C++17 with no external dependencies.

The code implements the algorithms of

> Gerth Stølting Brodal, Rolf Fagerberg and Christian N. S. Pedersen.
> *Efficient Algorithms for Computing the Triplet and Quartet Distance
> Between Trees of Arbitrary Degree.* Proceedings of the 24th Annual
> ACM-SIAM Symposium on Discrete Algorithms (SODA 2013), pp. 1814–1832.
> <https://doi.org/10.1137/1.9781611973105.130>

which compute the triplet distance between two rooted trees in
O(n log n) time and the quartet distance between two unrooted trees in
O(d n log n) time, where n is the number of leaves and d is the maximum
degree of the input trees. Both distances are exact integer counts; a
brute-force oracle is included only for testing.

## What is computed

Both trees must have the same set of n leaf labels. Every subset of three
leaves (a *triplet*) or four leaves (a *quartet*) induces a topology in
each tree, which is either *resolved* (binary) or *unresolved* (a star).
Comparing the two induced topologies puts each triplet or quartet into one
of five cells (Figure 3 of the paper):

| Cell | Meaning |
|------|---------|
| A | resolved in both trees, same topology |
| B | resolved in both trees, different topologies |
| C | resolved in the first tree only |
| D | resolved in the second tree only |
| E | unresolved in both trees |

The program reports all five counts together with

- the **distance** `B + C + D`, the number of triplets or quartets whose
  induced topologies differ, and
- the **parameterized distance** `B + p (C + D)` for a user-chosen
  `0 <= p <= 1` (Bansal et al.), which lets a resolved-versus-unresolved
  disagreement count less than a resolved-versus-resolved one.

## Repository layout

The repository contains two independent, self-contained implementations
of the same paper. They share no code and can be built separately.

| Directory | Contents |
|-----------|----------|
| `fable5_1/` | Library plus a command-line tool `soda13`. Built with the included `Makefile`. Totals are 128-bit integers where the compiler provides `__int128`, so quartet counts remain exact beyond n ≈ 123 000. Namespace `soda13`. |
| `astra6/` | Library with a header-level API and detailed per-run statistics (counter work, number of contractions, all nine resolved/unresolved cases). Compiled directly with `g++`; `astra6/run_checks.sh` builds and runs every test, and `astra6/README.md` documents the API in detail. Namespace `soda2013`. |

Both directories carry their own tests, including exhaustive comparison
against a naive enumeration oracle on all small trees and randomized
comparison on larger ones.

## Requirements

- A C++17 compiler (`g++` 9 or later, or a recent `clang++`).
- `make` for `fable5_1/`.

Nothing else is required.

## Quick start: command-line tool

```sh
cd fable5_1
make            # builds build/soda13 and the test binaries
make test       # runs all unit tests and the CLI regression test
```

Usage:

```
build/soda13 [--triplet|--quartet|--both] [-p P] [--naive] T1.nwk T2.nwk
```

| Option | Effect |
|--------|--------|
| `--triplet` | compute the triplet distance only |
| `--quartet` | compute the quartet distance only |
| `--both` | compute both (default) |
| `-p P` | parameter for `B + p (C + D)`; default `1` |
| `--naive` | use the brute-force oracle instead of the HDT algorithm (small trees only) |

Example. With `t1.nwk` containing `((a,b),(c,d),e);` and `t2.nwk`
containing `((a,c),(b,d),e);`:

```sh
$ build/soda13 -p 0.5 t1.nwk t2.nwk
triplets: A=0 B=4 C=2 D=2 E=2 distance=8 parameterized(p=0.5)=6
quartets: A=0 B=5 C=0 D=0 E=0 distance=5 parameterized(p=0.5)=5
```

Exit status is `0` on success, `1` for invalid input (for example,
differing leaf sets) and `2` for a usage or file error. Error messages go
to standard error.

## Input format

Trees are read in Newick format, one tree per file, terminated by a
semicolon.

- Internal nodes may have any number of children. Nodes with a single
  child are contracted, since they do not affect any induced topology.
- Leaf labels must be unique and non-empty. Labels may be single-quoted.
- Branch lengths, internal node labels and comments are accepted and
  ignored.
- Both trees must have exactly the same set of leaf labels.
- For the quartet distance the trees are treated as unrooted; the root in
  the Newick string is irrelevant. For the triplet distance the trees are
  rooted at the Newick root.

## Using the libraries from C++

### `fable5_1`

```cpp
#include "distance.hpp"   // in fable5_1/

using namespace soda13;

Tree t1 = Tree::from_newick("((a,b),(c,d),e);");
Tree t2 = Tree::from_newick("((a,c),(b,d),e);");

CaseCounts q = quartet_case_counts(t1, t2);   // q.A ... q.E
count_t    d = quartet_distance(t1, t2);      // B + C + D, 128-bit
double     w = parameterized_quartet_distance(t1, t2, 0.5);

CaseCounts t = triplet_case_counts(t1, t2);
```

Link against `tree.cpp single_tree.cpp naive_distance.cpp hdt.cpp count.cpp`
(the `LIB_SRCS` of the `Makefile`). Invalid trees or mismatched leaf sets
throw `std::invalid_argument`.

### `astra6`

```cpp
#include "triplet_distance.hpp"   // in astra6/
#include "quartet_distance.hpp"

using namespace soda2013;

const auto t1 = Tree::from_newick("((a,b),(c,d),e);");
const auto t2 = Tree::from_newick("((a,c),(b,d),e);");

DistanceCounts tc = triplet_counts(t1, t2);         // tc.A ... tc.E
Count          td = tc.distance();                  // B + C + D
long double    tw = tc.parameterized(0.5L);         // B + p (C + D)

QuartetDistanceStats stats;
DistanceCounts qc = quartet_counts(t1, t2, &stats); // stats: work counts
```

Compile an application together with

```
tree.cpp component_tree.cpp hdt.cpp colored_hdt.cpp compression.cpp
distance_statistics.cpp triplet_counters.cpp triplet_distance.cpp
quartet_shared_counters.cpp quartet_contributions.cpp
quartet_disagreement_counters.cpp quartet_disagreement_contributions.cpp
quartet_hdt.cpp quartet_distance.cpp
```

Counts use checked 64-bit arithmetic and throw `std::overflow_error` if a
value is not representable. `naive_distance.cpp` is a test oracle and is
not needed by the production API.

## Running the tests

```sh
# fable5_1: unit tests, oracle comparison, CLI regression
cd fable5_1 && make test

# astra6: compiles and runs every test with g++ (same commands as in
# astra6/README.md)
cd astra6 && sh run_checks.sh
```

The `astra6` checks take a few minutes: the quartet tests compare against
the oracle exhaustively on all small trees and on thousands of random
pairs, and then measure the counter work on trees with up to 10 000
leaves.

## Benchmark

`fable5_1/tests/bench.cpp` measures the running time on random trees of
growing size and prints time divided by `n log n`, so the O(n log n) and
O(d n log n) bounds can be checked empirically:

```sh
cd fable5_1
g++ -std=c++17 -O2 -I. tests/bench.cpp tree.cpp single_tree.cpp \
    naive_distance.cpp hdt.cpp count.cpp -o build/bench
build/bench tq 100000 b     # mode: t | q | tq; max n; shape: b | m | s | S | c
```

Shapes are `b` (near-binary), `m` (mixed degrees), `s` (high degree),
`S` (star of about sqrt(n) binary subtrees) and `c` (caterpillar).

## Notes on the implementation

Both implementations follow the paper closely: a hierarchical
decomposition tree (HDT) of the second tree, the counters of Figures
10–16 stored in every HDT node, and the `Count` traversal of the first
tree with coloring, `extract` and `contract` so that the compressed HDT
stays proportional to the current subtree. The transcription of the
paper's counter tables that the code was checked against, including a
few corrected typos, is kept as a local working note and is not part of
the repository.
