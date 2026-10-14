#!/bin/sh
# Build fable5_1 and astra6 and compare their triplet / quartet counts on
# edge cases, random trees, large trees and the naive oracle.
#
#   sh crosscheck/run.sh          full run (a few minutes)
#   sh crosscheck/run.sh --quick  reduced run
#
# Exit status is 0 only if every case agrees.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
root=$(dirname "$here")
mkdir -p "$here/build"

make -C "$root/fable5_1" build/soda13 >/dev/null

A="$root/astra6"
g++ -std=c++17 -O2 -I"$A" "$here/astra_cli.cpp" \
    "$A/tree.cpp" "$A/component_tree.cpp" "$A/hdt.cpp" "$A/colored_hdt.cpp" \
    "$A/compression.cpp" "$A/distance_statistics.cpp" "$A/triplet_counters.cpp" \
    "$A/triplet_distance.cpp" "$A/quartet_shared_counters.cpp" \
    "$A/quartet_contributions.cpp" "$A/quartet_disagreement_counters.cpp" \
    "$A/quartet_disagreement_contributions.cpp" "$A/quartet_hdt.cpp" \
    "$A/quartet_distance.cpp" -o "$here/build/astra_cli"

# Deeply nested caterpillar trees need a large stack in the C++ parsers.
ulimit -s unlimited 2>/dev/null || ulimit -s "$(ulimit -H -s)" 2>/dev/null || true

python3 "$here/crosscheck.py" "$root/fable5_1/build/soda13" "$here/build/astra_cli" "$@"
