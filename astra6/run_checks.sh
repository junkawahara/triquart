set -eu
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
